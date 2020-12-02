#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <algorithm>
#include <vector>
#include <iostream>

#include <jack/jack.h>

// Include FFTW header
#include <complex> //needs to be included before fftw3.h for compatibility
#include <fftw3.h>

// Eigen include
#include <Eigen/Eigen>

// Include unordered map
#include<unordered_map>

std::complex<double> *i_fft, *i_time, *o_fft, *o_time;
fftw_plan i_forward, o_inverse;
int ports_number = 3;
int OUT_SYSTEM_PORTS = 2;
jack_port_t **input_port;
jack_port_t **output_port;
jack_default_audio_sample_t **in, **out;

jack_client_t *client;

double sample_rate;
double *buffer, *fft1, *fft2, *hann;

unsigned int buffer_size;
unsigned int fft_size;

unsigned int in_i, fft1_i, fft2_i, fft1_offset, fft2_offset, out2_offset;

double *freqs;
std::complex<double> *delay_vector;
double delay_seconds;
const double MIC_DIST = 0.18;
const double SOUND_SPEED = 343.0;
const std::complex<double> COMP_I(0.0, 1.0);

void filter(double *data) {
	int i;

	// Hann
	for (i = 0; i < fft_size; ++i) {
		i_time[i] = data[i] * hann[i];
	}

	// FFT
	fftw_execute(i_forward);

	// Filter
	for (i = 0; i < fft_size; ++i) {
		o_fft[i] = i_fft[i] * delay_vector[i];
	}

	// ifft
	fftw_execute(o_inverse);

	// return
	for(i = 0; i < fft_size; i++){
		data[i] = std::real(o_time[i])/(double)fft_size; //fftw3 requiere normalizar su salida real de esta manera
	}
}

/**
 * The process callback for this JACK application is called in a
 * special realtime thread once for each audio cycle.
 *
 * This client does nothing more than copy data from its input
 * port to its output port. It will exit when stopped by 
 * the user (e.g. using Ctrl-C on a unix-ish operating system)
 */
int jack_callback (jack_nframes_t nframes, void *arg){
	int i;
	
	int testing_port_id = 0;

	for (int port_i = 0; port_i < ports_number; ++port_i) {
        in[port_i] = (jack_default_audio_sample_t *)jack_port_get_buffer (input_port[port_i], nframes);
	    	out[port_i] = (jack_default_audio_sample_t *)jack_port_get_buffer (output_port[port_i], nframes);
        for (i = 0; i < nframes; ++i) {
            out[port_i][i] = in[port_i][i];
        }
    }
	
	// copy in to las nframes of buffer
	for (i = 0; i < nframes; ++i, ++in_i) {
		if (in_i >= buffer_size) in_i = 0;
		buffer[in_i] = in[testing_port_id][i];
	}

	// Fill fft1 y fft2
	unsigned int fft1_i_local = fft1_i; 
	unsigned int fft2_i_local = fft2_i; 
	for (i = 0; i < fft_size; ++i, ++fft1_i_local, ++fft2_i_local) {
		if (fft1_i_local >= buffer_size) fft1_i_local = 0;
		fft1[i] = buffer[fft1_i_local];

		if (fft2_i_local >= buffer_size) fft2_i_local = 0;
		fft2[i] = buffer[fft2_i_local];
	}

	fft1_i += nframes;
	if (fft1_i >= buffer_size) fft1_i = 0;
	
	fft2_i += nframes;
	if (fft2_i >= buffer_size) fft2_i = 0;

	filter(fft1);
	filter(fft2);
	
	for (i = 0; i < nframes; ++i) {
		out[0][i] = fft1[fft1_offset + i] + fft2[fft2_offset + i];
	}

	// unphased output
	int out_i = in_i - out2_offset;
	if (out_i < 0) out_i += buffer_size;
	for (i = 0; i < nframes; ++i, ++out_i) {
		if (out_i >= buffer_size) out_i = 0;
		out[1][i] = buffer[out_i];
	}


	
	return 0;
}


/**
 * JACK calls this shutdown_callback if the server ever shuts down or
 * decides to disconnect the client.
 */
void jack_shutdown (void *arg){
	exit (1);
}


int main (int argc, char *argv[]) {
	const char *client_name = "das_simple";
	jack_options_t options = JackNoStartServer;
	jack_status_t status;
	int i;

	delay_seconds = 0.020;
	
	/* open a client connection to the JACK server */
	client = jack_client_open (client_name, options, &status);
	if (client == NULL){
		/* if connection failed, say why */
		printf ("jack_client_open() failed, status = 0x%2.0x\n", status);
		if (status & JackServerFailed) {
			printf ("Unable to connect to JACK server.\n");
		}
		exit (1);
	}
	
	/* if connection was successful, check if the name we proposed is not in use */
	if (status & JackNameNotUnique){
		client_name = jack_get_client_name(client);
		printf ("Warning: other agent with our name is running, `%s' has been assigned to us.\n", client_name);
	}
	
	/* tell the JACK server to call 'jack_callback()' whenever there is work to be done. */
	jack_set_process_callback (client, jack_callback, 0);
	
	
	/* tell the JACK server to call 'jack_shutdown()' if it ever shuts down,
	   either entirely, or if it just decides to stop calling us. */
	jack_on_shutdown (client, jack_shutdown, 0);
	
	
	/* display the current sample rate. */
	printf ("Sample rate: %d\n", jack_get_sample_rate (client));
	printf ("Window size: %d\n", jack_get_buffer_size (client));
	
	sample_rate = (double)jack_get_sample_rate(client);
	int nframes = jack_get_buffer_size (client);
	
	
	// Preparing overlap and add sizes
	buffer_size = 6 * nframes;
	fft_size = 4 * nframes;
	in_i = 5 * nframes;
	fft1_i = 0;
	fft2_i = 2 * nframes;
	fft1_offset = (int)(2.5 * (double)nframes);
	fft2_offset = (int)(0.5 * (double)nframes);
	out2_offset = (int)(3.5 * (double)nframes);

	// printf("%d %d %d\n", fft1_offset, fft2_offset, out2_offset);

	
	//preparing FFTW3 buffers
	i_fft = (std::complex<double> *) fftw_malloc(sizeof(std::complex<double>) * fft_size);
	i_time = (std::complex<double> *) fftw_malloc(sizeof(std::complex<double>) * fft_size);
	o_fft = (std::complex<double> *) fftw_malloc(sizeof(std::complex<double>) * fft_size);
	o_time = (std::complex<double> *) fftw_malloc(sizeof(std::complex<double>) * fft_size);
	
	i_forward = fftw_plan_dft_1d(fft_size, reinterpret_cast<fftw_complex*>(i_time), reinterpret_cast<fftw_complex*>(i_fft) , FFTW_FORWARD, FFTW_MEASURE);
	o_inverse = fftw_plan_dft_1d(fft_size, reinterpret_cast<fftw_complex*>(o_fft) , reinterpret_cast<fftw_complex*>(o_time), FFTW_BACKWARD, FFTW_MEASURE);
	
	
	buffer = (double *) calloc(buffer_size, sizeof(double));
	fft1 = (double *) calloc(fft_size, sizeof(double));
	fft2 = (double *) calloc(fft_size, sizeof(double));
	hann = (double *) calloc(fft_size, sizeof(double));
	delay_vector = (std::complex<double> *) malloc(fft_size * sizeof(std::complex<double>));
	freqs = (double *) calloc(fft_size, sizeof(double));

	// Calculate hann window
	// std::cout << "[" << std::endl;
	for (i = 0; i < fft_size;++i) {
		hann[i] = 0.5 * (1.0 - cos( 2.0 * M_PI * ((double)i)/((double)(fft_size - 1)) ));
		// std::cout << hann[i] << ",";
	}
	// std::cout << "]" << std::endl;
	// vector for frequency indexes
	freqs[0] = 0.0;
	freqs[fft_size/2] = sample_rate/2.00;
	int lim = fft_size/2;
	for (i = 1; i < lim; ++i) {
		freqs[i] = ((double)i)*(sample_rate/(double)(fft_size));
		freqs[fft_size-i] = -freqs[i];
	}

	// Exponencials for delay
	for (i = 0; i < fft_size; ++i) {
		std::complex<double> freq_res_complex = std::complex<double>(freqs[i],0.0);
		std::complex<double> delay_coefficient = (-MIC_DIST/SOUND_SPEED)*cos((-90 - (-30))*M_PI/180.0);
		std::cout << delay_coefficient << "  ++++" << std::endl;
		delay_vector[i] = std::exp(-COMP_I*2.0*M_PI*freq_res_complex*delay_coefficient);
	}

	input_port = (jack_port_t**) malloc(ports_number * sizeof (jack_port_t*));
  output_port = (jack_port_t**) malloc(ports_number * sizeof (jack_port_t*));
	in = (jack_default_audio_sample_t**) malloc(ports_number * sizeof (jack_default_audio_sample_t*));
  out = (jack_default_audio_sample_t**) malloc(ports_number * sizeof (jack_default_audio_sample_t*));
	char port_name[10];

	for (i = 0; i < ports_number; ++i) {
		
  	/* create the agent input ports */
    sprintf(port_name, "input_%d", (i+1));
    input_port[i] = jack_port_register (client, port_name, JACK_DEFAULT_AUDIO_TYPE,JackPortIsInput, 0);

    /* create the agent output port */
    sprintf(port_name, "output_%d", (i+1));
	  output_port[i] = jack_port_register (client, port_name, JACK_DEFAULT_AUDIO_TYPE,JackPortIsOutput, 0);    

  }

	// /* check that both ports were created succesfully */
	// if ((input_port == NULL) || (output_port == NULL) || (output_port2 == NULL)) {
	// 	printf("Could not create agent ports. Have we reached the maximum amount of JACK agent ports?\n");
	// 	exit (1);
	// }
	
	/* Tell the JACK server that we are ready to roll.
	   Our jack_callback() callback will start running now. */
	if (jack_activate (client)) {
		printf ("Cannot activate client.");
		exit (1);
	}
	
	printf ("Agent activated.\n");
	
	/* Connect the ports.  You can't do this before the client is
	 * activated, because we can't make connections to clients
	 * that aren't running.  Note the confusing (but necessary)
	 * orientation of the driver backend ports: playback ports are
	 * "input" to the backend, and capture ports are "output" from
	 * it.
	 */
	printf ("Connecting ports... ");
	 
	/* Assign our input port to a server output port*/
	// Find possible output server port names
	const char **serverports_names;
	serverports_names = jack_get_ports (client, NULL, NULL, JackPortIsPhysical|JackPortIsOutput);
	if (serverports_names == NULL) {
		printf("No available physical capture (server output) ports.\n");
		exit (1);
	}
	// // Connect the first available to our input port
	// if (jack_connect (client, serverports_names[0], jack_port_name (input_port))) {
	// 	printf("Cannot connect input port.\n");
	// 	exit (1);
	// }
	// free serverports_names variable for reuse in next part of the code
	free (serverports_names);
	
	
	/* Assign our output port to a server input port*/
	// Find possible input server port names
	serverports_names = jack_get_ports (client, NULL, NULL, JackPortIsPhysical|JackPortIsInput);
	if (serverports_names == NULL) {
		printf("No available physical playback (server input) ports.\n");
		exit (1);
	}
	for (i = 0; i < OUT_SYSTEM_PORTS; ++i) {

		// Connect the first available to our output port
		if (jack_connect (client, jack_port_name (output_port[i]), serverports_names[i])) {
			printf ("Cannot connect output ports.\n");
			exit (1);
		}

	}
	// free serverports_names variable, we're not going to use it again
	free (serverports_names);
	
	
	printf ("done.\n");
	/* keep running until stopped by the user */
	sleep (-1);
	
	
	/* this is never reached but if the program
	   had some other way to exit besides being killed,
	   they would be important to call.
	*/
	jack_client_close (client);
	exit (0);
}
