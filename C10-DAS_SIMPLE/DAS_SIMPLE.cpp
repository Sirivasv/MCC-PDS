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
// Include deque
#include <deque>

// Overlapping Split Buffer
template <class T>
class OSB {
		
		int total_size;
		std::deque<T> inner_buffer;

  public:
    
		OSB() {
		
			total_size = 0;

			inner_buffer = std::deque<T>();

		}

		OSB(int total_size_) {
		
			total_size = total_size_;

			inner_buffer = std::deque<T>(total_size);

		}

    void insert_element(T new_element_) {
			inner_buffer.push_back(new_element_);
			inner_buffer.pop_front();	
		}

		void legacy_fill(T *data, int fill_size, int start_position) {
			int lim_iter = start_position + fill_size;
			for(int b_i = start_position, i = 0; b_i < lim_iter; ++b_i, ++i){
				data[i] = inner_buffer[b_i];
			}
		}

};

std::complex<double> *i_fft, *i_time, *o_fft, *o_time;
fftw_plan i_forward, o_inverse;
const int ports_number = 3;
int OUT_SYSTEM_PORTS = 2;

OSB<double> osb_array[ports_number];

jack_port_t **input_port;
jack_port_t **output_port;
jack_default_audio_sample_t **in, **out;

jack_client_t *client;

double sample_rate;
double *hann;
double ***filtered_ffts;
std::complex<double> ***delayed_ffts;
std::complex<double> **delay_and_sum_fft;

double *doa_value;

double **delay_and_sum_output_ifft;

unsigned int buffer_size;
unsigned int fft_size;

unsigned int in_i, fft1_i, fft2_i, fft1_offset, fft2_offset, undelayed_offset;

double *freqs;
std::complex<double> *delay_vector;
double delay_seconds;
const double MIC_DIST = 0.18;
const double SOUND_SPEED = 343.0;
const std::complex<double> COMP_I(0.0, 1.0);

void delay_fft(double *data_in, std::complex<double> *data_out) {
	
	// Hann
	for (int i = 0; i < fft_size; ++i) {
		i_time[i] = data_in[i] * hann[i];
	}

	// FFT
	fftw_execute(i_forward);

	// Filter
	for (int i = 0; i < fft_size; ++i) {
		data_out[i] = i_fft[i] * delay_vector[i];
	}

}

void get_inverse_fft(std::complex<double> *data_in, double *data_out) {
	
	for (int i = 0; i < fft_size; ++i) {
	 o_fft[i] = data_in[i];
	}

	// ifft
	fftw_execute(o_inverse);

	// return
	for(int i = 0; i < fft_size; i++){
		data_out[i] = std::real(o_time[i])/(double)fft_size; //fftw3 requiere normalizar su salida real de esta manera
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
		
	
	for (int port_i = 0; port_i < ports_number; ++port_i) {
		in[port_i] = (jack_default_audio_sample_t *)jack_port_get_buffer (input_port[port_i], nframes);
		out[port_i] = (jack_default_audio_sample_t *)jack_port_get_buffer (output_port[port_i], nframes);
  }

	// Reset to zero
	for (int i = 0; i < fft_size; ++i) {
		delay_and_sum_fft[0][i] = 0;
		delay_and_sum_fft[1][i] = 0;
	}

	for (int current_mic = 0; current_mic < ports_number; current_mic++) {

		// copy in to las nframes of buffer
		for (int i = 0; i < nframes; ++i) {
			osb_array[current_mic].insert_element(in[current_mic][i]);
		}

		// Fill fft1 y fft2
		osb_array[current_mic].legacy_fill(filtered_ffts[current_mic][0], fft_size, fft1_i);
		osb_array[current_mic].legacy_fill(filtered_ffts[current_mic][1], fft_size, fft2_i);

		// Apply Filter to a given angle
		// Direction angle
		std::complex<double> direction_angle_coefficient = (-MIC_DIST/SOUND_SPEED)*cos((doa_value[current_mic])*M_PI/180.0);
		for (int i = 0; i < fft_size; ++i) {
			std::complex<double> freq_res_complex = std::complex<double>(freqs[i],0.0);
			delay_vector[i] = std::exp(-COMP_I*2.0*M_PI*freq_res_complex*direction_angle_coefficient);
		}

		delay_fft(filtered_ffts[current_mic][0], delayed_ffts[current_mic][0]);
		delay_fft(filtered_ffts[current_mic][1], delayed_ffts[current_mic][1]);

		for (int i = 0; i < fft_size; ++i) {
			delay_and_sum_fft[0][i] = delay_and_sum_fft[0][i] + delayed_ffts[current_mic][0][i];
			delay_and_sum_fft[1][i] = delay_and_sum_fft[1][i] + delayed_ffts[current_mic][1][i];
		}

	}
	for (int i = 0; i < fft_size; ++i) {
			delay_and_sum_fft[0][i] = delay_and_sum_fft[0][i] / std::complex<double>(3.0, 3.0);
			delay_and_sum_fft[1][i] = delay_and_sum_fft[1][i] / std::complex<double>(3.0, 3.0);
	}
	
	get_inverse_fft(delay_and_sum_fft[0], delay_and_sum_output_ifft[0]);
	get_inverse_fft(delay_and_sum_fft[1], delay_and_sum_output_ifft[1]);
	
	for (int i = 0; i < nframes; ++i) {
		//delay_and_sum_output_ifft
		out[0][i] = delay_and_sum_output_ifft[0][fft1_offset + i] + delay_and_sum_output_ifft[1][fft2_offset + i];
	}

	// for (int i = 0; i < nframes; ++i) {
	// 	out[0][i] = filtered_ffts[current_mic][0][fft1_offset + i] + filtered_ffts[current_mic][1][fft2_offset + i];
	// }

	// UNDELAYED output
	// osb_array[testing_port_id].legacy_fill((double*)out[1], nframes, undelayed_offset);
	
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
	fft1_i = 0;
	fft2_i = 2 * nframes;
	fft1_offset = (int)(2.5 * (double)nframes);
	fft2_offset = (int)(0.5 * (double)nframes);
	undelayed_offset = (int)(2.5 * (double)nframes);

	// printf("%d %d %d\n", fft1_offset, fft2_offset, out2_offset);

	
	//preparing FFTW3 buffers
	i_fft = (std::complex<double> *) fftw_malloc(sizeof(std::complex<double>) * fft_size);
	i_time = (std::complex<double> *) fftw_malloc(sizeof(std::complex<double>) * fft_size);
	o_fft = (std::complex<double> *) fftw_malloc(sizeof(std::complex<double>) * fft_size);
	o_time = (std::complex<double> *) fftw_malloc(sizeof(std::complex<double>) * fft_size);
	
	i_forward = fftw_plan_dft_1d(fft_size, reinterpret_cast<fftw_complex*>(i_time), reinterpret_cast<fftw_complex*>(i_fft) , FFTW_FORWARD, FFTW_MEASURE);
	o_inverse = fftw_plan_dft_1d(fft_size, reinterpret_cast<fftw_complex*>(o_fft) , reinterpret_cast<fftw_complex*>(o_time), FFTW_BACKWARD, FFTW_MEASURE);
	
	for (i = 0; i < ports_number; ++i) {
		osb_array[i] = OSB<double>(buffer_size);
	}
	
	filtered_ffts = (double ***) malloc(ports_number * sizeof(double **));
	for (i = 0; i < ports_number; ++i){
		filtered_ffts[i] = (double **) malloc(2 * sizeof(double *));
		filtered_ffts[i][0] = (double *) calloc(fft_size, sizeof(double));
		filtered_ffts[i][1] = (double *) calloc(fft_size, sizeof(double));
	}

	delayed_ffts = (std::complex<double> ***) fftw_malloc(ports_number * sizeof(std::complex<double> **));
	for (i = 0; i < ports_number; ++i){
		delayed_ffts[i] = (std::complex<double> **) fftw_malloc(2 * sizeof(std::complex<double> *));
		delayed_ffts[i][0] = (std::complex<double> *) fftw_malloc(fft_size * sizeof(std::complex<double>));
		delayed_ffts[i][1] = (std::complex<double> *) fftw_malloc(fft_size * sizeof(std::complex<double>));
	}
	delay_and_sum_fft = (std::complex<double> **) fftw_malloc(sizeof(std::complex<double>*) * 2);
	delay_and_sum_fft[0] = (std::complex<double> *) fftw_malloc(sizeof(std::complex<double>) * fft_size);
	delay_and_sum_fft[1] = (std::complex<double> *) fftw_malloc(sizeof(std::complex<double>) * fft_size);
	
	delay_and_sum_output_ifft = (double **) malloc(2 * sizeof(double*));
	delay_and_sum_output_ifft[0] = (double *) malloc(fft_size * sizeof(double));
	delay_and_sum_output_ifft[1] = (double *) malloc(fft_size * sizeof(double));
	
	doa_value = (double *) malloc(ports_number * sizeof(double));
	doa_value[0] = 180.0;
	doa_value[1] = -90.0 - (-30.0);
	doa_value[2] = -150.0 - (-30.0);
	

	hann = (double *) calloc(fft_size, sizeof(double));
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
	delay_vector = (std::complex<double> *) malloc(fft_size * sizeof(std::complex<double>));

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
