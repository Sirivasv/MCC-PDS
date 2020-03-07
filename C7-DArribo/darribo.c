/**
 * A simple example of how to do FFT with FFTW3 and JACK.
 */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <jack/jack.h>

// Include FFTW header
#include <complex.h> //needs to be included before fftw3.h for compatibility
#include <fftw3.h>

double complex *i_fft, *i_time, *o_fft, *o_time;
fftw_plan i_forward, o_inverse;

jack_port_t **input_port;
jack_port_t **output_port;
jack_default_audio_sample_t **in, **out;
jack_client_t *client;

double sample_rate;
double *hann;
const double MY_PI = 3.14159265358979323846264338327950288;

// CCV Requirements
double complex **signal_FFT, **signal_IFFT, **signal_CONJ, *signal_CCVFFT;
double *signal_NORM, **signal_CENTER, *resultCCV, *signal_CCVIFFT;

unsigned int ports_number, sub_buffer_size;
double *freqs;

void centerSignal(jack_default_audio_sample_t *data, double *return_data) {
    int i;

    // First we get the mean of the data
    // double mean_of_signal = 0;

    // for (i = 0; i < sub_buffer_size; ++i) {
    //     mean_of_signal += (double)data[i];
    // }

    // mean_of_signal /= (double)sub_buffer_size;

    // for (i = 0; i < sub_buffer_size; ++i) {
    //     return_data[i] -= mean_of_signal;
    // }

    for (i = 0; i < sub_buffer_size; ++i) {
        return_data[i] -= (double)data[i];
    }

}

void FFTSignal(double *data, double complex *return_data) {
    int i;

	for (i = 0; i < sub_buffer_size; ++i) {
		i_time[i] = data[i];
	}

	// FFT
	fftw_execute(i_forward);

	// return
	for (i = 0; i < sub_buffer_size; ++i) {
		return_data[i] = i_fft[i];
	}

}

void IFFTSignal(double complex *data, double *return_data) {
    int i;

	for (i = 0; i < sub_buffer_size; ++i) {
		o_fft[i] = data[i];
	}

	// IFFT
	fftw_execute(o_inverse);

	// return
	for (i = 0; i < sub_buffer_size; ++i) {
		return_data[i] = creal(o_time[i])/(double)sub_buffer_size;
	}

}

void conjugateCompSignal(double complex *data, double complex *return_data) {
    int i;

	for (i = 0; i < sub_buffer_size; ++i) {
		return_data[i] = conj(data[i]);
	}

}

void pointProduct(double complex *xf, double complex *y_conj_f, double complex *return_data) {

    int i;

    for (i = 0; i < sub_buffer_size; ++i) {
		return_data[i] = xf[i] * y_conj_f[i];
	}

}

double getNormSignalInR(double *data) {

    int i;
	double return_data;

    for (i = 0; i < sub_buffer_size; ++i) {
		return_data += data[i] * data[i];
	}

	return sqrt(return_data);
}

// void filter(double *data) {
// 	int i;

// 	// Hann
// 	for (i = 0; i < sub_buffer_size; ++i) {
// 		i_time[i] = data[i] * hann[i];
// 	}

// 	// FFT
// 	fftw_execute(i_forward);

// 	// Filter
// 	for (i = 0; i < sub_buffer_size; ++i) {
// 		o_fft[i] = i_fft[i] * delay_vector[i];
// 	}

// 	// ifft
// 	fftw_execute(o_inverse);

// 	// return
// 	for(i = 0; i < sub_buffer_size; i++){
// 		data[i] = creal(o_time[i])/(double)sub_buffer_size; //fftw3 requiere normalizar su salida real de esta manera
// 	}
// }

/**
 * The process callback for this JACK application is called in a
 * special realtime thread once for each audio cycle.
 *
 * This client does nothing more than copy data from its input
 * port to its output port. It will exit when stopped by 
 * the user (e.g. using Ctrl-C on a unix-ish operating system)
 */
int jack_callback (jack_nframes_t nframes, void *arg){
	int i, port_i;
    for (port_i = 0; port_i < ports_number; ++port_i) {
        in[port_i] = (jack_default_audio_sample_t *)jack_port_get_buffer (input_port[port_i], nframes);
	    out[port_i] = (jack_default_audio_sample_t *)jack_port_get_buffer (output_port[port_i], nframes);
        for (i = 0; i < nframes; ++i) {
            out[port_i][i] = in[port_i][i];
        }
    }
    
    
    for (port_i = 0; port_i < ports_number; ++port_i) {
        
        // We center the signals 
        centerSignal(in[port_i], signal_CENTER[port_i]); 
    
        // We get the centererd signals fourier transform
        FFTSignal(signal_CENTER[port_i], signal_FFT[port_i]);

        // We get the norm of the signals
        signal_NORM[port_i] = getNormSignalInR(signal_CENTER[port_i]);

    }

    // We get the conjugate from the second signal
    conjugateCompSignal(signal_FFT[1], signal_CONJ[1]);

    // We get the point product of fft from first signal and the conjugate of the second
    pointProduct(signal_FFT[0], signal_CONJ[1], signal_CCVFFT);

	// We get the inverse transform of CCVF
	IFFTSignal(signal_CCVFFT, signal_CCVIFFT);

	// We point divide the Inverse by the norm product
	double norm_product = signal_NORM[0] * signal_NORM[1];
	double max_val = signal_CCVIFFT[0] / norm_product;
	int max_i = 0;
	for (i = 0; i < sub_buffer_size; ++i) {
		resultCCV[i] = signal_CCVIFFT[i] / norm_product;
		if (resultCCV[i] > max_val) {
			max_i = i;
			max_val = resultCCV[i];
		}
	}
	
	printf("******************************************\n");
	printf("max_val = %lf /// max_i = %d\n", max_val, max_i);
	int delay;
	if ( max_i < (sub_buffer_size / 2)) {
		delay = max_i - 1;
	} else {
		delay = sub_buffer_size - max_i - 1;
	}

	printf("DESFASE = %d\n", delay);
	double delay_in_seconds = ((double)delay) / 44100.0;
	printf("DESFASE en segundos = %lf\n", delay_in_seconds);
	double angle_of_direction = (asin((343.0 * delay_in_seconds) / 0.18) * 180.0) / MY_PI;
	printf("DIRECCION DE ARRIBO en grados = %lf\n", angle_of_direction);
	fflush(stdout);
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
	const char *client_name = "darribo";
	jack_options_t options = JackNoStartServer;
	jack_status_t status;
	int i;
	
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
	sub_buffer_size = nframes;
    ports_number = 2;


	//preparing FFTW3 buffers
	i_fft = (double complex *) fftw_malloc(sizeof(double complex) * sub_buffer_size);
	i_time = (double complex *) fftw_malloc(sizeof(double complex) * sub_buffer_size);
	o_fft = (double complex *) fftw_malloc(sizeof(double complex) * sub_buffer_size);
	o_time = (double complex *) fftw_malloc(sizeof(double complex) * sub_buffer_size);
	
	i_forward = fftw_plan_dft_1d(sub_buffer_size, i_time, i_fft , FFTW_FORWARD, FFTW_MEASURE);
	o_inverse = fftw_plan_dft_1d(sub_buffer_size, o_fft , o_time, FFTW_BACKWARD, FFTW_MEASURE);

    // CCV
	signal_CONJ = (double complex **) fftw_malloc(sizeof(double complex*) * ports_number);
	signal_FFT = (double complex **) fftw_malloc(sizeof(double complex*) * ports_number);
	signal_IFFT = (double complex **) fftw_malloc(sizeof(double complex*) * ports_number);
    signal_CENTER = (double **) fftw_malloc(sizeof(double*) * ports_number);
    signal_NORM = (double *) fftw_malloc(sizeof(double*) * ports_number);
    
    for (i = 0; i < ports_number; ++i) {

        // CCV
        signal_CONJ[i] = (double complex *) fftw_malloc(sizeof(double complex) * sub_buffer_size);
        signal_FFT[i] = (double complex *) fftw_malloc(sizeof(double complex) * sub_buffer_size);
        signal_IFFT[i] = (double complex *) fftw_malloc(sizeof(double complex) * sub_buffer_size);
        signal_CENTER[i] = (double *) fftw_malloc(sizeof(double) * sub_buffer_size);
        signal_CCVFFT = (double complex *) fftw_malloc(sizeof(double complex) * sub_buffer_size);
        signal_CCVIFFT = (double *) fftw_malloc(sizeof(double) * sub_buffer_size);
        resultCCV = (double *) fftw_malloc(sizeof(double) * sub_buffer_size);

    }

	hann = (double *) calloc(sub_buffer_size, sizeof(double));
	freqs = (double *) calloc(sub_buffer_size, sizeof(double));


	// Calculate hann window
	for (i = 0; i < sub_buffer_size;++i)
		hann[i] = 0.5 * (1 - cos( (2.0*MY_PI*i)/(sub_buffer_size - 1) ));

	// vector for frequency indexes
	freqs[0] = 0.0;
	freqs[sub_buffer_size/2] = sample_rate/2.00;
	int lim = sub_buffer_size/2;
	for (i = 1; i < lim; ++i) {
		freqs[i] = ((double)i)*(sample_rate/(double)(sub_buffer_size));
		freqs[sub_buffer_size-i] = -freqs[i];
	}

	// Create agent ports
    input_port = malloc(ports_number * sizeof (jack_port_t*));
    output_port = malloc(ports_number * sizeof (jack_port_t*));
    in = malloc(ports_number * sizeof (jack_default_audio_sample_t*));
    out = malloc(ports_number * sizeof (jack_default_audio_sample_t*));
    char port_name[10];

    for (i = 0; i < ports_number; ++i) {

        /* create the agent input ports */
        sprintf(port_name, "input_%d", (i+1));
        input_port[i] = jack_port_register (client, port_name, JACK_DEFAULT_AUDIO_TYPE,JackPortIsInput, 0);

        /* create the agent output port */
        sprintf(port_name, "output_%d", (i+1));
	    output_port[i] = jack_port_register (client, port_name, JACK_DEFAULT_AUDIO_TYPE,JackPortIsOutput, 0);    

    } 
	
	/* check that both ports were created succesfully */
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

    // We are not connecting our inputs to the system mics

	// serverports_names = jack_get_ports (client, NULL, NULL, JackPortIsPhysical|JackPortIsOutput);
	// if (serverports_names == NULL) {
	// 	printf("No available physical capture (server output) ports.\n");
	// 	exit (1);
	// }
    // for (i = 0; i < ports_number; ++i) {
    //     // Connect the first available to our input port
    //     if (jack_connect (client, serverports_names[i], jack_port_name (input_port[i]))) {
    //         printf("Cannot connect input port.\n");
    //         exit (1);
    //     }
    // }
	// // free serverports_names variable for reuse in next part of the code
	// free (serverports_names);
	
	
	/* Assign our output port to a server input port*/
	// Find possible input server port names
	serverports_names = jack_get_ports (client, NULL, NULL, JackPortIsPhysical|JackPortIsInput);
	if (serverports_names == NULL) {
		printf("No available physical playback (server input) ports.\n");
		exit (1);
	}

    for (i = 0; i < ports_number; ++i) {

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
