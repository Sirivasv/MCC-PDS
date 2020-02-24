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

double complex  **is_time, **is_fft, **os_time, **os_fft;
fftw_plan *is_forward, *os_inverse;

jack_port_t **input_ports;
jack_port_t **output_ports;
jack_client_t *client;

double sample_rate;

double *freqs, **temp_outs, *hann_vals;
double complex  *exponentials;

const int NUM_BUFFS = 6;
const int NUM_WINDOWS = 4;
const int NUM_PORTS = 2;
const double MY_PI = 3.14159265358979323846;
int EXT_BUFF_SIZE;
double DELAY_TIME;
int idz = 0;
int *vids;

void hanning(double complex *vin) {
	for (int i = 0; i < EXT_BUFF_SIZE; ++i) {
		vin[i] *= hann_vals[i];
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
	jack_default_audio_sample_t **in, **out;
	int i, j;
	idz %= NUM_BUFFS;

	for (i = 0; i < NUM_BUFFS; ++i) {
		vids[i] = (idz + i) % NUM_BUFFS;
	}

	in = malloc(NUM_PORTS * sizeof(jack_default_audio_sample_t*));
    out = malloc(NUM_PORTS * sizeof(jack_default_audio_sample_t*));
    for (j = 0; j < NUM_PORTS; ++j) {
        in[j] = jack_port_get_buffer (input_ports[j], nframes);
        out[j] = jack_port_get_buffer (output_ports[j], nframes);
    }
	
	// Obteniendo la transformada de Fourier de este periodo
	for(i = 0; i < nframes; i++){
		// GUARDAMOS en 4tas posiciones de V<2>
		is_time[vids[2]][3072 + i] = in[0][i];

		// GUARDAMOS en 3ras posiciones de V<3>
		is_time[vids[3]][2048 + i] = in[0][i];

		// GUARDAMOS en 2das posiciones de V<4>
		is_time[vids[4]][1024 + i] = in[0][i];

		// GUARDAMOS en primeras posiciones de V<5>
		is_time[vids[5]][i] = in[0][i];
	}

	// HACEMOS HANN EN TEMP_OUTS V<2>
	hanning(is_time[vids[2]]);

	// TRANSORMAMOS V<2>
	fftw_execute(is_forward[vids[2]]);
	
	// PASAMOS A OUT FFT <2>
	for(i = 0; i < EXT_BUFF_SIZE; i++){
		os_fft[vids[2]][i] = is_fft[vids[2]][i];
	}

	// HACEMOS DESFASE en OUT FFT V<2>
	for(i = 0; i < EXT_BUFF_SIZE; i++){
		os_fft[vids[2]][i] *= exponentials[i];
	}

	// REGRESAMOS A TIEMPO
	fftw_execute(os_inverse[vids[2]]);

	// APLICAMOS LA DIVISION A LA SALIDA
	for(i = 0; i < EXT_BUFF_SIZE; i++){
		temp_outs[vids[2]][i] = creal(os_time[vids[2]][i])/(double)EXT_BUFF_SIZE;
	}

	// HACEMOS SUMA DE LOS INDICES INTERESANTES Y LOS GUARDAMOS EN OUT
	for(i = 0; i < nframes; i++){
		out[0][i] = temp_outs[vids[0]][i + 2559] + temp_outs[vids[2]][i + 511]; //fftw3 requiere normalizar su salida real de esta manera
		out[1][i] = creal(is_time[vids[0]][i + 2559]) + creal(is_time[vids[2]][i + 511]); //fftw3 requiere normalizar su salida real de esta manera
	}
	
	idz++;
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
	const char *client_name = "jack_fft";
	jack_options_t options = JackNoStartServer;
	jack_status_t status;
	int i;

	if (argc < 2) {
		printf("Need the delay time in seconds (0 to 0.011).\n");
		exit(1);
	}

	DELAY_TIME = atof(argv[1]);
	
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
	
	// vector ids
	vids = (int*) fftw_malloc(sizeof(int) * NUM_BUFFS);
	//preparing FFTW3 buffers
	is_fft = (double complex **) fftw_malloc(sizeof(double complex *) * NUM_BUFFS);
	is_time = (double complex **) fftw_malloc(sizeof(double complex *) * NUM_BUFFS);
	os_fft = (double complex **) fftw_malloc(sizeof(double complex *) * NUM_BUFFS);
	os_time = (double complex **) fftw_malloc(sizeof(double complex *) * NUM_BUFFS);

	temp_outs = (double **) fftw_malloc(sizeof(double *) * NUM_BUFFS);

	EXT_BUFF_SIZE = nframes * NUM_WINDOWS;

	// HAN VALUES
	hann_vals = (double *) fftw_malloc(sizeof(double) * EXT_BUFF_SIZE);
	int lim_HAN = EXT_BUFF_SIZE - 1;
	for (int i = 0; i <= lim_HAN; ++i) {
        double val = 0.5 * (1.0 - cos( 2.0 * MY_PI * ((double)i/(double)lim_HAN) ));
        hann_vals[i] = val;
    }

	for (i = 0; i < NUM_BUFFS; ++i) { 
		is_fft[i] = (double complex *) fftw_malloc(sizeof(double complex) * EXT_BUFF_SIZE );
		is_time[i] = (double complex *) fftw_malloc(sizeof(double complex) * EXT_BUFF_SIZE );
		os_fft[i] = (double complex *) fftw_malloc(sizeof(double complex) * EXT_BUFF_SIZE);
		os_time[i] = (double complex *) fftw_malloc(sizeof(double complex) * EXT_BUFF_SIZE);
		
		temp_outs[i] = (double *) fftw_malloc(sizeof(double) * EXT_BUFF_SIZE);
	}
	
	// PLANS
	is_forward = (fftw_plan *) fftw_malloc(sizeof(fftw_plan) * NUM_BUFFS);
	os_inverse = (fftw_plan *) fftw_malloc(sizeof(fftw_plan) * NUM_BUFFS);

	for (i = 0; i < NUM_BUFFS; ++i) {
		is_forward[i] = fftw_plan_dft_1d(EXT_BUFF_SIZE, is_time[i], is_fft[i] , FFTW_FORWARD, FFTW_MEASURE);
		os_inverse[i] = fftw_plan_dft_1d(EXT_BUFF_SIZE, os_fft[i], os_time[i], FFTW_BACKWARD, FFTW_MEASURE);
	}
	
	// vector for frequency indexes
	freqs = (double*)malloc(sizeof(double) * EXT_BUFF_SIZE);
	freqs[0] = 0.0;
	for (i = 1; i < EXT_BUFF_SIZE/2; ++i) {
		freqs[i] = ((double)i)*(sample_rate/(double)(EXT_BUFF_SIZE));
		freqs[EXT_BUFF_SIZE-i] = -freqs[i];
	}
	freqs[i] = ((double)i)*(sample_rate/(double)(EXT_BUFF_SIZE));

	// Vector for exponentials
	exponentials = (complex double *)malloc(sizeof(complex double) * EXT_BUFF_SIZE);
	for(i = 0; i < EXT_BUFF_SIZE; i++){
		exponentials[i] = cexp(-I*2.0*MY_PI*freqs[i]*DELAY_TIME);
	}

	/* allocate memory for ports */
    input_ports = malloc(NUM_PORTS * sizeof (jack_port_t*));
    output_ports = malloc(NUM_PORTS * sizeof (jack_port_t*));
    char port_name[10];

	/* create the agents for input port */
	for (int i = 0; i < NUM_PORTS; ++i) {
        sprintf(port_name, "input_%d", i);
        input_ports[i] = jack_port_register (client, port_name, JACK_DEFAULT_AUDIO_TYPE,JackPortIsInput, 0);
    }
	
	/* create the agents output port */
    for (int i = 0; i < NUM_PORTS; ++i) {
        sprintf(port_name, "output_%d", i);
        output_ports[i] = jack_port_register (client, port_name,JACK_DEFAULT_AUDIO_TYPE,JackPortIsOutput, 0);
    }
	
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
	
	
    for (int i = 0; i < NUM_PORTS; ++i) {

        // Connect the available ports to our input ports
        if (jack_connect (client, serverports_names[i], jack_port_name (input_ports[i]))) {
            printf("Cannot connect input port %d.\n", i);
            exit (1);
        }

    }

	// free serverports_names variable for reuse in next part of the code
	free (serverports_names);
	
	
	/* Assign our output port to a server input port*/
	// Find possible input server port names
	serverports_names = jack_get_ports (client, NULL, NULL, JackPortIsPhysical|JackPortIsInput);
	if (serverports_names == NULL) {
		printf("No available physical playback (server input) ports.\n");
		exit (1);
	}
	
	for (int i = 0; i < NUM_PORTS; ++i) {

        // Connect available to our output ports 
        if (jack_connect (client, jack_port_name (output_ports[i]), serverports_names[i])) {
            printf ("Cannot connect output port 1.\n");
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
