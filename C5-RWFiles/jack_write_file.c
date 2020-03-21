/**
 * A simple wav reader, using libsndfile.
 */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <time.h>

#include <jack/jack.h>

//To read audio files
#include <sndfile.h>

jack_port_t *input_port;
jack_client_t *client;

//sndfile stuff
SNDFILE * audio_file;
SF_INFO audio_info;
unsigned int audio_position = 0;

double sample_rate;

/**
 * The process callback for this JACK application is called in a
 * special realtime thread once for each audio cycle.
 *
 * This client does nothing more than copy data from its input
 * port to its output port. It will exit when stopped by 
 * the user (e.g. using Ctrl-C on a unix-ish operating system)
 */
int jack_callback (jack_nframes_t nframes, void *arg){
	jack_default_audio_sample_t *in;
	float *write_buffer;
	int write_count;
	int i;
	
	//Preparing buffers
	in = (jack_default_audio_sample_t *)jack_port_get_buffer (input_port, nframes);
	write_buffer = (float *)malloc(nframes*sizeof(float));
	
	//Input from channels
	for (i=0;i<nframes;i++){
		write_buffer[i] = in[i];
	}
	
	//Read from file
	write_count = sf_write_float(audio_file,write_buffer,nframes);
	
	//Print Audio position
	audio_position += write_count;
	printf("\rAudio file in position: %d (%0.2f secs)", audio_position, (double)audio_position/sample_rate);
	
	//Check for writing error
	if(write_count != nframes){
		printf("\nEncountered I/O error. Exiting.\n");
		sf_close(audio_file);
		jack_client_close (client);
		exit (1);
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
	
	if(argc != 2){
		printf ("Audio File Path not provided.\n");
		exit(1);
	}
	
	//printf("%s %s\n",argv[1],argv[2]);
	char audio_file_path[100];
	strcpy(audio_file_path,argv[1]);
	
	const char *client_name = "write_audio_file";
	jack_options_t options = JackNoStartServer;
	jack_status_t status;
	
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
	sample_rate = (double)jack_get_sample_rate(client);
	printf ("Engine sample rate: %0.0f\n", sample_rate);
	
	/* create the agent input port */
	input_port = jack_port_register (client, "input_1", JACK_DEFAULT_AUDIO_TYPE,JackPortIsInput, 0);
	
	printf("Trying to open audio File: %s\n",audio_file_path);
	audio_info.samplerate = sample_rate;
	audio_info.channels = 1;
	audio_info.format = SF_FORMAT_WAV | SF_FORMAT_PCM_32;
	
	audio_file = sf_open (audio_file_path,SFM_WRITE,&audio_info);
	if(audio_file == NULL){
		printf("%s\n",sf_strerror(NULL));
		exit(1);
	}else{
		printf("Audio file info:\n");
		printf("\tSample Rate: %d\n",audio_info.samplerate);
		printf("\tChannels: %d\n",audio_info.channels);
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
	const char **serverports_names;
	 
	/* Assign our output port to a server input port*/
	// Find possible input server port names
	serverports_names = jack_get_ports (client, NULL, NULL, JackPortIsPhysical|JackPortIsOutput);
	if (serverports_names == NULL) {
		printf("no available physical record (server output) ports.\n");
		exit (1);
	}
	
	// Connect the first available to our output ports
	if (jack_connect (client, serverports_names[0], jack_port_name (input_port))) {
		printf ("cannot connect input port\n");
		exit(1);
	}

	// free serverports_names variable, we're not going to use it again
	free (serverports_names);
	
	printf ("done.\n");
	
	/* keep running until finished reading file */
	sleep(-1);
	
	jack_client_close (client);
	exit (0);
}
