/**
 * A simple 1-input to 1-output JACK client.
 */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

#include <jack/jack.h>


// Doble apuntador
jack_port_t **input_ports;
jack_port_t **output_ports;
jack_client_t *client;

const unsigned int NUM_PORTS = 2;

/**
 * The process callback for this JACK application is called in a
 * special realtime thread once for each audio cycle.
 *
 * This client does nothing more than copy data from its input
 * port to its output port. It will exit when stopped by 
 * the user (e.g. using Ctrl-C on a unix-ish operating system)
 */
int jack_callback (jack_nframes_t nframes, void *arg) {

	int i, j;
    jack_default_audio_sample_t **in, **out;	
    
    in = malloc(NUM_PORTS * sizeof(jack_default_audio_sample_t*));
    out = malloc(NUM_PORTS * sizeof(jack_default_audio_sample_t*));
    for (j = 0; j < NUM_PORTS; ++j) {
        in[j] = jack_port_get_buffer (input_ports[j], nframes);
        out[j] = jack_port_get_buffer (output_ports[j], nframes);
    }
	
    for (j = 0; j < NUM_PORTS; ++j) {
        for (i = 0; i < nframes; ++i) {
            out[j][i] = in[j][i];
        }
    }
	
    // memcpy (out_1, in_1, nframes * sizeof (jack_default_audio_sample_t));
	
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
	const char *client_name = "in_to_out_mod_1";
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
	printf ("Engine sample rate: %d\n", jack_get_sample_rate (client));
	
	/* allocate memory for ports */
    input_ports = malloc(NUM_PORTS * sizeof (jack_port_t*));
    output_ports = malloc(NUM_PORTS * sizeof (jack_port_t*));
    char port_name[10];

	/* create the agents for input port */
	for (int i = 0; i < NUM_PORTS; ++i) {
        sprintf(port_name, "input_%d", (i+1));
        input_ports[i] = jack_port_register (client, port_name, JACK_DEFAULT_AUDIO_TYPE,JackPortIsInput, 0);
    }
	
	/* create the agents output port */
    for (int i = 0; i < NUM_PORTS; ++i) {
        sprintf(port_name, "output_%d", (i+1));
        output_ports[i] = jack_port_register (client, port_name,JACK_DEFAULT_AUDIO_TYPE,JackPortIsOutput, 0);
    }
	
	/* check that both ports were created succesfully */
    // Put this if on the fors
	// if ((input_port_1 == NULL) || (input_port_2 == NULL) || (output_port_1 == NULL) || (output_port_2 == NULL)) {
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
