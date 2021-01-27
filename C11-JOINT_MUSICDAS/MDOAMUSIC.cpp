/**
 * A simple example of how to do FFT with FFTW3 and JACK.
 */

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
// To read/write files
#include <sndfile.h>


std::complex<double> *i_fft, *i_time, *o_fft, *o_time;
fftw_plan i_forward, o_inverse;

jack_port_t **input_port;
jack_port_t **output_port;
jack_default_audio_sample_t **in, **out;
jack_client_t *client;

double sample_rate;
double *hann;
const double MY_PI = 3.14159265358979323846264338327950288;
const std::complex<double> COMP_I(0.0, 1.0);
const int WINDOW_MEM_LEN = 5;
const double MIC_DIST = 0.18;
const double SOUND_SPEED = 343.0;

const double FREQ_LOW_RANGE = 20;
const double FREQ_HIGH_RANGE = 800;
const int FREQS_RES_LEN = 16;

const int OUT_SYSTEM_PORTS = 2;
const int INPUT_SIGNALS = 2;
const int NOISE_SPACES = 1;
const double ANGLE_LOW_RANGE = -90.0;
const double ANGLE_HIGH_RANGE = 100.0;
const double ANGLE_STEP = 5.0;

// FULL
// const double FREQ_LOW_RANGE = 0;
// const double FREQ_HIGH_RANGE = 22050;
// const int FREQS_RES_LEN = 100;

// VOICE
// const double FREQ_LOW_RANGE = 100;
// const double FREQ_HIGH_RANGE = 260;
// const int FREQS_RES_LEN = 4;

// MUSIC  Requirements
std::complex<double> ***signal_FFT_MEM;
int window_mem_current_id = 0;
std::vector<int> freqs_res_values;
std::vector<double> angle_values;
std::vector<std::vector<std::complex<double>>> angle_delays2;
std::vector<std::vector<std::complex<double>>> angle_delays3;

std::complex<double> **signal_FFT, **signal_IFFT, **signal_CONJ, *signal_CONJPROD, *signal_ABSCONJPROD, *signal_CCVPHATFFT;
double *signal_NORM, **signal_CENTER, **signal_WINDOW, *resultCCV, *signal_CCVIFFT;

unsigned int sub_buffer_size;
const int ports_number = 3;
double *freqs;
double *freqs_energies;

/* DAS */
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

double *doa_angle, *doa_value;
OSB<double> osb_array[ports_number];
double *freqs_das;
unsigned int fft_size_das;
std::complex<double> **delay_vector;
std::complex<double> **delay_and_sum_fft;
double ***input_time_values;
unsigned int fft1_i, fft2_i, fft1_offset, fft2_offset, undelayed_offset;
std::complex<double> ***window_ffts;
std::complex<double> ***delayed_ffts;
std::complex<double> *i_fft_das, *i_time_das, *o_fft_das, *o_time_das;
fftw_plan i_forward_das, o_inverse_das;
double *hann_das;
double **delay_and_sum_output_ifft;
// Write file required data
//sndfile stuff
SNDFILE ** audio_file;
SF_INFO *audio_info;
unsigned int * audio_position;
unsigned int buffer_size;


void get_fft(double *data_in, std::complex<double> *data_out) {
	
	// Hann
	for (int i = 0; i < fft_size_das; ++i) {
		i_time_das[i] = data_in[i] * hann_das[i];
	}

	// FFT
	fftw_execute(i_forward_das);

	// Filter
	for (int i = 0; i < fft_size_das; ++i) {
		data_out[i] = i_fft_das[i];
	}

}

void get_delayed_fft(int current_mic, std::complex<double> *data_in, std::complex<double> *data_out) {
	
	for (int i = 0; i < fft_size_das; ++i) {
		data_out[i] = data_in[i] * delay_vector[current_mic][i];
	}

}

void get_inverse_fft(std::complex<double> *data_in, double *data_out) {
	
	for (int i = 0; i < fft_size_das; ++i) {
	 o_fft_das[i] = data_in[i];
	}

	// ifft
	fftw_execute(o_inverse_das);

	// return
	for(int i = 0; i < fft_size_das; i++){
		data_out[i] = std::real(o_time_das[i])/(double)fft_size_das; //fftw3 requiere normalizar su salida real de esta manera
	}

}

/*MUSIC*/

void applyHann(jack_default_audio_sample_t *data, double *return_data) {
    int i;

    for (i = 0; i < sub_buffer_size; ++i) {
        return_data[i] = data[i] * hann[i];
    }

}

void centerSignal(double *data, double *return_data) {
    int i;

    // First we get the mean of the data
    double mean_of_signal = 0;

    for (i = 0; i < sub_buffer_size; ++i) {
        mean_of_signal += (double)data[i];
    }

    mean_of_signal /= (double)sub_buffer_size;

    for (i = 0; i < sub_buffer_size; ++i) {
        return_data[i] = data[i] - mean_of_signal;
    }

	// Why there is no much difference?
    // for (i = 0; i < sub_buffer_size; ++i) {
    //     return_data[i] = (double)data[i];
    // }

}

void FFTSignal(double *data, std::complex<double> *return_data) {
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

void IFFTSignal(std::complex<double> *data, double *return_data) {
    int i;

	for (i = 0; i < sub_buffer_size; ++i) {
		o_fft[i] = data[i];
	}

	// IFFT
	fftw_execute(o_inverse);

	// return
	for (i = 0; i < sub_buffer_size; ++i) {
		return_data[i] = std::real(o_time[i])/(double)sub_buffer_size;
	}

}

void conjugateCompSignal(std::complex<double> *data, std::complex<double> *return_data) {
    int i;

	for (i = 0; i < sub_buffer_size; ++i) {
		return_data[i] = std::conj(data[i]);
	}

}

void pointProduct(std::complex<double> *xf, std::complex<double> *y_conj_f, std::complex<double> *return_data) {

    int i;

    for (i = 0; i < sub_buffer_size; ++i) {
		return_data[i] = xf[i] * y_conj_f[i];
	}

}

void pointDivision(std::complex<double> *dividend_s, std::complex<double> *divisor_s, std::complex<double> *return_data) {

    int i;

    for (i = 0; i < sub_buffer_size; ++i) {
		return_data[i] = (std::real(dividend_s[i]) / std::real(divisor_s[i]) ) + ( (std::imag(dividend_s[i])) * COMP_I);
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

void getAbs(std::complex<double> *data, std::complex<double> *return_data) {
	int i;

    for (i = 0; i < sub_buffer_size; ++i) {
		return_data[i] = std::abs(data[i]);
	}

}

std::vector<double> get_peaks_in_array(int many_peaks, std::vector<double> value_array) {
	
	int lim_varray = value_array.size();
	
	std::vector<double> peaks_to_return;
	for (int i = 0; i < many_peaks; ++i) {
		
		// Get max value in array to determine peak
		double greatest_music_value = -1e25;
		int greatest_music_index = 0;
		for (int pos_i = 0; pos_i < lim_varray; ++pos_i) {

			if (value_array[pos_i] > greatest_music_value) {
				greatest_music_value = value_array[pos_i];
				greatest_music_index = pos_i;
			}

		}

		// Insert peak at vector
		peaks_to_return.push_back(angle_values[greatest_music_index]);

		// Delete peak values to the left
		double previous_value = value_array[greatest_music_index];
		for (int pos_i = greatest_music_index - 1; pos_i >= 0; --pos_i) {
			if (value_array[pos_i] <= previous_value) {
				previous_value = value_array[pos_i];
				value_array[pos_i] = -1;
			} else {
				break;
			}
		}

		// Delete peak values to the right
		previous_value = value_array[greatest_music_index];
		for (int pos_i = greatest_music_index + 1; pos_i < lim_varray; ++pos_i) {
			if (value_array[pos_i] <= previous_value) {
				previous_value = value_array[pos_i];
				value_array[pos_i] = -1;
			} else {
				break;
			}
		}

		// supress peak
		value_array[greatest_music_index] = -1;
	}
	
	return peaks_to_return;
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
	int i, j, port_i;
	for (port_i = 0; port_i < ports_number; ++port_i) {
		in[port_i] = (jack_default_audio_sample_t *)jack_port_get_buffer (input_port[port_i], nframes);
		out[port_i] = (jack_default_audio_sample_t *)jack_port_get_buffer (output_port[port_i], nframes);
		// for (i = 0; i < nframes; ++i) {
		// 	out[port_i][i] = in[port_i][i];
		// }
	}
    
	for (port_i = 0; port_i < ports_number; ++port_i) {
			
		// We Apply a Hann Window
		applyHann(in[port_i], signal_WINDOW[port_i]);

		// We center the signals 
		centerSignal(signal_WINDOW[port_i], signal_CENTER[port_i]); 

		// We get the centered signals fourier transform
		FFTSignal(signal_CENTER[port_i], signal_FFT[port_i]);

		// We save the FFT Matrix in memeory
		// for (int j = 0; j < freqs_res_values.size(); ++j) {
		for (int j = 0; j < (sub_buffer_size/2 + 1); ++j) {
			signal_FFT_MEM[window_mem_current_id][port_i][j] = i_fft[j];
			// freqs_energies[j] += abs(i_fft[j]);
		}

	}

	// We save the FFT Matrix in memeory
	// for (int j = 0; j < freqs_res_values.size(); ++j) {
	for (int j = 0; j < (sub_buffer_size/2 + 1); ++j) {
		for (int window_iter = 0; window_iter < WINDOW_MEM_LEN; ++window_iter) {	
			for (port_i = 0; port_i < ports_number; ++port_i) {
				freqs_energies[j] += abs(signal_FFT_MEM[window_iter][port_i][j]);
			}
		}
	}
		


	// Get relevant freqs
	std::vector<std::pair<double, int>> freqs_energies_to_sort;
	for (int j = 0; j < (sub_buffer_size/2 + 1); ++j) {
		freqs_energies_to_sort.push_back(std::pair<double, int>(freqs_energies[j], j));
		freqs_energies[j] = 0.0;
	}
	std::sort (freqs_energies_to_sort.begin(), freqs_energies_to_sort.end());
	std::reverse(freqs_energies_to_sort.begin(),freqs_energies_to_sort.end());
	freqs_res_values.clear();
	// std::cout << "++ FREQS\n";
	// std::cout << "[";
	for (int freq_res_i = 0; freq_res_i < freqs_energies_to_sort.size(); ++freq_res_i) {
		double curr_freq_value = freqs[freqs_energies_to_sort[freq_res_i].second];
		if ((curr_freq_value >= FREQ_LOW_RANGE) && (curr_freq_value <= FREQ_HIGH_RANGE)) {
			freqs_res_values.push_back(freqs_energies_to_sort[freq_res_i].second);
			// std::cout << curr_freq_value << ',';
			if (freqs_res_values.size() == FREQS_RES_LEN) {
				break;
			}
		}	
	}
	// std::cout << "]\n";
	// std::cout << "++ ANGLES\n";
	// std::cout << "[";
	// for (int i = 0; i < (int)angle_values.size(); ++i) {
		// std::cout << angle_values[i] << ", ";
	// }
	// std::cout << "]\n";

	std::vector<std::vector<double>> music_spectrum;
	double mat_val = 0.0;
	double mint_val = 0.0;
	double mt_cnt = 0.0;

	for (int freq_i = 0; freq_i < freqs_res_values.size(); ++freq_i){

		Eigen::MatrixXcd X_MATR(ports_number,WINDOW_MEM_LEN);
		for (int mic_i = 0; mic_i < ports_number; ++mic_i) {
			
			for (j = 0; j < WINDOW_MEM_LEN; ++j) {

				int curr_saved_window_id = window_mem_current_id - j;
				if (curr_saved_window_id < 0) {
					curr_saved_window_id += WINDOW_MEM_LEN;
				}

			
				X_MATR(mic_i,curr_saved_window_id) = signal_FFT_MEM[window_mem_current_id][mic_i][freqs_res_values[freq_i]];
			}
		}

		Eigen::MatrixXcd R_MAT(ports_number, ports_number);
		R_MAT = X_MATR * X_MATR.adjoint();

		Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es_c(R_MAT);
		// std::cout << "--- The eigenvalues of C are:\n" << es_c.eigenvalues() << std::endl << std::endl;
		// std::cout << "--- The abs of C(0) are:\n" << abs(es_c.eigenvalues()(0))<< std::endl << std::endl;
		// std::cout << "--- The abs of C(1) are:\n" << abs(es_c.eigenvalues()(1))<< std::endl << std::endl;
		// std::cout << "--- The abs of C(2) are:\n" << abs(es_c.eigenvalues()(2))<< std::endl << std::endl;
		// std::cout << "--- The eigenvectors of C are (one vector per column):\n" << es_c.eigenvectors() << std::endl << std::endl;
		// std::cout << "--- The first eigenvector:\n" << es_c.eigenvectors().col(0) << std::endl << std::endl;
		
		std::vector<std::pair<double, int>> eigen_values_to_sort;
		for (int mic_i = 0; mic_i < ports_number; ++mic_i) {
			eigen_values_to_sort.push_back(std::pair<double, int>(abs(es_c.eigenvalues()(mic_i)), mic_i));
		}
		std::sort (eigen_values_to_sort.begin(), eigen_values_to_sort.end());
		std::reverse(eigen_values_to_sort.begin(),eigen_values_to_sort.end());
		// std::cout << "eigen_values_to_sort contains:";
		// for (std::vector<std::pair<double, int>>::iterator it=eigen_values_to_sort.begin(); it!=eigen_values_to_sort.end(); ++it){
			// std::cout << "( " << (*it).first << ' ' << (*it).second << " )";
		// }
		// std::cout << '\n';
		Eigen::MatrixXcd Qs(ports_number, INPUT_SIGNALS);
		for (int i = 0; i < ports_number; ++i) {
			for (int j = 0; j < INPUT_SIGNALS; ++j) {
				Qs(i,j) = es_c.eigenvectors().col(eigen_values_to_sort[j].second)(i);
			}
		}
		Eigen::MatrixXcd Qn(ports_number, NOISE_SPACES);
		for (int i = 0; i < ports_number; ++i) {
			for (int j = 0; j < NOISE_SPACES; ++j) {
				Qn(i,j) = es_c.eigenvectors().col(eigen_values_to_sort[((int)eigen_values_to_sort.size()) - j - 1].second)(i);
			}
		}
		
		// std::cout << "--- Complex Number Matrix Qs:\n" << Qs << std::endl << std::endl;
		// std::cout << "--- Complex Number Matrix Qn:\n" << Qn << std::endl << std::endl;
		// std::cout << "--- Complex Number Matrix X_MATR:\n" << X_MATR << std::endl << std::endl;

		// compute MUSIC spectrum
		std::vector<double> current_spectrum;
		
		//std::cout << "[";
		std::vector<double> best_aod;
		double previous_music_value = -1e25;
		double previous_angle = 0.0;
		int flag = 1;
		for (int i = 0; i < (int)angle_values.size(); ++i) {
			Eigen::MatrixXcd a1(ports_number, 1);
			a1(0,0) = std::complex<double>(1.0, 0.0);
			a1(1,0) = angle_delays2[freqs_res_values[freq_i]][i];
			a1(2,0) = angle_delays3[freqs_res_values[freq_i]][i];
			
			
			double current_music_value = std::real((a1.adjoint() * a1)(0,0)) / std::real((a1.adjoint()*Qn*Qn.adjoint()*a1)(0,0));
			
			if (flag == 1) { // scaling up
				if (current_music_value < previous_music_value) {
					best_aod.push_back(previous_angle);
					flag = 0;
				}
			} else { // scaling down
				if (current_music_value > previous_music_value) {
					flag = 1;
				}
			}
			previous_music_value = current_music_value;
			previous_angle = angle_values[i];
			//std::cout << current_music_value << ", ";
			current_spectrum.push_back(current_music_value);
		}
		if (best_aod.size() == 0) {
			best_aod.push_back(previous_angle);
			best_aod.push_back(previous_angle);
		}
		if (best_aod.size() == 1) {
			best_aod.push_back(previous_angle);
		}
		// std::cout << "]\n";
		// printf("DIRECCION DE ARRIBO 1 en grados = %lf\n", best_aod[0]);
		// printf("DIRECCION DE ARRIBO 2 en grados = %lf\n", best_aod[1]);
		// printf("FRECUENCIA ACTUAL = %lf\n", freqs[freqs_res_values[freq_i]]);
		music_spectrum.push_back(current_spectrum);
		if (best_aod[0] < best_aod[1]) {
			mat_val += best_aod[1];
			mint_val += best_aod[0];
		} else {
			mat_val += best_aod[0];
			mint_val += best_aod[1];
		}
		mt_cnt += 1.0;
	}

	printf("DIRECCION DE ARRIBO 1 en grados = %lf\n", mint_val / mt_cnt);
	printf("DIRECCION DE ARRIBO 2 en grados = %lf\n", mat_val / mt_cnt);
	window_mem_current_id += 1;
	if (window_mem_current_id >= WINDOW_MEM_LEN) window_mem_current_id = 0;
	

	// std::vector<double> compressed_music;
	// int lim_spectrums = music_spectrum.size();
	// int lim_per_spectrum = music_spectrum[0].size();
	// // std::cout << "[";
	// for (int pos_i = 0; pos_i < lim_per_spectrum; ++pos_i) {
	// 	double current_value_i = 0.0;
	// 	for (int spectrum_i = 0; spectrum_i < lim_spectrums; ++spectrum_i) {
	// 		current_value_i += music_spectrum[spectrum_i][pos_i];
	// 	}
	// 	compressed_music.push_back(current_value_i / (double)lim_spectrums);
	// 	compressed_music.push_back(current_value_i);
	// 	// std::cout << compressed_music[pos_i] << ", ";
	// }
	// // std::cout << "]\n";
	
	// std::vector<double> peaks = get_peaks_in_array(2, compressed_music);
	
	// printf("DIRECCION DE ARRIBO MEAN 1 en grados = %lf\n", peaks[0]);
	// printf("DIRECCION DE ARRIBO MEAN 2 en grados = %lf\n", peaks[1]);

	
	/* DAS */
	
	// doa_angle[0] = mint_val / mt_cnt;
	// doa_angle[1] = mat_val / mt_cnt;
	
	doa_angle[0] = 90.0;
	doa_angle[1] = -30.0;

	// Write file buffer
	float **write_buffer_doa;
	int *write_count_doa;
	write_count_doa = (int*)malloc(2 * sizeof(int));
	write_buffer_doa = (float**)malloc(2 * sizeof(float*));
	write_buffer_doa[0] = (float *)malloc(nframes*sizeof(float));
	write_buffer_doa[1] = (float *)malloc(nframes*sizeof(float));


	for (int current_mic = 0; current_mic < ports_number; current_mic++) {
		// copy in to las nframes of buffer
		for (int i = 0; i < nframes; ++i) {
			osb_array[current_mic].insert_element(in[current_mic][i]);
		}
	}

	
	// MAIN LOOP
	for (int curr_angle = 0; curr_angle < 2; ++curr_angle) {
		
		// Direction value
		doa_value[0] = (-MIC_DIST/SOUND_SPEED)*cos((-90.0 - doa_angle[curr_angle])*M_PI/180.0);
		doa_value[1] = (-MIC_DIST/SOUND_SPEED)*cos((-150.0 - doa_angle[curr_angle])*M_PI/180.0);

		// Calculate delay vector
		for (int current_mic = 1; current_mic < ports_number; ++current_mic){
			
			for (int i = 0; i < fft_size_das; ++i) {
				std::complex<double> freq_res_complex = std::complex<double>(freqs_das[i],0.0);
				delay_vector[current_mic][i] = std::exp(-COMP_I*2.0*M_PI*freq_res_complex*doa_value[current_mic - 1]);
			}

		}
		

		// Reset to zero
		for (int i = 0; i < fft_size_das; ++i) {
			delay_and_sum_fft[0][i] = 0;
			delay_and_sum_fft[1][i] = 0;
		}


		for (int current_mic = 0; current_mic < ports_number; current_mic++) {

			// Fill fft1 y fft2
			osb_array[current_mic].legacy_fill(input_time_values[current_mic][0], fft_size_das, fft1_i);
			osb_array[current_mic].legacy_fill(input_time_values[current_mic][1], fft_size_das, fft2_i);

			// Get FFT and delay input per direction
			get_fft(input_time_values[current_mic][0], window_ffts[current_mic][0]);
			get_delayed_fft(current_mic, window_ffts[current_mic][0], delayed_ffts[current_mic][0]);

			get_fft(input_time_values[current_mic][1], window_ffts[current_mic][1]);
			get_delayed_fft(current_mic, window_ffts[current_mic][1], delayed_ffts[current_mic][1]);
			
			// Sum FFTS
			for (int i = 0; i < fft_size_das; ++i) {
				delay_and_sum_fft[0][i] = delay_and_sum_fft[0][i] + delayed_ffts[current_mic][0][i];
				delay_and_sum_fft[1][i] = delay_and_sum_fft[1][i] + delayed_ffts[current_mic][1][i];
			}

		}

		for (int i = 0; i < fft_size_das; ++i) {
			delay_and_sum_fft[0][i] = delay_and_sum_fft[0][i] / std::complex<double>(3.0, 3.0);
			delay_and_sum_fft[1][i] = delay_and_sum_fft[1][i] / std::complex<double>(3.0, 3.0);
		}

		get_inverse_fft(delay_and_sum_fft[0], delay_and_sum_output_ifft[0]);
		get_inverse_fft(delay_and_sum_fft[1], delay_and_sum_output_ifft[1]);
		
		for (int i = 0; i < nframes; ++i) {
			//delay_and_sum_output_ifft
			out[curr_angle][i] = delay_and_sum_output_ifft[0][fft1_offset + i] + delay_and_sum_output_ifft[1][fft2_offset + i];
			write_buffer_doa[curr_angle][i] = out[curr_angle][i];
		}

		

		// Write to file
		write_count_doa[curr_angle] = sf_write_float(audio_file[curr_angle], write_buffer_doa[curr_angle], nframes);

		//Print Audio position
		audio_position[curr_angle] += write_count_doa[curr_angle];
		// printf("\rAudio file in position: %d (%0.2f secs)", audio_position[curr_angle], (double)audio_position[curr_angle]/sample_rate);

		//Check for writing error
		if(write_count_doa[curr_angle] != nframes){
			printf("\nEncountered I/O error in file %d. Exiting.\n", curr_angle);
			sf_close(audio_file[curr_angle]);
			jack_client_close (client);
			exit (1);
		}

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


	//preparing FFTW3 buffers
	i_fft = (std::complex<double> *) fftw_malloc(sizeof(std::complex<double>) * sub_buffer_size);
	i_time = (std::complex<double> *) fftw_malloc(sizeof(std::complex<double>) * sub_buffer_size);
	o_fft = (std::complex<double> *) fftw_malloc(sizeof(std::complex<double>) * sub_buffer_size);
	o_time = (std::complex<double> *) fftw_malloc(sizeof(std::complex<double>) * sub_buffer_size);
	
	i_forward = fftw_plan_dft_1d(sub_buffer_size, reinterpret_cast<fftw_complex*>(i_time), reinterpret_cast<fftw_complex*>(i_fft), FFTW_FORWARD, FFTW_MEASURE);
	o_inverse = fftw_plan_dft_1d(sub_buffer_size, reinterpret_cast<fftw_complex*>(o_fft) , reinterpret_cast<fftw_complex*>(o_time), FFTW_BACKWARD, FFTW_MEASURE);

    // CCV
	signal_CONJ = (std::complex<double> **) fftw_malloc(sizeof(std::complex<double>*) * ports_number);
	signal_FFT = (std::complex<double> **) fftw_malloc(sizeof(std::complex<double>*) * ports_number);
	signal_IFFT = (std::complex<double> **) fftw_malloc(sizeof(std::complex<double>*) * ports_number);
  signal_CENTER = (double **) fftw_malloc(sizeof(double*) * ports_number);
  signal_WINDOW = (double **) fftw_malloc(sizeof(double*) * ports_number);
  signal_NORM = (double *) fftw_malloc(sizeof(double*) * ports_number);

  signal_CCVPHATFFT = (std::complex<double> *) fftw_malloc(sizeof(std::complex<double>) * sub_buffer_size);
	signal_CONJPROD = (std::complex<double> *) fftw_malloc(sizeof(std::complex<double>) * sub_buffer_size);
	signal_ABSCONJPROD = (std::complex<double> *) fftw_malloc(sizeof(std::complex<double>) * sub_buffer_size);
	signal_CCVIFFT = (double *) fftw_malloc(sizeof(double) * sub_buffer_size);
	resultCCV = (double *) fftw_malloc(sizeof(double) * sub_buffer_size);

    for (i = 0; i < ports_number; ++i) {

			// CCV
			signal_CONJ[i] = (std::complex<double> *) fftw_malloc(sizeof(std::complex<double>) * sub_buffer_size);
			signal_FFT[i] = (std::complex<double> *) fftw_malloc(sizeof(std::complex<double>) * sub_buffer_size);
			signal_IFFT[i] = (std::complex<double> *) fftw_malloc(sizeof(std::complex<double>) * sub_buffer_size);
			signal_CENTER[i] = (double *) fftw_malloc(sizeof(double) * sub_buffer_size);
			signal_WINDOW[i] = (double *) fftw_malloc(sizeof(double) * sub_buffer_size);

    }

	hann = (double *) calloc(sub_buffer_size, sizeof(double));
	freqs = (double *) calloc(sub_buffer_size, sizeof(double));
	freqs_energies = (double *) calloc(sub_buffer_size, sizeof(double));

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

	// Frequency indexes for DAS
	fft_size_das = 4 * nframes;
	doa_angle = (double *) malloc(2 * sizeof(double));
	doa_value = (double *) malloc(ports_number * sizeof(double));
	delay_and_sum_output_ifft = (double **) malloc(2 * sizeof(double*));
	delay_and_sum_output_ifft[0] = (double *) malloc(fft_size_das * sizeof(double));
	delay_and_sum_output_ifft[1] = (double *) malloc(fft_size_das * sizeof(double));
	freqs_das = (double *) calloc(fft_size_das, sizeof(double));

	freqs_das[0] = 0.0;
	freqs_das[fft_size_das/2] = sample_rate/2.00;
	int lim_das = fft_size_das/2;
	for (i = 1; i < lim_das; ++i) {
		freqs_das[i] = ((double)i)*(sample_rate/(double)(fft_size_das));
		freqs_das[fft_size_das-i] = -freqs_das[i];
	}

	hann_das = (double *) calloc(fft_size_das, sizeof(double));
	for (i = 0; i < fft_size_das;++i) {
		hann_das[i] = 0.5 * (1.0 - cos( 2.0 * M_PI * ((double)i)/((double)(fft_size_das - 1)) ));
		// std::cout << hann_das[i] << ",";
	}

	// Exponencials for delay
	delay_vector = (std::complex<double> **) malloc(ports_number * sizeof(std::complex<double>));
	for (i = 0; i < ports_number; ++i) {
		delay_vector[i] = (std::complex<double> *) malloc(fft_size_das * sizeof(std::complex<double>));
	}
	for (i = 0; i < fft_size_das;++i) {
		delay_vector[0][i] = std::complex<double>(1.0, 0.0);
	}
	
	delay_and_sum_fft = (std::complex<double> **) fftw_malloc(sizeof(std::complex<double>*) * 2);
	delay_and_sum_fft[0] = (std::complex<double> *) fftw_malloc(sizeof(std::complex<double>) * fft_size_das);
	delay_and_sum_fft[1] = (std::complex<double> *) fftw_malloc(sizeof(std::complex<double>) * fft_size_das);

	input_time_values = (double ***) malloc(ports_number * sizeof(double **));
	for (i = 0; i < ports_number; ++i){
		input_time_values[i] = (double **) malloc(2 * sizeof(double *));
		input_time_values[i][0] = (double *) calloc(fft_size_das, sizeof(double));
		input_time_values[i][1] = (double *) calloc(fft_size_das, sizeof(double));
	}

	window_ffts = (std::complex<double> ***) fftw_malloc(ports_number * sizeof(std::complex<double> **));
	for (i = 0; i < ports_number; ++i){
		window_ffts[i] = (std::complex<double> **) fftw_malloc(2 * sizeof(std::complex<double> *));
		window_ffts[i][0] = (std::complex<double> *) fftw_malloc(fft_size_das * sizeof(std::complex<double>));
		window_ffts[i][1] = (std::complex<double> *) fftw_malloc(fft_size_das * sizeof(std::complex<double>));
	}

	delayed_ffts = (std::complex<double> ***) fftw_malloc(ports_number * sizeof(std::complex<double> **));
	for (i = 0; i < ports_number; ++i){
		delayed_ffts[i] = (std::complex<double> **) fftw_malloc(2 * sizeof(std::complex<double> *));
		delayed_ffts[i][0] = (std::complex<double> *) fftw_malloc(fft_size_das * sizeof(std::complex<double>));
		delayed_ffts[i][1] = (std::complex<double> *) fftw_malloc(fft_size_das * sizeof(std::complex<double>));
	}

	fft1_i = 0;
	fft2_i = 2 * nframes;
	fft1_offset = (int)(2.5 * (double)nframes);
	fft2_offset = (int)(0.5 * (double)nframes);
	undelayed_offset = (int)(2.5 * (double)nframes);
	
	//preparing FFTW3 buffers
	i_fft_das = (std::complex<double> *) fftw_malloc(sizeof(std::complex<double>) * fft_size_das);
	i_time_das = (std::complex<double> *) fftw_malloc(sizeof(std::complex<double>) * fft_size_das);
	o_fft_das = (std::complex<double> *) fftw_malloc(sizeof(std::complex<double>) * fft_size_das);
	o_time_das = (std::complex<double> *) fftw_malloc(sizeof(std::complex<double>) * fft_size_das);

	i_forward_das = fftw_plan_dft_1d(fft_size_das, reinterpret_cast<fftw_complex*>(i_time_das), reinterpret_cast<fftw_complex*>(i_fft_das) , FFTW_FORWARD, FFTW_MEASURE);
	o_inverse_das = fftw_plan_dft_1d(fft_size_das, reinterpret_cast<fftw_complex*>(o_fft_das) , reinterpret_cast<fftw_complex*>(o_time_das), FFTW_BACKWARD, FFTW_MEASURE);

	buffer_size = 6 * nframes;
	for (i = 0; i < ports_number; ++i) {
		osb_array[i] = OSB<double>(buffer_size);
	}

	// MUSIC
	// printf("FREQSRES: %d\n", FREQS_RES_LEN);
	int cnt_tst = 0;
	for (i = 0; i <= lim; ++i) {
		if (freqs[i] >= FREQ_LOW_RANGE && freqs[i] <= FREQ_HIGH_RANGE){
			// printf("I: %d // CNT: %d // FREQ = %lf\n", i, cnt_tst, freqs[i]);
			cnt_tst++;
		}
	}
	int mod_array_to_select = cnt_tst / FREQS_RES_LEN;
	cnt_tst = 0;
	for (i = 0; i <= lim; ++i) {
		if (freqs[i] >= FREQ_LOW_RANGE && freqs[i] <= FREQ_HIGH_RANGE){
			if (cnt_tst == 0){
				// printf("I: %d // CNT: %d // FREQ = %lf\n", i, cnt_tst, freqs[i]);
				freqs_res_values.push_back(i);
			}
			cnt_tst++;
			if (cnt_tst >= mod_array_to_select){
				cnt_tst = 0;
			}
		}
	}
	// printf("SIZEVECT: %d\n", (int)freqs_res_values.size());

	// MUSIC
	signal_FFT_MEM = (std::complex<double> ***) fftw_malloc(sizeof(std::complex<double>**) * WINDOW_MEM_LEN);

    for (i = 0; i < WINDOW_MEM_LEN; ++i) {
		// MUSIC
		signal_FFT_MEM[i] = (std::complex<double> **) fftw_malloc(sizeof(std::complex<double>*) * ports_number);

		for (int j = 0; j < ports_number; ++j){
			// signal_FFT_MEM[i][j] = (std::complex<double> *) fftw_malloc(sizeof(std::complex<double>) * (int)freqs_res_values.size());
			signal_FFT_MEM[i][j] = (std::complex<double> *) fftw_malloc(sizeof(std::complex<double>) * sub_buffer_size);
		}
	}

	// ANGLES
	for (int j = 0; j < sub_buffer_size; ++j) {
		double curr_angle = ANGLE_LOW_RANGE;
		std::vector<std::complex<double>> current_delays2;
		std::vector<std::complex<double>> current_delays3;
		while (curr_angle <= ANGLE_HIGH_RANGE) {

			if (j == 0) {
				angle_values.push_back(curr_angle);
			}
			
			std::complex<double> freq_res_complex = std::complex<double>(freqs[j],0.0);
			current_delays2.push_back(std::exp(-COMP_I*2.0*M_PI*(freq_res_complex)*(-MIC_DIST/SOUND_SPEED)*cos((-90 - curr_angle)*MY_PI/180.0)));
			current_delays3.push_back(std::exp(-COMP_I*2.0*M_PI*(freq_res_complex)*(-MIC_DIST/SOUND_SPEED)*cos((-150.0 - curr_angle)*MY_PI/180.0)));
			curr_angle += ANGLE_STEP;

		}
		angle_delays2.push_back(current_delays2);
		angle_delays3.push_back(current_delays3);
	}
	
	// Create agent ports
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
	
		// Info for files
	char audio_file_path[100];

	audio_info = (SF_INFO*)malloc(2 * sizeof(SF_INFO));
	audio_position = (unsigned int*)malloc(2 * sizeof(unsigned int));
	audio_position[0] = 0;
	audio_position[1] = 0;
	audio_file = (SNDFILE**)malloc(2 * sizeof(SNDFILE*));

	for (int num_in_mics = 0; num_in_mics < 2; ++num_in_mics) {

		sprintf(audio_file_path, "./out_file_doa_%d.wav", num_in_mics+1);
		printf("Trying to open audio File: %s\n", audio_file_path);
		audio_info[num_in_mics].samplerate = sample_rate;
		audio_info[num_in_mics].channels = 1;
		audio_info[num_in_mics].format = SF_FORMAT_WAV | SF_FORMAT_PCM_32;
		
		audio_file[num_in_mics] = sf_open (audio_file_path,SFM_WRITE,&audio_info[num_in_mics]);
		if(audio_file[num_in_mics] == NULL){
			printf("%s\n",sf_strerror(NULL));
			exit(1);
		}else{
			printf("Audio file info:\n");
			printf("\tSample Rate: %d\n",audio_info[num_in_mics].samplerate);
			printf("\tChannels: %d\n",audio_info[num_in_mics].channels);
		}

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

    // for (i = 0; i < ports_number; ++i) {
    for (i = 0; i < OUT_SYSTEM_PORTS; ++i) {

        // Connect the first available to our output port
        if (jack_connect (client, jack_port_name (output_port[i + 1]), serverports_names[i])) {
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
