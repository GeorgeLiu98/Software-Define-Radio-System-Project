/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Copyright by Nicola Nicolici
Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>

#define queue_blocks 5


void readStdinBlcokData(unsigned int num_samples, \
												unsigned int block_id,\
												std::vector<float> &block_data)
{
	std::vector<char> raw_data(num_samples);
	std::cin.read(reinterpret_cast<char*>(&raw_data[0]), num_samples*sizeof(char));

	for(unsigned int k = 0; k < num_samples; k++){
		block_data[k] = float(((unsigned char)raw_data[k]-128)/128.0);
	}
}

void IQ_class(const std::vector<float> &block_data,\
	 			std::vector<float> &Is,\
				std::vector<float> &Qs){
	for(int i = 0; i < block_data.size(); i++){
		if(i%2==0){
			Is.push_back(block_data[i]);
		}
		else{
			Qs.push_back(block_data[i]);
		}
	}
}
//---------------------------------------
//thread for RF front end
void RF_front_end(std::queue<std::vector <float>> &my_queue,std::mutex &my_mutex,\
									std::condition_variable &my_cvar,std::vector<float> &rf_coeff,int block_size,\
								  std::vector<float> &saved_I,std::vector<float> &saved_Q,std::vector<float> &filter_sig_I,\
									std::vector<float> &filter_sig_Q,std::vector<float> &prev_I,std::vector<float> &prev_Q,\
									std::vector<float> &fm_demod,std::vector<float> &Is,std::vector<float> &Qs,int rf_decim,int mode){

	//for reading data and performing calculation
	for (unsigned int block_id = 0; ; block_id++){
			//read block data in
			std::cerr << "Read block " << block_id << std::endl;
			std::vector<float> block_data(block_size);
			readStdinBlcokData(block_size, block_id, block_data);

			if((std::cin.rdstate())!=0){
					std::cerr << "--------------end of mode "<< mode<<" calculation and streaming------------------"<< std::endl;
					std::cerr << "End of input stream reached" << '\n';
					//break;
					exit(1);
			}
			//function to seperate I/Q data samples and initialize to empty vectors
			Is.resize(0.0);
			Qs.resize(0.0);
			filter_sig_I.resize(0.0);
			filter_sig_Q.resize(0.0);

			//function to obtain I/Q samples
			IQ_class(block_data, Is, Qs);
			if (mode==0){
				//RF front end for mode 1

				//pass RF front end filter, perfrom both filtering and downSampling for I/Q
				block_proc_downSamp(Is,rf_coeff,saved_I,filter_sig_I,rf_decim);
				block_proc_downSamp(Qs,rf_coeff,saved_Q,filter_sig_Q,rf_decim);

				//perform demodulation
				self_demod(filter_sig_I,filter_sig_Q,prev_Q,prev_I,block_id,fm_demod);

			}

			else if (mode==1){
				//RF front end for mode 1
				//perfrom both filtering and downSampling for I/Q
				block_proc_downSamp(Is,rf_coeff,saved_I,filter_sig_I,rf_decim);
				block_proc_downSamp(Qs,rf_coeff,saved_Q,filter_sig_Q,rf_decim);

				//perform demodulation
				self_demod(filter_sig_I,filter_sig_Q,prev_Q,prev_I,block_id,fm_demod);


			}
				//push into the queue after RF front end calculation is done
				std::unique_lock<std::mutex>my_lock(my_mutex);
				if (my_queue.size()==queue_blocks){
					my_cvar.wait(my_lock);
				}

				my_queue.push(fm_demod);

				my_lock.unlock();
				my_cvar.notify_one();


		}

}
//---------------------------------------
//thread for audio
void audio_compute(std::queue <std::vector<float>> &my_queue,std::mutex &my_mutex,\
									std::condition_variable &my_cvar,std::vector<float> &audio_coeff,std::vector<float> &audio_save,\
									int audio_decim,std::vector<float> &audio_filt,std::vector<short int> &audio_data,int mode){

		std::vector<float> block_up_samp;
		while(true){
			//std::cerr << "enter audio_compute" << '\n';

			//poping from the queue
			std::unique_lock<std::mutex>my_lock(my_mutex);
			if (my_queue.empty()){
				my_cvar.wait(my_lock);

			}
			//std::cerr << "pop testing" << '\n';
			std::vector<float> read_demod=my_queue.front();
			my_queue.pop();


			audio_filt.resize(0.0);
			audio_data.resize(0.0);
			//block_up_samp.resize(0.0);

			//--------------------------------------------------------------------
			if (mode==0){
			//perfrom filtering for audio
				block_proc_downSamp(read_demod,audio_coeff,audio_save,audio_filt,audio_decim);
			}

			else if (mode==1){
				int up_samp_factor=24;
				//std::cerr << "enter mode 1 compute_audio" << '\n';

				//upsampling by a factor of 24
				expender(read_demod,block_up_samp,up_samp_factor);

				//perfrom both filtering and downSampling for the demodulated signal
				block_proc_downSamp(block_up_samp,audio_coeff,audio_save,audio_filt,audio_decim);

			}



			//--------------------------------------------------------------------
			//set a magnify factor for mode 1
			int mag;
			if (mode==0)mag=1;
			else if (mode==1) mag=100;
			//--------------------------------------------------------------------

			//streaming output audio
			for(unsigned int k = 0; k < audio_filt.size(); k++){
				//if(std::isnan(block_data[k])) audio_data[k] = 0;
				//else
				//{
					audio_data.push_back(static_cast<short int>(audio_filt[k]*16384*mag));
				//}
			}

			//output audio block by block
			fwrite(&audio_data[0], sizeof(short int), audio_data.size(), stdout);

			my_lock.unlock();
			my_cvar.notify_one();
	}
}



//--------------------------------------------------
int main(int argc, char* argv[])
{
	int mode = 0;

	//--------------------------------------------------
	//declare mode 0 parameters

	float rf_Fs=2.4e6;
	float rf_Fc=100e3;
	int rf_taps=151;
	int rf_decim=10;

	float audio_Fc=16e3;
	int audio_taps=151;
	int audio_decim=5;

	//--------------------------------------------------
	//declare mode 1 parameters
	float rf_Fs1=2.5e6;
	float rf_Fc1=100e3;
	int rf_taps1=151;
	int rf_decim1=10;
	int up_samp_factor=24;

	float audio_Fc1=16e3;
	int audio_taps1=151;
	int audio_decim1=125;

	float m1_Fs = 250e3;

	//--------------------------------------------------
	//declare stereo parameters

	//cutoff frequency for stereo carrier recovery
/*	float pilot_low = 18.5e3;
	float pilot_high = 19.5e3;

	//Cutoff frequency for stereo audio extraction
	float stereo_low = 22e3;
	float stereo_high = 54e3

	float pllFreq=19e3

	// nco States inital values
	std::vector<float> state_nco={0.0,0.0,1.0,0.0,1.0,0.0};*/

	//--------------------------------------------------



	int block_size;

	if(argc < 2){
		std::cerr << "Operating in default mode 0" << std::endl;
	}
	else if(argc == 2){
		mode = atoi(argv[1]);
		if(mode != 1){
			std::cerr << "Wrong mode" <<mode<< std::endl;
			exit(1);
		}
	}
	else{
				std::cerr << "Usage: " << argv[0] << std::endl;
				std::cerr << "or " << std::endl;
				std::cerr << "Usage: " << argv[0] << " 1" << std::endl;
				exit(1);
	}
	//obtain coeffs for two filters based on either it is mode 1 or 0
	std::vector<float> rf_coeff;
	// mono
	std::vector<float> audio_coeff;
	// stereo
	std::vector<float> carrier_coeff;
	std::vector<float> streo_coeff;

	if (mode ==0){
			std::cerr << "--------- mode "<< mode<<" filter coeffs---------------"<< std::endl;

			block_size = 1024*rf_decim*audio_decim*2;
			// RF front end filter
			impulseResponseLPF(rf_Fs, rf_Fc,rf_taps, rf_coeff);
			// Mono data filter
			impulseResponseLPF(rf_Fs/rf_decim, audio_Fc,audio_taps, audio_coeff);

			//---------------------------------------------------------------------
			// Stereo carrier recovery BPF
			//impulseResponseBPF(rf_Fs, pilot_low, pilot_high, audio_taps, carrier_coeff);
			// Stereo Channel Extraction BPF
			//impulseResponseBPF(rf_Fs, stereo_low, stereo_high, audio_taps, stereo_coeff);

	}
	else if (mode==1){
			std::cerr << "--------- mode "<< mode<<" filter coeffs---------------"<< std::endl;

			block_size = rf_decim1*audio_decim1*up_samp_factor*2;
			impulseResponseLPF(rf_Fs1, rf_Fc1,rf_taps1, rf_coeff);
			impulseResponseLPF(rf_Fs1/rf_decim1*up_samp_factor, audio_Fc1,audio_taps1, audio_coeff);

	}

	//set up vecotrs for demod function
	std::vector<float>  prev_Q;
	std::vector<float>  prev_I;

	prev_I.resize(1);
	prev_Q.resize(1);
	prev_I[0]=0.0;
	prev_Q[0]=0.0;

	//vectors for I and Q samples
	std::vector<float> Is;
	std::vector<float> Qs;


	//two vectors to store saved inputs and initalize them two zeros
	std::vector<float> saved_I;
	saved_I.resize(rf_taps-1,0.0);
	std::vector<float> saved_Q;
	saved_Q.resize(rf_taps-1,0.0);

	//declare two vecotrs to store the filtered signal for both I/Q
	std::vector<float> filter_sig_I;
	filter_sig_I.resize(0.0);
	std::vector<float> filter_sig_Q;
	filter_sig_Q.resize(0.0);

	//declare a vector to store the demodulated signal
	std::vector<float> fm_demod;
	fm_demod.resize(0.0);

	//declare two vecotrs to store values for audio_save and audio_filt
	std::vector<float> audio_save;
	audio_save.resize(audio_taps-1,0.0);
	std::vector<float> audio_filt;
	audio_filt.resize(0.0);

	std::vector<float> stereo_save;
	stereo_save.resize(audio_taps-1,0.0);
	std::vector<float> stereo_filt;
	stereo_filt.resize(0.0);

	//vector to store value for output audio data
	std::vector<short int> audio_data;
	audio_data.resize(0.0);

	//for mode1 block data upsampling
	std::vector<float> block_up_samp;
	//------------------------------------------------------------
	//threading included

	//declare queue, mutex and condition_variable for threading
	std::queue<std::vector<float>>my_queue;
	std::mutex my_mutex;
	std::condition_variable my_cvar;

	/*if (mode==0)audio_decim=5;
	else */if (mode==1)audio_decim=audio_decim1;

	std::thread tp=std::thread(RF_front_end,std::ref(my_queue),std::ref(my_mutex),std::ref(my_cvar),std::ref(rf_coeff),block_size,\
														std::ref(saved_I),std::ref(saved_Q),std::ref(filter_sig_I),std::ref(filter_sig_Q),\
														std::ref(prev_I),std::ref(prev_Q),std::ref(fm_demod),std::ref(Is),std::ref(Qs),rf_decim,mode);


	std::thread tc=std::thread(audio_compute,std::ref(my_queue),std::ref(my_mutex),std::ref(my_cvar),std::ref(audio_coeff),\
														std::ref(audio_save),audio_decim,std::ref(audio_filt),std::ref(audio_data),mode);


	tp.join();
	tc.join();


	//------------------------------------------------------------
	//for reading data and performing calculation
	/*for (unsigned int block_id = 0; ; block_id++){
		std::vector<float> block_data(block_size);
		readStdinBlcokData(block_size, block_id, block_data);

		if((std::cin.rdstate())!=0){
			std::cerr << "--------------end of mode "<< mode<<"calculation and streaming------------------"<< std::endl;
			std::cerr << "End of input stream reached" << '\n';
			break;
		}
		std::cerr << "Read block " << block_id << std::endl;

		//function to seperate I/Q data samples and initialize to empty vectors
		Is.resize(0.0);
		Qs.resize(0.0);
		IQ_class(block_data, Is, Qs);

		//------------------------------------------------------------------------
		//mode 0 calculation
		if (mode==0){
			std::cerr << "-----------------enter mode 0 calculation-----------------------" <<  std::endl;
			//--------------------------RF Front End---------------------------------

			//perfrom both filtering and downSampling for I/Q
			block_proc_downSamp(Is,rf_coeff,saved_I,filter_sig_I,rf_decim);

			block_proc_downSamp(Qs,rf_coeff,saved_Q,filter_sig_Q,rf_decim);

			//perform demodulation
			self_demod(filter_sig_I,filter_sig_Q,prev_Q,prev_I,block_id,fm_demod);

			//---------------------------Mono processing-------------------------------

			//perfrom both filtering and downSampling for the demodulated signal
			block_proc_downSamp(fm_demod,audio_coeff,audio_save,audio_filt,audio_decim);

			//-------------------------Stereo processing-------------------------------
			// Stereo channel extraction filtering
			// putting in 1 means no downsampling
			//block_proc_downSamp(fm_demod,stereo_coeff,stereo_filt, , 1);

	}
	//mode 1 calculation
	else if (mode ==1){
			std::cerr << "-----------------enter mode 1 calculation-----------------------" <<  std::endl;
			//--------------------------RF Front End---------------------------------

			//perfrom both filtering and downSampling for I/Q
			block_proc_downSamp(Is,rf_coeff,saved_I,filter_sig_I,rf_decim1);

			block_proc_downSamp(Qs,rf_coeff,saved_Q,filter_sig_Q,rf_decim1);

			//perform demodulation
			self_demod(filter_sig_I,filter_sig_Q,prev_Q,prev_I,block_id,fm_demod);

			//---------------------------Mono processing-------------------------------

			//upsampling by a factor of 24
			expender(fm_demod,block_up_samp,up_samp_factor);

			//perfrom both filtering and downSampling for the demodulated signal
			block_proc_downSamp(block_up_samp,audio_coeff,audio_save,audio_filt,audio_decim1);

	}

		//--------------------------------------------------------------------
		//set a magnify factor for mode 1
		int mag;
		if (mode==0)mag=1;
		else if (mode==1) mag=100;
		//--------------------------------------------------------------------

		for(unsigned int k = 0; k < audio_filt.size(); k++){
			if(std::isnan(block_data[k])) audio_data[k] = 0;
			else
			{

				audio_data.push_back(static_cast<short int>(audio_filt[k]*16384*mag));

			}
		}

		//[L,R,L,R,L,R,L,R]????

		//output audio block by block
		fwrite(&audio_data[0], sizeof(short int), audio_data.size(), stdout);
		// fwrite(samples, 1, numSamples*CHANNELS, outFile);

		//reset vecotors that store intermediate values
		block_up_samp.resize(0.0);
		audio_data.resize(0.0);
		filter_sig_I.resize(0.0);
		filter_sig_Q.resize(0.0);
		audio_filt.resize(0.0);



	}*/
	//finish all blocks and output the audio file as one
	//fwrite(&audio_data[0], sizeof(short int), audio_data.size(), stdout);
	return 0;
}



//-------------------------------------------------------------------------------------------
//main function given (default code of lab3)

//int main()
//{
	// binary files can be generated through the
	// Python models from the "../model/" sub-folder
	/*const std::string in_fname = "../data/fm_demod_10.bin";
	std::vector<float> bin_data;
	readBinData(in_fname, bin_data);

	// generate an index vector to be used by logVector on the X axis
	std::vector<float> vector_index;
	genIndexVector(vector_index, bin_data.size());
	// log time data in the "../data/" subfolder in a file with the following name
	// note: .dat suffix will be added to the log file in the logVector function
	logVector("demod_time", vector_index, bin_data);

	// take a slice of data with a limited number of samples for the Fourier transform
	// note: NFFT constant is actually just the number of points for the
	// Fourier transform - there is no FFT implementation ... yet
	// unless you wish to wait for a very long time, keep NFFT at 1024 or below
	std::vector<float> slice_data = \
		std::vector<float>(bin_data.begin(), bin_data.begin() + NFFT);
	// note: make sure that binary data vector is big enough to take the slice

	// declare a vector of complex values for DFT
  std::vector<std::complex<float>> Xf;
	// ... in-lab ...
	// compute the Fourier transform
	// the function is already provided in fourier.cpp

	// compute the magnitude of each frequency bin
	// note: we are concerned only with the magnitude of the frequency bin
	// (there is NO logging of the phase response, at least not at this time)
	std::vector<float> Xmag;
	// ... in-lab ...
	// compute the magnitude of each frequency bin
	// the function is already provided in fourier.cpp

	// log the frequency magnitude vector
	vector_index.clear();
	genIndexVector(vector_index, Xmag.size());
	logVector("demod_freq", vector_index, Xmag);   */ // log only positive freq

	// for your take-home exercise - repeat the above after implementing
	// your OWN function for PSD based on the Python code that has been provided
	// note the estimate PSD function should use the entire block of "bin_data"
	//
	// ... complete as part of the take-home ...
	//
	//--------------------------------------------------------------
	//testing filter coeff function
	//the filter coeff can be checked use the python function we implemented in lab 1 and C++ in lab2
	//go to filter.cpp for impulseResponseLPF function which generates the filter coeffs

	/*std::vector<float> h;

	impulseResponseLPF(2.4e6/10, 16e3, 151, h);
	for (int i=0;i<151;i++){
		std::cout << "h coeff(0 to 150)=" << h[i] << "\n";

	}*/
	//--------------------------------------------------------------
	//testing block processing function
	//we can compare the block processing resutls with regular convolution calculator online
	//go to filter.cpp to see the own_block_proc_func function

/*	std::vector<float> h{3,4,5,6,7,8};
	int num_taps=6;
	std::vector<float> block_data{45,42.54,534.54,342,5345,832,54,23.32};
	std::vector<float> saved_input;
	saved_input.resize(num_taps-1,0.0);
	std::vector<float> filter_sig;
	filter_sig.resize(0.0);


	own_block_proc_func(block_data,h,saved_input,filter_sig);

	std::vector<float> block_data2 {34,655,43,45,6,54352,324.43,352.5};

	own_block_proc_func(block_data2,h,saved_input,filter_sig);

	for (int i=0;i<filter_sig.size();i++){
		std::cout << "filtered sig=" <<filter_sig[i] << "\n";

	}*/

	//--------------------------------------------------------------
	//testing for demod function, comparing the results with lab 3 python implementation
	//this test uses two different I/Q data sets to simulate the behaviour of two blocks of data: block 0 and 1
	//I/Q data can be changed for more test cases

	//1. prev_I and prev_Q need to be declared as vecotrs with length 1 in top-level function and initial them to be zero
	//for block 0 (very first block) and two values need to be updated in the function for the next block demod calculation
	//2. vector fm_demod needs to be declared in the top-level. self_demod function sets it back to empty vector at the beginning
	//of each time when the function is called
	//3. go to filter.cpp to see self_demod function
	/*std::vector<float> fm_demod;

	//initalize two vectors with length of one each to store values for prev_I and prev_Q
	std::vector<float> prev_I;
	std::vector<float> prev_Q;
	prev_I.resize(1);
	prev_Q.resize(1);
	prev_I[0]=0.0;
	prev_Q[0]=0.0;

	std::vector<float> I{3,4,5,6,7,238,239,12,432,42,5,3};
	std::vector<float> Q{12,32,54,643,732,432,232,13,23,12,31,314};
	int block_count=0;

	self_demod(I,Q,prev_Q,prev_I,block_count,fm_demod);

	std::cout << "----------------first fm_demod-----------------" << "\n";
  for (int i=0;i<fm_demod.size();i++){
			std::cout << "demod function result=" <<fm_demod[i] << "\n";
	}
	block_count=1;

	I={13,34,25,36,47,28,39,132,222,452,23,376};
	Q={232,332,534,6243,7424,4322,2332,133,232,132,313,313};

	std::cout << "prev_I=" <<prev_I[0] << "\n";
	std::cout << "prev_Q=" <<prev_Q[0] << "\n";


	self_demod(I,Q,prev_Q,prev_I,block_count,fm_demod);

	std::cout << "----------------second fm_demod-----------------" << "\n";
	for (int i=0;i<fm_demod.size();i++){
			std::cout << "demod function result=" <<fm_demod[i] << "\n";
	}
	std::cout << "prev_I=" <<prev_I[0] << "\n";
	std::cout << "prev_Q=" <<prev_Q[0] << "\n"; */

	//--------------------------------------------------------------
	//testing for block processing function with downsampling capability

/*	std::vector<float> h{3,4,5,6,7,8};
	int num_taps=6;
	std::vector<float> block_data{45,42.54,534.54,342,5345,832,54,23.32,23,4,3.23,32,23.32};
	std::vector<float> saved_input;
	saved_input.resize(num_taps-1,0.0);
	std::vector<float> filter_sig;
	filter_sig.resize(0.0);
	int dw_sample=2;


	block_proc_downSamp(block_data,h,saved_input,filter_sig,dw_sample);

	for (int i=0;i<filter_sig.size();i++){
		std::cout << "filtered sig=" <<filter_sig[i] << "\n";

	}

	std::cout << "---------------- block_proc with downsampling -----------------" << "\n";
	std::vector<float> block_data2 {34,655,43,45,6,54352,324.43,352.5,23,4,3.23,32,23.32};

	filter_sig.resize(0.0);
	block_proc_downSamp(block_data2,h,saved_input,filter_sig,dw_sample);

	for (int i=0;i<filter_sig.size();i++){
		std::cout << "filtered sig=" <<filter_sig[i] << "\n";

	}*/

	//--------------------------------------------------------------
	//testing for expanding function

	/*std::vector<float> block_in{3,4,5,6,7,8,34,123,43,12,32,32,121,342};
	int up_samp_factor=24;

	std::vector<float> block_out;

	expender(block_in,block_out,up_samp_factor);

	int count=0;

	for (int i=0;i<block_out.size();i++){
		std::cout << "block_out value=" <<block_out[i] << "\n";

	}*/

	//--------------------------------------------------------------
	//testing for passband filter function

	/*float Fs=240e3;
	float Fb=22e3;
	float Fe=54e3;
	int num_taps=151;
	std::vector<float> h;

	impulseResponseBPF(Fs,Fb,Fe,num_taps,h);

	std::ofstream myfile;
	myfile.open("bpf.txt");
	for(int i = 0; i<h.size(); i++) {
		myfile << h[i] <<"\n";
		std::cout << h[i] << std::endl;
	}

  myfile.close();

	for (int i=0;i<h.size();i++){
		std::cout << "BPF coeff value=" <<h[i] << "\n";

	}*/

	//--------------------------------------------------------------
	//testing for PLLNCO function

	/*std::vector<float> state={0.0,0.0,1.0,0.0,1.0,0.0};

	float freq = 19e3;
	float Fs = 240e3;

	std::vector<float> pllIn= {1.0,2.0,3.0,4.0,5.0}; // test input
	std::vector<float> ncoOut;


	PLLNCO(pllIn, ncoOut, freq, Fs, state, 2.0, 0.0, 0.01);

	for (int i=0;i<ncoOut.size();i++){
						std::cout << "output nco value=" <<ncoOut[i] << "\n";

	}*/
	//--------------------------------------------------------------
	//testing for mixing function
	/*std::vector<float> nco_out={33,23,32.5,34,54,34.5,34.43,1222.32,23.34};
	std::vector<float> stereo_filt={56.7,546.5,32.565,34.453,54.534,3453.5345,3534.43,1543.32};

	std::vector<float> rtn;
	rtn.resize(8,0.0);


	mixing(nco_out,stereo_filt,rtn);

	for (int i=0;i<rtn.size();i++){
						std::cout << "rtn[i]=" <<rtn[i] << "\n";

	}*/

	//--------------------------------------------------------------
	//testing for stereoLeftRight & interleave functions
	/*std::vector<float> mono_data={23.4,234.5,45.56,65.45};
	std::vector<float> stereo_data={123.5,43.5,232.4,24.54};

	std::vector<float> stereo_left;
	stereo_left.resize(stereo_data.size(),0.0);
	std::vector<float> stereo_right;
	stereo_right.resize(stereo_data.size(),0.0);

	std::vector<float>rtn;
	rtn.resize(stereo_data.size()*2);


	stereoLeftRight(mono_data,stereo_data,stereo_left,stereo_right );
	for (int i=0;i<stereo_right.size();i++){
				std::cout << "stereo_left[i]=" <<stereo_left[i] << "\n";
				std::cout << "stereo_right[i]=" <<stereo_right[i] << "\n";
	}

	std::cout << "-------------------" << "\n";

	interleave(stereo_left, stereo_right,rtn);

	for (int i=0;i<rtn.size();i++){
				std::cout << "rtn[i]=" <<rtn[i] << "\n";

	}*/


		//--------------------------------------------------------------

	// if you wish to write some binary files, see below example
	// const std::string out_fname = "../data/outdata.bin";
	// writeBinData(out_fname, bin_data);

	// naturally, you can comment the line below once you are comfortable to run gnuplot
	//std::cout << "Run: gnuplot -e 'set terminal png size 1024,768' example.gnuplot > ../data/example.png\n";

//	return 0;
//}
