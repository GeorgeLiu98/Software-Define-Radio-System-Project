/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"
#include <math.h>       /* atan */
#include <cmath>

// function to compute the impulse response "h" based on the sinc function
/*void impulseResponseLPF(float Fs, float Fc, unsigned short int num_taps, std::vector<float> &h)
{
	// bring your own functionality
} */

//--------------------------------------------------------------
//function to generate the low pass filter coeff
void impulseResponseLPF(float Fs, float Fc, unsigned short int num_taps, std::vector<float> &h)
{
	// allocate memory and intital coeff array to all float 0
	//h.resize(num_taps,0.0);
	h.resize(0.0);

  float Norm_Fc=Fc/(Fs/2);

	for (auto i=0;i<num_taps;i++){
		if(i==((num_taps-1)/2)){
			h.push_back(Norm_Fc);
		}
		else{
			float push_result=Norm_Fc*std::sin(PI*Norm_Fc*(i-(num_taps-1)/2))/(PI*Norm_Fc*(i-(num_taps-1)/2));
			h.push_back(push_result);

		}
		h[i]*=pow(std::sin((i*PI)/(num_taps)),2.0);


	}
}

//--------------------------------------------------------------
//function for block processing
void own_block_proc_func(const std::vector<float> &block_data,\
												const std::vector<float> &h,std::vector<float> &saved_input,\
												std::vector<float> &filter_sig){

			std::vector<float> block_filtered;
			block_filtered.resize(block_data.size(),0.0);

			auto Num_taps=h.size();
			//std::cout << "Num_taps"<< Num_taps<<"\n";
			auto input_length=block_data.size();
			auto index_save=Num_taps-2;
			auto temp=0.0;

			for (auto i=0;i<input_length;i++){

				if ((i<=Num_taps-2) && (i!=0)){
					for (auto k=1;k<saved_input.size();k++){
						saved_input[k-1]=saved_input[k];

					}
					saved_input[Num_taps-2]=block_data[i-1];

				}
				for (auto j=0;j<Num_taps;j++){
						if (i<=Num_taps-2){
							if (j==0){
									block_filtered[i]=block_filtered[i]+h[j]*block_data[i];

							}
							else{
									block_filtered[i]=block_filtered[i]+h[j]*saved_input[index_save];

									index_save=index_save-1;

							}

						}
						else if (i-j>=0){
							block_filtered[i]=block_filtered[i]+h[j]*block_data[i-j];

						}

				}
				index_save=Num_taps-2;
				filter_sig.push_back(block_filtered[i]);

			}

			for (auto l=0;l<saved_input.size();l++){
				saved_input[saved_input.size()-l-1]=block_data[block_data.size()-l-1];

			}
}
//--------------------------------------------------------------
//function for block processing with downsampling capability
//this block processing only performs multiplication to every nth sample of the input block
void block_proc_downSamp(const std::vector<float> &block_data,\
												const std::vector<float> &h,std::vector<float> &saved_input,\
												std::vector<float> &filter_sig, int dw_sample){

			//std::cerr << "enter block_proc_downSamp" << std::endl;
			//std::cerr << "len of block_data="<<block_data.size() << std::endl;

			//--------------------------------------------------
			filter_sig.resize(block_data.size()/dw_sample);
			//------------------------------------------------
			//std::vector<float> block_filtered;
			//block_filtered.resize(block_data.size()/dw_sample);
			//block_filtered.resize((block_data.size()/dw_sample),0.0);
			float temp_filtered=0.0;

			auto Num_taps=h.size();
			auto input_length=block_data.size();
			auto index_save=Num_taps-2;
			auto temp=0.0;
			int next=0;


			for (auto i=0;i<input_length;i++){
				//block_filtered[i]=0.0;  //intital i_th filtered data to be 0
				temp_filtered=0.0;
				if ((i<=Num_taps-2) && (i!=0)){
						for (auto k=1;k<saved_input.size();k++){
							saved_input[k-1]=saved_input[k];

						}
						saved_input[Num_taps-2]=block_data[i-1];

				}

				//perform multiplication for every nth element
				if (i==next){

					for (auto j=0;j<Num_taps;j++){
							if (i<=Num_taps-2){
									if (j==0){

											//block_filtered[i/dw_sample]=block_filtered[i/dw_sample]+h[j]*block_data[i];
											temp_filtered=temp_filtered+h[j]*block_data[i];

									}
									else{
											//block_filtered[i/dw_sample]=block_filtered[i/dw_sample]+h[j]*saved_input[index_save];
											temp_filtered=temp_filtered+h[j]*saved_input[index_save];

											index_save=index_save-1;

									}

							}
							else if (i-j>=0){
									//block_filtered[i/dw_sample]=block_filtered[i/dw_sample]+h[j]*block_data[i-j];
									temp_filtered=temp_filtered+h[j]*block_data[i-j];

							}

					}

					index_save=Num_taps-2;
					//filter_sig.push_back(block_filtered[i]);
					//---------------------------------
					//filter_sig[i/dw_sample]=block_filtered[i/dw_sample];
					filter_sig[i/dw_sample]=temp_filtered;
					//---------------------------------
					/*std::cerr << "temp_filtered="<<temp_filtered<< std::endl;
					std::cerr << "filter_sig[i/dw_sample]="<<i/dw_sample<<"   "<<filter_sig[i/dw_sample] << std::endl;
					std::cerr << "-----------------------"<< std::endl;*/

					next=next+dw_sample;

			 	}

			}

			//update saved input for next block iteration
			for (auto l=0;l<saved_input.size();l++){
				saved_input[saved_input.size()-l-1]=block_data[block_data.size()-l-1];

			}

}

//--------------------------------------------------------------
//function for demodulation
//please go to https://www.embedded.com/dsp-tricks-frequency-demodulation-algorithms/ equation 13-117 for
//the equation used in this demnod function
void self_demod(const std::vector<float> &I,\
												const std::vector<float> &Q,std::vector<float> & prev_Q,std::vector<float> & prev_I,\
                        int block_count,std::vector<float> &fm_demod){

      float deriv_I;
      float deriv_Q;
      fm_demod.resize(0.0); //initial it back to empty vecotor with type float every time the function is called for every block

      // if part is for when block_count=0 and it prevents division zero error
      //the first demod result is set to be zero when block_count=0, for loop starts at position 1
      if (block_count==0){

        fm_demod.push_back(0.0);

        for (int k=1;k<I.size();k++){
          deriv_I=I[k]-prev_I[0];
          deriv_Q=Q[k]-prev_Q[0];

          float temp=(I[k]*deriv_Q-Q[k]*deriv_I)/(pow(I[k],2.0)+pow(Q[k],2.0));
          fm_demod.push_back(temp);

          prev_I[0]=I[k];
          prev_Q[0]=Q[k];

        }
      }
      //else statement is for all other blocks except when block_count=0
      //for loop starts at postion 0
      else {

        for (int k=0;k<I.size();k++){
          deriv_I=I[k]-prev_I[0];
          deriv_Q=Q[k]-prev_Q[0];

          float temp=(I[k]*deriv_Q-Q[k]*deriv_I)/(pow(I[k],2.0)+pow(Q[k],2.0));
          fm_demod.push_back(temp);

          prev_I[0]=I[k];
          prev_Q[0]=Q[k];

        }

      }
}
//--------------------------------------------------------------
//function for upsampling and expanding
void expender(const std::vector<float> &block_in,\
								std::vector<float> &block_out,\
								int up_samp_factor){

		block_out.resize(up_samp_factor*block_in.size(),0.0);

		for (int i=0;i<block_out.size() ;i++){

				//check if the position in the expanded array is a multiple of the up_samp_factor
				if (i%up_samp_factor==0){
						block_out[i]=block_in[i/up_samp_factor];
				}
		}
}

//--------------------------------------------------------------
//function to generate the bandpass filter coeff (Fb: begining of passband, Fc: end of passband)
void impulseResponseBPF(float Fs, float Fb, float Fe, unsigned short int num_taps, std::vector<float> &h)
{
	// allocate memory and intital coeff array to all float 0
	//h.resize(num_taps, 0.0);
	h.resize(0.0);

	float Norm_center=((Fb+Fe)/2)/(Fs/2);
	float Norm_pass=(Fe-Fb)/(Fs/2);

	for (auto i=0;i<num_taps;i++){
		if(i==((num_taps-1)/2)){
			h.push_back(Norm_pass);

		}
		else{
			float push_result = Norm_pass*sin(PI*(Norm_pass/2)*(i-(num_taps-1)/2))/(PI*(Norm_pass/2)*(i-(num_taps-1)/2));
			std::cout << "i="<<i << '\n';
			std::cout << "push_back="<<push_result << '\n';

			h.push_back(push_result);

		}

		h[i]=h[i]*cos(i*PI*Norm_center);
		h[i]=h[i]*pow(sin((i*PI)/(num_taps)),2.0);
		std::cout << "h[i]="<<h[i]<< '\n';
		std::cout << "---------------------------" << '\n';

	}
}

void PLLNCO (std::vector<float> &pllIn, std::vector<float> &ncoOut, float freq, float Fs,\
	 						std::vector<float> state, float ncoScale = 2.0, float phaseAdjust = 0.0, float normBandwith = 0.01)
{
		/*

		pllIn 	 				array of floats
										input signal to the PLL (assume known frequency)

		freq 						float
										reference frequency to which the PLL locks

		Fs  						float
										sampling rate for the input/output signals

		ncoScale				float
										frequency scale factor for the NCO output

		phaseAdjust			float
										phase adjust to be added to the NCO only

		normBandwidth		float
										normalized bandwidth for the loop filter
										(relative to the sampling rate)

		ncoOut 					array of floats
										output after the PLL and NCO

		state

		*/

		float Cp = 2.666;
		float Ci = 3.555;

		float Kp = normBandwith*Cp;
		float Ki = normBandwith*normBandwith*Ci;

#define PI 3.14159265358979323846
		ncoOut.resize(pllIn.size()+1);

		float integrator = state[0];
		float phaseEst = state[1];
		float feedbackI = state[2];
		float feedbackQ = state[3];
		ncoOut[0] = state[4];
		float trigOffset = state[5];

		for (auto k = 0;k<pllIn.size();k++){

				float errorI=pllIn[k]*(+feedbackI);
				float errorQ=pllIn[k]*(-feedbackQ);

				float errorD = atan2(errorQ, errorI);

				integrator=integrator+Ki*errorD;
				phaseEst=phaseEst+Kp*errorD+integrator;

				float trigArg = 2*PI*(freq/Fs)*(trigOffset+k+1) + phaseEst;
				feedbackI = cos(trigArg);
				feedbackQ = sin(trigArg);
				ncoOut[k+1] = cos(trigArg*ncoScale + phaseAdjust);
		}
}






//--------------------------------------------------------------
// function to compute the filtered output "y" by doing the convolution
// of the input data "x" with the impulse response "h"
/*void convolveFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h)
{
	// bring your own functionality
}*/

/*void convolveFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h)
{
	// allocate memory for the output (filtered) data
	//y.resize(x.size()+h.size()-1, 0.0);
	y.resize(0.0);
	auto return_size=h.size()+x.size()-1;
	std::cout << "size of x=" << x.size() << "\n";
	std::cout << "size of h= " << h.size() << "\n";

	for (auto n=0;n<return_size;n++){
		float temp=0.0;

		for (auto k=0;k<h.size();k++){
			if ((n-k)>=0 && (n-k)<=return_size-h.size()){
				temp+=x[n-k]*h[k];
			}

		}
		y.push_back(temp);
	}
	std::cout << "size of return_size=" << return_size << "\n";


	// the rest of the code in this function is to be completed by you
	// based on your understanding and the Python code from the first lab
}*/

void mixing(std::vector<float> &nco_out, std::vector<float> &stereo_filt, std::vector<float> &rtn)	{
	//Pointwise multiplication
	//rtn & its size to be declared and passed in from experiment.cpp
	//a is the NCO_output
	//b is the stereo filtered data
	for(int i = 0; i < stereo_filt.size(); i++) {
		rtn[i] = (stereo_filt[i] * nco_out[i+1]);
	}
}

void stereoLeftRight(std::vector<float> &mono_data, std::vector<float> &stereo_data,\
	std::vector<float> &stereo_left, std::vector<float> &stereo_right ) {
		//This function will proceses and generated a block of stereo right and left
		//Set the size of stereo_left and stereo_right to equal the size of mono_data in experiment.cpp
		for(int i = 0; i < mono_data.size(); i++) {
				stereo_left[i] = ((mono_data[i] + stereo_data[i])/2);
				stereo_right[i] = ((mono_data[i] - stereo_data[i])/2);
		}
	}

void interleave(std::vector<float> &left, std::vector<float> &right, std::vector<float> &rtn) {
	//Set the size of rtn is equal to the size of left + right
	//combine left and right (two arrays into one)
	auto j=0;

	for(int i = 0; i < left.size(); i++) {
		rtn[j] = left[i];
		rtn[j+1] = right[i];
		j=j+2;
	}
}
