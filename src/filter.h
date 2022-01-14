/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_FILTER_H
#define DY4_FILTER_H

// add headers as needed
#include <iostream>
#include <vector>

// declaration of a function prototypes
void impulseResponseLPF(float, float, unsigned short int, std::vector<float> &);

void own_block_proc_func(const std::vector<float> &,\
												const std::vector<float> &,std::vector<float> &,\
												std::vector<float> &);

void self_demod(const std::vector<float> &,\
								const std::vector<float> &,\
								std::vector<float> &, std::vector<float> &, int,\
								std::vector<float> &);

void block_proc_downSamp(const std::vector<float> &block_data,\
											const std::vector<float> &h,std::vector<float> &saved_input,\
											std::vector<float> &filter_sig, int dw_sample);

void expender(const std::vector<float> &block_in,\
								std::vector<float> &block_out,\
								int up_samp_factor);

void convolveFIR(std::vector<float> &, const std::vector<float> &, const std::vector<float> &);

void impulseResponseBPF(float, float,float, unsigned short int, std::vector<float> &);

void PLLNCO (std::vector<float> &, std::vector<float> &, float, float,\
	 						std::vector<float> , float, float , float );

void mixing(std::vector<float> &, std::vector<float> &, std::vector<float> &);

void stereoLeftRight(std::vector<float> &, std::vector<float> &,\
	std::vector<float> &, std::vector<float> &);

void interleave(std::vector<float> &, std::vector<float> &, std::vector<float> &);

#endif // DY4_FILTER_H
