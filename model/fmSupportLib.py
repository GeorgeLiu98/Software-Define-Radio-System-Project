#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import numpy as np
import math, cmath

#
# you should use the demodulator based on arctan given below as a reference
#
# in order to implement your OWN FM demodulator without the arctan function,
# a very good and to-the-point description is given by Richard Lyons at:
#
# https://www.embedded.com/dsp-tricks-frequency-demodulation-algorithms/
#
# the demodulator boils down to implementing equation (13-117) from above, where
# the derivatives are nothing else but differences between consecutive samples
#
# needless to say, you should not jump directly to equation (13-117)
# rather try first to understand the entire thought process based on calculus
# identities, like derivative of the arctan function or derivatives of ratios

#----------------------------------------------------------------------
#function to gerate own filter coeff
def filterCoeff(Fs,Fc,N_taps):
    print ("enter filterCoeff fucntion")
    Norm_CF=Fc/(Fs/2)
    
    result=np.zeros(N_taps)
	
	#for loop for FIR filter
    for i in range (0,N_taps):
        
        if i==(N_taps-1)/2:
            
            result[i]=Norm_CF
        else:
            
            result[i]=Norm_CF*math.sin(math.pi*Norm_CF*(i-(N_taps-1)/2))/(math.pi*Norm_CF*(i-(N_taps-1)/2))

        result[i]=result[i]*(math.sin((i*math.pi)/N_taps)**2)
        
    return result

#----------------------------------------------------------------------

#----------------------------------------------------------------------
#function to perform convoltuion 
def convolution(x,h):
    print ("enter convolution fucntion")
    #prod=[]
    k=0
    finalL=len(x)+len(h)-1
    prod=np.zeros(finalL,dtype='float_')
    
    for n in range (0, finalL):
        res=0
        for k in range (0,len(h)):
            if ((n-k)>=0 and (n-k)<=finalL - len(h)):
                res=res+x[n-k]*h[k]
        #prod.append(res)
        prod[n]=res

    return prod
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#own function for demodulation

def self_demod(I,Q,prev_I=0.0, prev_Q=0.0,block_count=0):
    print ("enter self_demod fucntion")
    fm_demod = np.empty(len(I))
    
    #initiate fm_demod[0] to 0 due to division equals to zero
    if (block_count==0):
        fm_demod[0]=0.0
    
        #for k in range(1,len(I)-1):
        for k in range(1,len(I)):
            # take the derivative of I
            
            #print ("k=",k)
            deriv_I=I[k]-prev_I
            
            # take the derivative of Q
            deriv_Q=Q[k]-prev_Q
            #print ('k=',k)
            
            #denom=math.pow(I[k],2)+math.pow(Q[k],2)
            #print ("I[k]= ",I[k])
            #print ("Q[k]=",Q[k])
            #print ("denom=",denom)
            #print ("------------------------")
            fm_demod[k]=(I[k]*deriv_Q-Q[k]*deriv_I)/(math.pow(I[k],2)+math.pow(Q[k],2))
            
            prev_I=I[k]
            prev_Q=Q[k]
    else:
        for k in range(0,len(I)):
            # take the derivative of I
            deriv_I=I[k]-prev_I
            
            # take the derivative of Q
            deriv_Q=Q[k]-prev_Q
            #print ('k=',k)
            
            #denom=math.pow(I[k],2)+math.pow(Q[k],2)
            #print ("I[k]= ",I[k])
            #print ("Q[k]=",Q[k])
            #print ("denom=",denom)
            #print ("------------------------")
            fm_demod[k]=(I[k]*deriv_Q-Q[k]*deriv_I)/(math.pow(I[k],2)+math.pow(Q[k],2))
            
            prev_I=I[k]
            prev_Q=Q[k]
        
    
    print ("len fm_demod= ",len(fm_demod))
    return fm_demod,prev_I,prev_Q

#----------------------------------------------------------------------
def own_block_proc_func(input_sig,coeff,saved_input):
	#filtered_sig = np.empty(shape = input_sig.shape)
	filtered_sig=np.zeros(len(input_sig),dtype='float_')

	Num_taps=len(coeff)
	input_lenth=len(input_sig)
	#print ("begining saved input=",saved_input)
	index_save=Num_taps-2

	for i in range (input_lenth):

		if i<=Num_taps-2 and i!=0:
			#shift all elements to the left by one and inset the newest input sig value to -1 position
			saved_input[0:Num_taps-2]=saved_input[1:]
			saved_input[Num_taps-2]=input_sig[i-1]
			#print ("intermediate saved input=",saved_input)
			#print ("i value=", i)

		for j in range (Num_taps):

			if i<=Num_taps-2:
				if j==0:
					filtered_sig[i]=filtered_sig[i]+coeff[j]*input_sig[i]
					#print ("filtered [i]=",filtered_sig[i])
				else:
					filtered_sig[i]=filtered_sig[i]+coeff[j]*saved_input[index_save]
					index_save=index_save-1
					#print ("filtered [i]=",filtered_sig[i])


			elif i-j>=0:
				filtered_sig[i]=filtered_sig[i]+coeff[j]*input_sig[i-j]

		index_save=Num_taps-2
			#else:
			#	break


	saved_input=input_sig[-(Num_taps-1):]
	#print (saved_input)


	return filtered_sig,saved_input
  
#----------------------------------------------------------------------
    
#
# use the four quadrant arctan function for phase detect between a pair of
# IQ samples; then unwrap the phase and take its derivative to demodulate
#
def fmDemodArctan(I, Q, prev_phase = 0.0):
#
# the default prev_phase phase is assumed to be zero, however
# take note in block processing it must be explicitly controlled

	# empty vector to store the demodulated samples
	fm_demod = np.empty(len(I))

	# iterate through each of the I and Q pairs
	for k in range(len(I)):

		# use the atan2 function (four quadrant version) to detect angle between
		# the imaginary part (quadrature Q) and the real part (in-phase I)
		current_phase = math.atan2(Q[k], I[k])

		# we need to unwrap the angle obtained in radians through arctan2
		# to deal with the case when the change between consecutive angles
		# is greater than Pi radians (unwrap brings it back between -Pi to Pi)
		[prev_phase, current_phase] = np.unwrap([prev_phase, current_phase])

		# take the derivative of the phase
		fm_demod[k] = current_phase - prev_phase

		# save the state of the current phase
		# to compute the next derivative
		prev_phase = current_phase

	# return both the demodulated samples as well as the last phase
	# (the last phase is needed to enable continuity for block processing)
	return fm_demod, prev_phase

# custom function for DFT that can be used by the PSD estimate
def DFT(x):

	# number of samples
	N = len(x)

	# frequency bins
	Xf = np.zeros(N, dtype='complex')

	# iterate through all frequency bins/samples
	for m in range(N):
		for k in range(N):
			Xf[m] += x[k] * cmath.exp(1j * 2 * math.pi * ((-k) * m) / N)

	# return the vector that holds the frequency bins
	return Xf

# custom function to estimate PSD based on the Bartlett method
# this is less accurate than the Welch method from matplotlib
# however, as the visual inspections confirm, the estimate gives
# the user a "reasonably good" view of the power spectrum
def estimatePSD(samples, NFFT, Fs):

	# rename the NFFT argument (notation consistent with matplotlib.psd)
	# to freq_bins (i.e., frequency bins for which we compute the spectrum)
	freq_bins = NFFT
	# frequency increment (or resolution)
	df = Fs/freq_bins

	# create the frequency vector to be used on the X axis
	# for plotting the PSD on the Y axis (only positive freq)
	freq = np.arange(0, Fs/2, df)

	# design the Hann window used to smoothen the discrete data in order
	# to reduce the spectral leakage after the Fourier transform
	hann = np.empty(freq_bins)
	for i in range(len(hann)):
		hann[i] = pow(math.sin(i*math.pi/freq_bins),2)

	# create an empty list where the PSD for each segment is computed
	psd_list = []

	# samples should be a multiple of frequency bins, so
	# the number of segments used for estimation is an integer
	# note: for this to work you must provide an argument for the
	# number of frequency bins not greater than the number of samples!
	no_segments = int(math.floor(len(samples)/float(freq_bins)))

	# iterate through all the segments
	for k in range(no_segments):

		# apply the hann window (using pointwise multiplication)
		# before computing the Fourier transform on a segment
		windowed_samples = samples[k*freq_bins:(k+1)*freq_bins] * hann

		# compute the Fourier transform using the built-in FFT from numpy
		Xf = np.fft.fft(windowed_samples, freq_bins)

		# note, you can check how MUCH slower is DFT vs FFT by replacing the
		# above function call with the one that is commented below
		#
		# Xf = DFT(windowed_samples)
		#
		# note, the slow impelementation of the Fourier transform is not as
		# critical when computing a static power spectra when troubleshooting
		# note: time permitting a custom FFT can be implemented

		# since input is real, we keep only the positive half of the spectrum
		# however, we will also add the signal energy of negative frequencies
		# to have a better a more accurate PSD estimate when plotting
		Xf = Xf[0:int(freq_bins/2)] # keep only positive freq bins
		psd_seg = 1/(Fs*freq_bins/2) * abs(Xf)**2 # compute signal power
		psd_seg = 2*psd_seg # add the energy from the negative freq bins

		# translate to the decibel (dB) scale
		for i in range(len(psd_seg)):
			psd_seg[i] = 10*math.log10(psd_seg[i])

		# append to the list where PSD for each segment is stored
		# in sequential order (first segment, followed by the second one, ...)
		psd_list.extend(psd_seg)

	# compute the estimate to be returned by the function through averaging
	psd_est = np.zeros(int(freq_bins/2))

	# iterate through all the frequency bins (positive freq only)
	# from all segments and average them (one bin at a time ...)
	for k in range(int(freq_bins/2)):
		# iterate through all the segments
		for l in range(no_segments):
			psd_est[k] += psd_list[k + l*int(freq_bins/2)]
		# compute the estimate for each bin
		psd_est[k] = psd_est[k] / no_segments

	# the frequency vector and PSD estimate
	return freq, psd_est

if __name__ == "__main__":

	# do nothing when this module is launched on its own
	#pass
	prev_I=0.0
	prev_Q=0.0
	
	I=[3,4,5,6,7,238,239,12,432,42,5,3]
	Q=[12,32,54,643,732,432,232,13,23,12,31,314]
	block_count=0
	state_phase=0
	
	fm_demod,prev_I,prev_Q=self_demod(I,Q,prev_I, prev_Q,block_count)
	#fm_demod, state_phase = fmDemodArctan(I, Q, state_phase)
	print("----------------first fm_demod-----------------")
	print("demod function result=",fm_demod)
	
	I=[13,34,25,36,47,28,39,132,222,452,23,376]
	Q=[232,332,534,6243,7424,4322,2332,133,232,132,313,313]
	block_count=1
	
	fm_demod,prev_I,prev_Q=self_demod(I,Q,prev_I, prev_Q,block_count)
	#fm_demod, state_phase = fmDemodArctan(I, Q, state_phase)
	print("----------------second fm_demod-----------------")
	print("demod function result=",fm_demod)
