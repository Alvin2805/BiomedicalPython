import numpy as np
import matplotlib.pyplot as plt


def freqfind(data, value):
	n = len(data)
	for i in range(0,n):
		if (data[i] == value):
			indexdata = i
			return indexdata
			break
		elif (data[i] > value):
			indexdata = i
			return indexdata
			break

def LinInterpolation(data,interval):
	size_series = len(data)
	peak_point = np.zeros([size_series])
	peak_point[0] = data[0]
	#print(peak_point.shape)
	
	for t in range(1,size_series):
		peak_point[t] = peak_point[t-1] + data[t-1]
		
	itprrlist = []
	#print(int(peak_point[0]))
	#print(int(max(peak_point)))
	for i in range(int(peak_point[0]), int(max(peak_point)), interval):
		itprrlist.append(np.interp(i,peak_point, data))
	
	return itprrlist

def FreqOri(rrdata,fs):
	interpRange = 1000 # in ms usually fs is 1 seconds (1000ms) for obtaining PSD range 0 - 0.5 Hz
	intpData = LinInterpolation(rrdata, interpRange)
	datainp = np.array(intpData)
	PS = np.fft.fft(datainp/fs,n=128)
	yabs = np.abs(PS)
	freq = np.linspace(0,1,len(yabs))
	freq_vlf = freqfind(freq,0.04) # < 0.04 Hz
	freq_lf = freqfind(freq,0.15) # 0.04 - 0.15 Hz
	freq_hf = freqfind(freq,0.4) # 0.15 - 0.4 Hz
	
	VLF = np.sum(yabs[1:freq_vlf])
	LF = np.sum(yabs[freq_vlf:freq_lf])
	HF = np.sum(yabs[freq_lf:freq_hf])
	TP = np.sum(yabs[1:freq_hf])
	LFHF = LF/HF
	#print(VLF,LF,HF)
	plt.figure()
	plt.plot(freq[0:round(len(yabs)/2)],yabs[0:round(len(yabs)/2)])
	plt.ylim([0,5])
	plt.show()
	
	return VLF,LF,HF,TP,LFHF
	
