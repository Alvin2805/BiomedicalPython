import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal

def filtered(non_filtered_data):
	n = 300
	fc = 15 # frequency cut off
	fx = 250 # frequency sampling
	b = signal.firwin(n, fc, width=None, window="hamming", pass_zero="lowpass", fs=fx)
	a = 1
	w,h = signal.freqz(b,a,fs = fx)
	filtered_signal = signal.filtfilt(b, a, non_filtered_data)
	return filtered_signal


def detectRR(dataECG):
	f = open("hasilRR.txt", "w")
	#meanpeak = np.mean(dataECG)
	meanpeak = 550
	#print(meanpeak)
	i = 0
	rridx = []
	rr = []
	rtemp = 0
	while i < len(dataECG):
		peak_test = dataECG[i]
		if peak_test > meanpeak:
			window = 100
			datacheck = dataECG[i:i + window]
			rpeak = max(datacheck)
			rindex = i + np.argmax(datacheck)
			interv = rindex - rtemp
			if interv > 80:
				rr.append(rpeak)
				rridx.append(rindex)
				i = i + window + 1
				rtemp = rindex
				f.write(str(rindex) + "\n")
			else:
				i = i + 1
		else:
			i = i + 1
	f.close()
	return rr, rridx

def main():
	raw = np.genfromtxt("S2_post2.txt", delimiter="\t")
	print(raw)
	data_filter = filtered(-raw[:,19])
	
	data_detect = signal.detrend(data_filter)
	
	print("Data Check ---> ")
	
	print(data_detect)

	#write the data
	rr, rridx = detectRR(data_detect)
	tulis = open("RRint.txt","w")
	for i in range(0,len(rridx)-1):
		tulis.write(str(rridx[i+1] - rridx[i]) + "\n")
	tulis.close()
	plt.figure()
	plt.plot(data_detect)
	plt.plot(rridx, rr, "ro")
	plt.show()

if __name__ == "__main__":
	main()
