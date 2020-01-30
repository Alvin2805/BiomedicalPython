import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
import sys
import codecs
from pylab import *
from spectrum import *
#import corren #only used for Correntropy Processing (experimental)

#===================================================================

def findmax(xaxix,vmax):
	n = len(xaxix)
	for search in range(0,n):
		if (xaxix[search] >= vmax):
			found = search
			return found
			break


def AR_Spectrum(dataAR,norder,fsamp):
	AR,P,k = aryule(dataAR,norder)
	PSD = arma2psd(AR,NFFT=len(dataAR))
	PSD = PSD[len(PSD):len(PSD)//2:-1]
	f = fsamp*(linspace(0,1,len(PSD)))
	PSDout = abs(PSD)*2./(2.*pi)
	return f,PSDout

def singleread(dataread,detype):
	if (detype==1):
		data = np.genfromtxt(dataread,delimiter="\t")
	elif (detype==2):
		data = np.genfromtxt(dataread,delimiter=" ")
	else:
		data = np.genfromtxt(dataread,delimiter=",")
	n = len(data)
	t = np.linspace(0,(1/fsamp)*len(data),len(data))
	print("Data Found, Please wait...")
	print("Length Data = %d Samples/data" %n)

	return data,t

def tgamread(dataread,eegdata):
	try:
		data = codecs.open(dataread,"r","utf-8",errors="ignore")
		print("Data Found, Please wait...")
		ch1=[]
		ch2=[]
		idx = 0
		delta = []
		theta = []
		alpha = []
		beta = []
		gamma = []

		for line in data:
			read = line.split(",")
			#print(read)
			TP = 0
			if (read[0]=="\r\n") or (len(read)==1) or (read[0]==" "):
				temp = read
			elif (len(read)>2):
				for j in range(3,11):
				    TP += int(read[j-1])
				TP_add = read[10].splitlines()
				TP += int(TP_add[0])
				delta.append(int(read[3])/(TP*1.0))
				theta.append(int(read[4])/(TP*1.0))
				alpha.append(int(read[5])/(TP*1.0)+int(read[6])/(TP*1.0))
				beta.append(int(read[7])/(TP*1.0)+int(read[8])/(TP*1.0))
				gamma.append(int(read[9])/(TP*1.0)+int(TP_add[0])/(TP*1.0))
			    #print(TP)
			else:
				baca = read[1].splitlines()
				#print(baca)
				ch1.append(int(read[0]))
				ch2.append(int(baca[0]))
				#print(int(read[0]),int(baca[0]))
			#print(idx)
			idx += 1
		print("Data Length : " + str(idx) + " Samples/data")
		data.close()

		# ch1 Rectified (EMG) ch2 (Raw/ECG/EGG)
		ch1 = np.asarray(ch1)
		ch2 = np.asarray(ch2)
		t1 = np.linspace(0,(1/500)*len(ch1),len(ch1))
		t2 = np.linspace(0,(1/500)*len(ch2),len(ch2))

		ttt = signal.firwin(200,20,width=None,window="hamming",pass_zero="lowpass",fs=500)
		chraw = signal.filtfilt(ttt,1,ch2)

		# EEG data
		delta = np.asarray(delta)
		theta = np.asarray(theta)
		alpha = np.asarray(alpha)
		beta = np.asarray(beta)
		gamma = np.asarray(gamma)

		if (eegdata == 1):
			printf("++ EEG data extracted")
			return ch1,ch2,delta,theta,alpha,beta,gamma
		else:
			return ch1,chraw,t1
	except IOError:
		print("SORRY, File not Existed, Please Re-Check")

def preprocessing(datainp,fsamp,fdown):
	f1 = 0.015
	f2 = 0.15
	numtap = 50

	#detrend
	datainp = signal.detrend(datainp)

	#h = signal.firwin(numtap,[f1,f2],width=None, window="hamming",pass_zero="bandpass",scale=False,fs=fsamp)
	b,a = signal.butter(6,f2,btype="lowpass",analog=False,output="ba",fs=fsamp)
	w,h = signal.freqz(b,a,fs=fsamp)
	plt.figure(5)
	plt.plot(w,abs(h))
	#ch_clean = signal.filtfilt(h,1,datainp)
	ch_clean = signal.filtfilt(b,a,datainp)
	t_clean = np.linspace(0,(1/fsamp)*len(ch_clean),len(ch_clean))

	#ch_decimate = signal.decimate(ch_clean,round(fsamp/fdown),n=None,ftype="fir",zero_phase=False)
	ch_decimate = downsampling(ch_clean,fsamp,fdown)
	
	return ch_decimate,ch_clean
	
def preprocessing2(datainp,fsamp,fdown):
	f1 = 0.015
	f2 = 0.15

	#detrend data (must and compulsory)
	datainpp = signal.detrend(datainp)
	
	ch_decimate = downsampling(datainpp,fsamp,fdown)

	b,a = signal.butter(6,[f1,f2],btype="bandpass",analog=False,output="ba",fs=fdown)
	w,h = signal.freqz(b,a,fs=fdown)
	plt.figure(5)
	plt.plot(w,abs(h))
	ch_clean = signal.filtfilt(b,a,ch_decimate)
	t_clean = np.linspace(0,(1/fsamp)*len(ch_clean),len(ch_clean))

	#ch_decimate = signal.decimate(ch_clean,round(fsamp/fdown),n=None,ftype="fir",zero_phase=False)  # based on package tool, using downsampling for scratch
	
	return ch_clean,ch_decimate

def downsampling(origin,forigin,flast):
	factor = round(forigin/flast)
	ndata = len(origin)
	cek = 0
	downdata = []
	while (cek < ndata):
		downdata.append(origin[cek])
		cek += factor
	return downdata

def freqanalysis(chaninp,fsamp):
	yfft = np.fft.fft(chaninp,n=len(chaninp))
	yabs = abs(yfft)
	f = np.linspace(0,fsamp,round(len(yabs)/2))
	half_F = round(len(chaninp)/2)

	return f[0:half_F],yabs[0:half_F]

if __name__ == "__main__":
	fsamp = 200
	fdown = 2
	tsegmen = 20 # minute
	
	#for extracting the EEG data from INSINAS device
	#ch1,ch2,t = tgamread(sys.argv[1],eegdata=0)
	
	chread,t = singleread(sys.argv[1],0)
	
	ch2 = chread[:,1]
	
	#nsegmen = fsamp*tsegmen*60
	nsegmen = len(ch2)
	
	nepoch = round(len(ch2)/nsegmen)
	print("Data segment found : %d" %(nepoch))
	colorchange = 0
	
	for seg in range(0,nepoch):
		start_seg = seg*nsegmen
		stop_seg = (seg+1)*nsegmen

		#print(len(ch2[start_seg:stop_seg]))
		#ch2 = corren.mean_correntropy(ch2[start_seg:stop_seg])

		ch_fix,ch = preprocessing2(ch2[start_seg:stop_seg],fsamp,fdown)
		f,y = freqanalysis(ch_fix,fdown)

		plt.figure(seg)
		plt.subplot(411)
		plt.plot(t,ch2)
		plt.xlim([0,30])
		#plt.ylim([480,530])
		plt.subplot(412)
		plt.plot(ch)
		#plt.ylim([-0.2,0.2])
		plt.subplot(413)
		tfix = np.linspace(0,(1/fdown)*len(ch_fix),len(ch_fix))
		plt.plot(ch_fix)
		#plt.ylim([-0.2,0.2])
		plt.subplot(414)
		plt.plot(f,y)
		#plt.xlim([0,0.2])
		#plt.ylim([0,2])

		far,psdar = AR_Spectrum(ch_fix,round(0.2*(len(ch_fix)-1)),fdown)

		TPpsdar = np.sum(psdar)
		psdar = psdar/TPpsdar

		TPy = np.sum (y)
		y = y/TPy
		

		#for plot set
		maxy1 = max(y)
		maxy2 = max(psdar)
		
		
		maxfind = findmax(far,0.1)
		fmax = np.argmax(psdar[:maxfind])
		print("Freq Max : %.4f with value is : %.4f" %(far[fmax],psdar[fmax]))

		plt.figure(1)
		
		plt.subplot(211)
		
		
		plt.plot(f,y)
		plt.fill_between(f,y,color="k")
		plt.ylim([0,0.1])
		plt.xlim([0,0.3])
		plt.subplot(212)
		plt.plot(far,psdar)		
		plt.fill_between(far,psdar,color="k")
		plt.ylim([0,0.1])
		plt.xlim([0,0.3])
		
		colorchange += 100
		
		#calculate the stomach PSD at interval frequency
		low_freq = findmax(far,0.03)
		high_freq = findmax(far,0.07)
		PSD = np.sum(psdar[low_freq:high_freq])
		print("PSD on stomach = %.4f" %PSD)
		
plt.show()
