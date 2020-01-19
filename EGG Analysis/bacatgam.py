import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
import sys
import codecs
from pylab import *
from spectrum import *

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
    print("Data Found, Please wait...")
    print("Length Data = %d Samples/data",%n)

    return data

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
    f1 = 0.0015
    f2 = 0.15
    numtap = 2000

    #detrend
    datainp = signal.detrend(datainp)

    #h = signal.firwin(numtap,[f1,f2],width=None, window="hamming",pass_zero="bandpass",scale=False,fs=fsamp)
    h = signal.butter(7,0.15/fsamp,btype="low",output="sos")
    #f,w = signal.freqz(b,a,fs=fsamp)
    #ch_clean = signal.filtfilt(h,1,datainp)
    ch_clean = signal.sosfilt(h,datainp)
    t_clean = np.linspace(0,(1/500)*len(ch_clean),len(ch_clean))

    ch_decimate = signal.decimate(ch_clean,round(fsamp/fdown),n=None,ftype="fir",zero_phase=False)
    #ch_decimat = signal.resample(ch_clean,round((fdown*len(datainp))/fsamp)) 

    return ch_decimate,ch_clean

def freqanalysis(chaninp,fsamp):
    yfft = np.fft.fft(chaninp,n=len(chaninp))
    yabs = abs(yfft)
    f = np.linspace(0,fsamp,round(len(yabs)/2))
    half_F = round(len(chaninp)/2)

    return f[0:half_F],yabs[0:half_F]

if __name__ == "__main__":
    fsamp = 500
    fdown = 4
    ch1,ch2,t = tgamread(sys.argv[1],eegdata=0)
    ch_fix,ch = preprocessing(ch2,fsamp,fdown)
    f,y = freqanalysis(ch_fix,fdown)
    
    plt.figure()
    plt.subplot(411)
    plt.plot(t,ch2)
    plt.xlim([0,30])
    plt.ylim([480,530])
    plt.subplot(412)
    plt.plot(ch)
    plt.ylim([-0.2,0.2])
    plt.subplot(413)
    tfix = np.linspace(0,(1/fdown)*len(ch_fix),len(ch_fix))
    plt.plot(ch_fix)
    plt.ylim([-0.2,0.2])
    plt.subplot(414)
    plt.plot(f,y)
    plt.ylim([0,25])
    
    far,psdar = AR_Spectrum(ch_fix,len(ch_fix)-1,fdown)

    plt.figure()
    plt.subplot(211)
    plt.plot(f,y)
    plt.ylim([0,70])
    plt.xlim([0,0.15])
    plt.subplot(212)
    plt.plot(far,psdar)
    #plt.ylim([0,70])
    #plt.xlim([0,0.15])
    
    
    plt.show()
