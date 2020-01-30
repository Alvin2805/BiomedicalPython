import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal

def filtered(non_filtered_data):
    n = 347
    fc = 15 # frequency cut off
    fx = 177 # frequency sampling
    a = signal.firwin(n, cutoff=[float(fc) / fx], window="hamming")
    filtered_signal = signal.filtfilt(a, 1, non_filtered_data)
    return filtered_signal


def detectRR(dataECG):
    f = open("hasilRR.txt", "w")
    #meanpeak = np.mean(dataECG)
    meanpeak = 525
    #print(meanpeak)
    i = 0
    rridx = []
    rr = []
    rtemp = 0
    while i < len(dataECG):
        peak_test = dataECG[i]
        if peak_test > meanpeak:
            window = 120
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
    nstart = 100
    nstop = 150000
    SubjectID = "S5r"
    folder_parent = "D:\\GoogleDriveAlvin\\MyResearchGroup\\SupineStandingECG\\Supine\\"
    #folder_parent = "//media//nomad2805//LENOVO//GoogleDriveAlvin//MyResearchGroup//SupineStandingECG//Standing//"
    raw = np.genfromtxt(folder_parent + SubjectID + ".txt", delimiter=" ")
    data_filter = filtered(raw)
    data_detect = data_filter[nstart:nstop]

    #print(data_detect)

    #write the data
    rr, rridx = detectRR(data_detect)
    tulis = open(folder_parent + "RRint" + SubjectID + ".txt","w")
    for i in range(0,len(rridx)-1):
        tulis.write(str(rridx[i+1] - rridx[i]) + "\n")
    tulis.close()
    plt.figure()
    plt.plot(data_detect)
    plt.plot(rridx, rr, "ro")
    plt.show()

if __name__ == "__main__":
    main()
