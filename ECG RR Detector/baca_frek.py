import numpy as np
import matplotlib.pyplot as plt

#baca data
#parent_folder = "D:\\GoogleDriveAlvin\\MyResearchGroup\\SupineStandingECG\\Supine\\"
parent_folder1 = "//media//nomad2805//LENOVO/GoogleDriveAlvin//MyResearchGroup//SupineStandingECG/Standing//" 
parent_folder2 = "//media//nomad2805//LENOVO/GoogleDriveAlvin//MyResearchGroup//SupineStandingECG/Supine//" 
file_name = "RRintS5r"
dataRR1 = np.genfromtxt(parent_folder1 + file_name + ".txt", delimiter = " ")
dataRR2 = np.genfromtxt(parent_folder2 + file_name + ".txt", delimiter = " ")
y1 = np.fft.fft(dataRR1)
y2 = np.fft.fft(dataRR2)
yabs1 = np.abs(y1)
yabs2 = np.abs(y2)
freq1 = np.linspace(0,200,len(yabs1))
freq2 = np.linspace(0,200,len(yabs2))
print(len(freq1),len(freq2),len(yabs1),len(yabs2))
plt.figure()
plt.plot(freq1[1:round(len(yabs1)/2)],yabs1[1:round(len(yabs1)/2)])
plt.plot(freq2[1:round(len(yabs2)/2)],yabs2[1:round(len(yabs2)/2)])
plt.show()
