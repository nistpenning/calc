import numpy as np
from numpy import sqrt
import matplotlib.pyplot as plt
import pylab
from scipy import stats, signal
import scipy.stats as stats


nm='batch_hist_1016_wODF_tau3000.pdf.txt'
newnm='ODF_Histogram.dat'
nbins=60
L=127

#Read in data with ODF
Isingstub='IsingHistogramN127ang'
IsingstubNoise='IsingHistogram_WithNoiseN127ang'
filehandle=open(nm,'r')
dat180=[]
dat176=[]
dat1746=[]
dat173=[]
dat90=[]
dat88=[]
dat86=[]

nsamples=0
for line in filehandle:
    nsamples+=1
    s=line.rsplit(',')
    sfl=map(float,s)
    dat180.append(sfl[0])
    dat176.append(sfl[1])
    dat1746.append(sfl[2])
    dat173.append(sfl[3])
    dat90.append(sfl[4])
    dat88.append(sfl[5])
    dat86.append(sfl[6])
print 'nsamples',nsamples
filehandle.close()

#create normalized histograms    
hist_180=np.histogram(dat180,nbins,(-1.0,1.0),normed=True)
hist_1746=np.histogram(dat1746,nbins,(-1.0,1.0),normed=True)
hist_90=np.histogram(dat90,nbins,(-1.0,1.0),normed=True)
hist_88=np.histogram(dat88,nbins,(-1.0,1.0),normed=True)


shotnoise=0.030* np.random.randn(nsamples)
histShot=np.histogram(shotnoise,nbins,(-1.0,1.0),normed=True)

shotfcn=lambda x: np.exp(-x*x/(2.0*0.030*0.030))/sqrt(2.0*np.pi*0.030*0.030)

#read in theory histograms with noise added in
nm90=IsingstubNoise+'90.000.dat'
flhandle=open(nm90,'r')
theory_90Noise=[]
for i in range(L+1):
    s=flhandle.readline()
    sfl=s.split()
    theory_90Noise.append(float(sfl[0]))
flhandle.close()

nm88=IsingstubNoise+'88.000.dat'
flhandle=open(nm88,'r')
theory_88Noise=[]
for i in range(L+1):
    s=flhandle.readline()
    sfl=s.split()
    theory_88Noise.append(float(sfl[0]))
flhandle.close()

nm1746=IsingstubNoise+'174.600.dat'
flhandle=open(nm1746,'r')
theory_1746Noise=[]
for i in range(L+1):
    s=flhandle.readline()
    sfl=s.split()
    theory_1746Noise.append(float(sfl[0]))
flhandle.close()

nm180=IsingstubNoise+'180.000.dat'
flhandle=open(nm180,'r')
theory_180Noise=[]
for i in range(L+1):
    s=flhandle.readline()
    sfl=s.split()
    theory_180Noise.append(float(sfl[0]))
flhandle.close()


npts=1000
sconv = np.linspace(-1.0,1.0,L+1)
svals = np.arange(L+1)
xvals= np.linspace(-1.0,1.0,npts)

#convolve discrete theory pdf with (continuous) photon shot noise
d1746conv=[]
for i in range(npts):
    val=0.0
    xi=xvals[i]
    for j in range(len(sconv)):
        xj=sconv[j]
        val=val+theory_1746Noise[j]*shotfcn(xi-xj)
    d1746conv.append(val)

d180conv=[]
for i in range(npts):
    val=0.0
    xi=xvals[i]
    for j in range(len(sconv)):
        xj=sconv[j]
        val=val+theory_180Noise[j]*shotfcn(xi-xj)
    d180conv.append(val)

d90conv=[]
for i in range(npts):
    val=0.0
    xi=xvals[i]
    for j in range(len(sconv)):
        xj=sconv[j]
        val=val+theory_90Noise[j]*shotfcn(xi-xj)
    d90conv.append(val)

d88conv=[]
for i in range(npts):
    val=0.0
    xi=xvals[i]
    for j in range(len(sconv)):
        xj=sconv[j]
        val=val+theory_88Noise[j]*shotfcn(xi-xj)
    d88conv.append(val)

f=open("ConvolvedTheoryPDF_N127.txt",'w')
for i in range(npts):
    f.write('%30.15E'%(xvals[i])+'%30.15E'%(d180conv[i])+'%30.15E'%(d1746conv[i])+'%30.15E'%(d90conv[i])+'%30.15E'%(d88conv[i])+'\n')    

f.close()


plt.figure(1)
plt.plot(xvals,[a for a in d1746conv])
plt.hist(dat1746,bins=nbins,range=(-1.0,1.0),normed=True)


plt.figure(2)
plt.plot(xvals,[a for a in d180conv])
plt.hist(dat180,bins=nbins,range=(-1.0,1.0),normed=True)


plt.figure(3)
plt.plot(xvals,[a for a in d90conv])
plt.hist(dat90,bins=nbins,range=(-1.0,1.0),normed=True)


plt.figure(4)
plt.plot(xvals,[a for a in d88conv])
plt.hist(dat88,bins=nbins,range=(-1.0,1.0),normed=True)

plt.show()

   
