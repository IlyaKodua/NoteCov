import numpy as np
from matplotlib import pyplot as plt
import librosa
import math



def get_spectrum(y, sr):
    dt = 1/sr
    duration = dt*len(y)  
    N = len(y)
    dF = sr/len(y)

    Y = np.fft.fftshift(np.fft.fft(y))          
    f = np.arange(-sr/2, sr/2,dF)
    return Y, f

def sound_gen(f, duration, sr, alph):
    t = np.arange(0,duration, 1/sr)
    sig = np.exp(2 *1j* np.pi * f * t) * np.exp(-alph * t)
    return sig

def sig_ln(f, duration, sr, alph):
    t = np.arange(0,duration, 1/sr)
    sig = np.real(np.log(1.3 - np.exp(2 *1j* np.pi * f * t)))* np.exp(-alph * t)   
    return sig


def sig_exp(f, duration, sr, alph):
    t = np.arange(0,duration, 1/sr)
    sig = np.real(np.exp(np.exp(2 *1j* np.pi * f * t)))* np.exp(-alph * t)
    return sig

def sig_div(f, duration, sr, alph):
    t = np.arange(0,duration, 1/sr)
    sig = 1 / (1.3- np.exp(2 *1j* np.pi * f * t))* np.exp(-alph * t) 
    return np.real(sig)

def sig_syntetic(f, duration, sr, alph):
    t = np.arange(0,duration, 1/sr)
    N = 100
    sig =  np.exp(2 * 1j *np.pi * f * t * N) - np.exp(-2 * 1j *np.pi * f * t * N) / ( np.exp(2 * 1j *np.pi * f * t) - 1)  * np.sinc(10*t)
    sig[0] = 10
    return sig


def correlation_freq(x,y, sr):
    x,_ = get_spectrum(x, sr)
    y,_ = get_spectrum(y, sr)


    x = np.fft.ifftshift(x)[0:int(len(x)/2)]
    y = np.fft.ifftshift(y)[0:int(len(y)/2)]
    # plt.plot(x)
    # plt.show()



    return np.abs(np.sum(x * np.conj(y))) / np.sqrt( np.sum(np.abs(x)**2) * np.sum(np.abs(y)**2))


def calc_corr_for_freqs(N, d, sr):
    mat = np.zeros(N+1)
    freq = np.zeros(N+1)
    for i in range(N+1):
        freq[i] = freq_by_note(i/d)
        x = sig_exp(freq_by_note(0), 2, sr, 5)
        y = sig_exp(freq[i], 2, sr, 5)
        mat[i] =  correlation_time(x, y, sr)
    return mat, freq


def calc_covq(N, sr):
    q_array = []
    cov_array = []
    for i in range(1,N):
        for j in range(1,N):
            if i > j and i/j <= 2 and math.gcd(i,j) == 1:
                x = sig_exp(440.0, 2, sr, 5)
                y = sig_exp(440.0 * i / j, 2, sr, 5)
                cov_array.append(correlation_time(x, y, sr))
                q_array.append(str(i) + " : " + str(j))
    cov_array = np.array(cov_array)
    q_array = np.array(q_array)
    inds = cov_array.argsort()
    q_array = q_array[inds]
    cov_array = cov_array[inds]
    return cov_array, q_array

def freq_by_note(x):
    return 440.0 * 2**(x)


def conv_time(x,y,sr):
    return np.max(np.convolve(x,y)) / np.sqrt( np.sum(np.abs(x)**2) * np.sum(np.abs(y)**2))

def correlation_time(x,y,sr):
    return np.real(np.sum(x * np.conj(y))) / np.sqrt( np.sum(np.abs(x)**2) * np.sum(np.abs(y)**2))



def calculate_freq_diff(f1, f2, N):

    x = f1 * np.arange(N).reshape((N,1))
    y = f2 * np.arange(N).reshape((1,N))
    mat_x = np.repeat(x,N, axis=1)
    mat_y = np.repeat(y,N, axis=0)
    return np.mean(1/((mat_x - mat_y)**2 + 1))