import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import find_peaks


def sig_div(f, duration, sr, alph):
    t = np.arange(0,duration, 1/sr)
    sig = 1 / (1.3- np.exp(2 *1j* np.pi * f * t))* np.exp(-alph * t) 
    return np.real(sig)




def calc_corr_for_freqs(N, d, sr):
    mat = np.zeros(N+1)
    freq = np.zeros(N+1)
    for i in range(N+1):
        freq[i] = freq_by_note(i/d)
        x = sig_div(freq_by_note(0), 2, sr, 5)
        y = sig_div(freq[i], 2, sr, 5)
        mat[i] =  correlation_time(x, y)
    return mat, freq


def freq_by_note(x):
    return 440.0 * 2**(x)




def correlation_time(x,y):
    return np.abs(np.sum(x * np.conj(y))) / np.sqrt( np.sum(np.abs(x)**2) * np.sum(np.abs(y)**2))


def plot_and_save1(x, y):

    plt.plot(x, y, 'b')
    plt.plot(x, y, 'rx')
    for xs,ys in zip(x,y):

        labely = "{:.2f}".format(ys)
        label = labely.replace('.',',')

        if float(labely) > 0.58:
            plt.annotate(label, (xs,ys), # these are the coordinates to position the label
                textcoords="offset points", # how to position the text
                bbox=dict(boxstyle='round', fc='gray', alpha=0.1),
                xytext=(0,10),
                # fontsize = 7,
                ha='center') # horizontal alignment can be left, right or center

    plt.ylim([0.57, 1.05])
    plt.xlabel("Frequency ratio")
    plt.ylabel("Normalized covariance")
    plt.savefig("figs/fig1.pdf", dpi=300, bbox_inches='tight')

    plt.show()
    
    pass


def plot_and_save2(x, y):

    plt.plot(x, y, 'b')

    peaks, _ = find_peaks(y, distance=20)
    plt.plot(x[peaks], y[peaks], 'rx')
    for xs,ys in zip(x[peaks],y[peaks]):

        labely = "{:.3f}".format(ys)
        label = labely.replace('.',',')

        if float(labely) >= 0.59:
            plt.annotate(label, (xs,ys), # these are the coordinates to position the label
                textcoords="offset points", # how to position the text
                bbox=dict(boxstyle='round', fc='gray', alpha=0.1),
                xytext=(0,10),
                # fontsize = 7,
                ha='center') # horizontal alignment can be left, right or cente

    plt.ylim([0.57, 1.05])
    plt.xlabel("Frequency ratio")
    plt.ylabel("Normalized covariance")
    plt.savefig("figs/fig2.pdf", dpi=300, bbox_inches='tight')

    plt.show()
    
    pass