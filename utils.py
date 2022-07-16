from operator import index
import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import find_peaks
from scipy import signal
from scipy import stats


def sig_div(f, duration, sr, alph):
    t = np.arange(0,duration, 1/sr)
    sig = 1 / (1.3- np.exp(2 *1j* np.pi * f * t))* np.exp(-alph * t) 
    return np.real(sig)


def sig_simple(f, duration, sr, alph):
    t = np.arange(0,duration, 1/sr)
    sig = np.sin(2 * np.pi * f * t) * np.exp(-alph * t) 
    return sig




def calc_corr_for_freqs_by_ebeling(N, d, sr):
    mat = np.zeros(N)
    freq = np.zeros(N)
    for i in range(N):
        freq[i] = freq_by_note(i/d) / freq_by_note(0)
        x = sig_div(freq_by_note(0), 2, sr, 5)
        y = sig_div(freq_by_note(i/d), 2, sr, 5)
        mat[i] =  generalized_coincidence_function(x, y)
    mat /= mat[0]
    mat -= np.median(mat)
    mat /= np.max(mat)
    return mat, freq



def calc_cov_for_freqs(N, d, sr):
    mat = np.zeros(N)
    freq = np.zeros(N)
    for i in range(N):
        freq[i] =  freq_by_note(0)
        x = sig_div(freq_by_note(0), 2, sr, 5)
        y = sig_div(freq_by_note(i/d), 2, sr, 5)
        mat[i] =  correlation_time(x, y)
    mat[i] /= mat[0]
    mat -= np.median(mat)
    mat /= np.max(mat)
    return mat, freq


def freq_by_note(x):
    return 440.0 * 2**(x)



def calc_dist_for_freqs(N, d, sr):
    mat = np.zeros(N)
    freq = np.zeros(N)
    for i in range(N):
        freq[i] = freq_by_note(i/d)
        freqs = np.concatenate([freq_by_note(i/d) * np.arange(0,10), freq_by_note(0) * np.arange(0,10)])
        ampls =  np.concatenate([np.exp( -np.log(1.3) * np.arange(0,10)), np.exp( -np.log(1.3) * np.arange(0,10))]) 
        
        mat[i] =  dissmeasure(freqs, ampls, "product")
    mat *= -1
    mat -= np.min(mat)
    mat /= np.max(mat)
    return mat, freq


def correlation_time(x,y):
    return np.abs(np.sum(x * np.conj(y))) / np.sqrt( np.sum(np.abs(x)**2) * np.sum(np.abs(y)**2))

def  generalized_coincidence_function(sig1, sig2):
    s = sig1 + sig2
    corr = signal.correlate(s, s)
    return np.sum(corr**2)


def teoretical(N, d):

    
    vals = np.zeros(N, dtype=np.complex64)
    freqs = np.zeros(N, np.float32)
    for n in range(N):
        for k in range(20):
            for m in range(20):
                vals[n] += 1.3**(-k - m) / (10 + 1j*(m*freq_by_note(n/d) - k * freq_by_note(0)))
        freqs[n] = freq_by_note(n/d) /  freq_by_note(0)
    get_vals = np.real(vals)
    get_vals -= np.median(get_vals)
    get_vals /= np.max(get_vals)
    return get_vals, freqs


def dissmeasure(fvec, amp, model='min'):
    #taken from https://gist.github.com/endolith/3066664
    """
    Given a list of partials in fvec, with amplitudes in amp, this routine
    calculates the dissonance by summing the roughness of every sine pair
    based on a model of Plomp-Levelt's roughness curve.
    The older model (model='product') was based on the product of the two
    amplitudes, but the newer model (model='min') is based on the minimum
    of the two amplitudes, since this matches the beat frequency amplitude.
    """
    # Sort by frequency
    sort_idx = np.argsort(fvec)
    am_sorted = np.asarray(amp)[sort_idx]
    fr_sorted = np.asarray(fvec)[sort_idx]

    # Used to stretch dissonance curve for different freqs:
    Dstar = 0.24  # Point of maximum dissonance
    S1 = 0.0207
    S2 = 18.96

    C1 = 5
    C2 = -5

    # Plomp-Levelt roughness curve:
    A1 = -3.51
    A2 = -5.75

    # Generate all combinations of frequency components
    idx = np.transpose(np.triu_indices(len(fr_sorted), 1))
    fr_pairs = fr_sorted[idx]
    am_pairs = am_sorted[idx]

    Fmin = fr_pairs[:, 0]
    S = Dstar / (S1 * Fmin + S2)
    Fdif = fr_pairs[:, 1] - fr_pairs[:, 0]

    if model == 'min':
        a = np.amin(am_pairs, axis=1)
    elif model == 'product':
        a = np.prod(am_pairs, axis=1)  # Older model
    else:
        raise ValueError('model should be "min" or "product"')
    SFdif = S * Fdif
    D = np.sum(a * (C1 * np.exp(A1 * SFdif) + C2 * np.exp(A2 * SFdif)))

    return D


# def get_spectrum(x):
    
#     spectrum = np.abs(np.fft.fft(x))
#     spectrum = spectrum[0:int(len(spectrum)/2)]

#     spectrum[1:len(spectrum)] *= 2
#     spectrum /= len(x)
#     return spectrum

# def get_ample_and_phase(x):
#     y = np.fft.fft(x)
#     y[0:int(len(y)/2)] = 0
#     return np.fft.ifft(y)

# def metric_of_harm(x1, x2):
#     x = x1 + x2
#     y = get_spectrum(x)
#     y = np.diff(y)
#     z = get_spectrum(y)

#     return np.max(z) / np.sum(np.abs(z)**2)**(1/2)

def new_func(y):
    plt.figure(1)
    plt.plot(y, 'b')

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
                fontsize = 7,
                ha='center') # horizontal alignment can be left, right or center

    plt.ylim([0.57, 1.05])
    plt.xlabel("Frequency ratio")
    plt.ylabel("Normalized covariance")
    plt.savefig("data/fig1.pdf", dpi=300, bbox_inches='tight')

    plt.show()
    
    pass



def plot_and_save2(x, y, name):

    plt.plot(x, y, 'b')


    # plt.ylim([0.57, 1.05])
    plt.xlabel("Frequency ratio")
    plt.ylabel(name)
    plt.savefig("data/fig2.pdf", dpi=300, bbox_inches='tight')

    plt.show()
    
    pass


def get_less_median(x, freq):
    x_new  = []
    freq_new = []
    for i in range(len(x)):
        if np.abs(x[i] - np.median(x)) > 1e-3 and x[i] < np.median(x):
            x_new.append(x[i])
            freq_new.append(freq[i])
    x_new = np.array(x_new)
    freq_new = np.array(freq_new)
    x_new, freq_new = sorting(x_new, freq_new)
    return x_new, freq_new


def get_more_median(x, freq):
    x_new  = []
    freq_new = []
    for i in range(len(x)):
        if np.abs(x[i] - np.median(x)) > 1e-2 and x[i] > np.median(x):
            x_new.append(x[i])
            freq_new.append(freq[i])
    x_new = np.array(x_new)
    freq_new = np.array(freq_new)
    # x_new, freq_new = sorting(x_new, freq_new)
    return x_new, freq_new

def sorting(x, freqs):

    freqs_uniq, indexes = np.unique(freqs.round(decimals=2), return_index=True)
    return_freqs = []
    return_x = []
    for i in range(len(freqs_uniq)):
        mask = freqs_uniq[i] == freqs.round(decimals=2)
        if np.sum(mask) >= 1:
            freqs_mask = freqs[mask]
            indx = np.argmax(x[mask])
            return_freqs.append(freqs_mask[indx])

    for i in range(len(return_freqs)):
        return_x.append(x[return_freqs[i] == freqs][0])
    
    return  return_x, return_freqs


def write_txt(x, freq):

    x_new, freq_new = get_less_median(x, freq)
    with open('data/less.txt', 'w') as f:
        for i in range(len(x_new)):
                txt_x = "{:.5f}".format(x_new[i])
                txt_freq = "{:.5f}".format(freq_new[i])
                f.write(txt_x + " " + txt_freq)
                f.write('\n')

    x_new, freq_new = get_more_median(x, freq)
    with open('data/more.txt', 'w') as f:
        for i in range(len(x_new)):
                txt_x = "{:.5f}".format(x_new[i])
                txt_freq = "{:.5f}".format(freq_new[i])
                f.write(txt_x + " " + txt_freq)
                f.write('\n')
        

def teoretical(N, d):

    
    vals = np.zeros(N, dtype=np.complex64)
    freqs = np.zeros(N, np.float32)

    for n in range(N):
        freqs[n] = freq_by_note(n/d) /     freq_by_note(0)
        for k in range(20):
            for m in range(20):
                vals[n] += 1.3**(-k - m) / (10 + 1j*(m*freq_by_note(n/d) - k * freq_by_note(0)))
    get_vals = np.real(vals)
    get_vals -= np.median(get_vals)
    get_vals /= np.max(get_vals)
    return get_vals, freqs


def compare(x,y,z, freq):
    plt.plot(freq, x, 'b',   label="Normalized covariance by formula", linewidth=0.5)
    plt.plot(freq, y, 'r',   label="Normalized covariance", linewidth=0.5)
    plt.plot(freq, z, 'black', label="GCF", linewidth=0.5)
    # plt.plot(freq, z, 'black', label="Teoretical curve", linewidth=0.5)
    plt.legend(loc="upper right")
    plt.ylim([-0.1, 1.1])
    plt.xlabel("Frequency ratio")
    # plt.ylabel(name)
    plt.savefig("data/fig3.pdf", dpi=300, bbox_inches='tight')


    plt.show()



# def dicconanse_plot(x,freq):

#     plt.plot(freq, x, 'b', linewidth=0.5)
#     _, freq_new = get_less_median(x, freq)
#     for f in freq_new:
#         plt.axvline(x = f, color = 'r', linewidth=0.5)
#     # _, freq_new = get_more_median(x, freq)
#     # for f in freq_new:
#     #     plt.axvline(x = f, color = 'b', linewidth=0.5)
#     # plt.legend(loc="upper right")
#     plt.ylim([-0.1, 1.1])
#     plt.xlabel("Frequency ratio")
#     plt.ylabel("Normalized covariance")
#     plt.savefig("data/fig2.pdf", dpi=300, bbox_inches='tight')


#     plt.show()


