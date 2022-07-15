from utils import *

sr = 16000

mat, freq = calc_corr_for_freqs(13, 12, sr)

print(np.median(mat))

plot_and_save("fig1.pdf", freq/440.0, mat)