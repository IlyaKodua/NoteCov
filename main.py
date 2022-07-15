from utils import *

sr = 16000

mat, freq = calc_corr_for_freqs(13, 12, sr)

print(np.median(mat))

plot_and_save1(freq/440.0, mat)


mat, freq = calc_corr_for_freqs(1010, 1000, sr)

plot_and_save2(freq/440.0, mat)
