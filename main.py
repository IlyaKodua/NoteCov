from utils import *
import os
import shutil


if os.path.exists("data") and os.path.isdir("data"):
    shutil.rmtree("data")

os.mkdir("data")
sr = 16000


# get_spectrum(sig_div(freq_by_note(0), 2, sr, 5))
mat1, freq = calc_cov_for_freqs(13, 12, sr)
# print(np.median(mat))

plot_and_save1(freq, mat1)



# mat1, freq = calc_corr_for_freqs_by_ebeling(1010, 1000, sr)
mat1, freq = teoretical(1010, 1000)

# plot_and_save2(freq/440.0, mat, )
mat2, freq = calc_cov_for_freqs(1010, 1000, sr)

mat3, freq = calc_corr_for_freqs_by_ebeling(1010, 1000, sr)


# write_txt(mat2, freq)

# dicconanse_plot(mat2, freq)

# mat3, freq = calc_dist_for_freqs(1010, 1000, sr)


# print(stats.pearsonr(mat2, mat3))
# print(stats.pearsonr(mat1, mat3))
compare(mat1,mat2,mat3, freq)

