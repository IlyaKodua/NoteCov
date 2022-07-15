from utils import *
import os
import shutil


if os.path.exists("data") and os.path.isdir("data"):
    shutil.rmtree("data")

os.mkdir("data")
sr = 16000

# mat, freq = calc_cov_for_freqs(13, 12, sr)
# print(np.median(mat))

# plot_and_save1(freq/440.0, mat)


mat1, freq = calc_corr_for_freqs_by_ebeling(1010, 1000, sr)

mat1 -= np.median(mat1)
mat1 /= np.max(mat1)
# plot_and_save2(freq/440.0, mat, )


mat2, freq = calc_cov_for_freqs(1010, 1000, sr)

mat2 -= np.median(mat2)
mat2 /= np.max(mat2)


compare(mat1,mat2,freq/440.0)

