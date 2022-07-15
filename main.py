from utils import *
from matplotlib import pyplot as plt
from scipy.io.wavfile import write
sr = 16000

mat, freq = calc_corr_for_freqs(201, 200, sr)


plt.plot(freq/ 440.0, mat, 'b')
plt.plot(freq/ 440.0, mat, 'rx')
print(np.median(mat))

# plt.savefig("graph.png")
# mask = mat > np.median(mat)
# print(mat[mask])

# plt.plot(freq[mask]/ 440.0, mat[mask], 'b+')

plt.xlabel("Frequency ratio")
plt.ylabel("Normalized covariance")
plt.savefig("graph.png")
plt.show()