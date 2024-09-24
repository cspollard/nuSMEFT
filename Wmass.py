import matplotlib.pyplot as plt
import numpy as np

f = open('Wmass.txt', 'r')
data = f.readlines()

data.pop(0)

a = np.array([])

for val in data:
    a = np.append(a, float(val.removesuffix('\n')))

r = np.max(a) - np.min(a)

counts, bins = np.histogram(a, bins = 500)
plt.stairs(counts, bins)

plt.savefig('Wmass.png')