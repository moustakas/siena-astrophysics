import matplotlib.pylab as plt

plt.hist(data['Z'],bins=1000,range=(0.0,1.0))
plt.hist(data[tot1]['Z'],bins=1000,range=(0.0,1.0))
plt.hist(data[tot2]['Z'],bins=1000,range=(0.0,1.0))
plt.hist(data[tot3]['Z'],bins=1000,range=(0.0,1.0))
plt.hist(data[tot4]['Z'],bins=1000,range=(0.0,1.0))
plt.hist(data[tot5]['Z'],bins=1000,range=(0.0,1.0))
plt.hist(data[tot6]['Z'],bins=1000,range=(0.0,1.0))
plt.hist(data[tot7]['Z'],bins=1000,range=(0.0,1.0))
plt.ylim(0,7500)
plt.show()
