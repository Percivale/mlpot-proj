import numpy as np
import matplotlib.pyplot as plt

#path = input("path to lcurve file: ")
path = ""
data = np.genfromtxt(path + "lcurve.out", names=True)

print(data.dtype.names)
print(data.shape)
plt.figure()
plt.title("Loss as a funtion of training length", size=18)
for name in data.dtype.names[3:-1]:
    plt.plot(data['step'][::10], data[name][::10], label=name)
plt.legend()
plt.xlabel('Step', size=16)
plt.ylabel('Loss', size=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xscale('symlog')
plt.xlim(0, 12**5)
plt.yscale('log')
plt.grid()
plt.show()
#plt.savefig('lcurve.png', transparent=False)

plt.figure()
plt.title("Loss as a funtion of training length", size=18)
for name in data.dtype.names[3:-1]:
    plt.plot(data['step'][::5], data[name][::5], label=name)
plt.legend()
plt.xlabel('Step', size=16)
plt.ylabel('Loss', size=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xscale('symlog')
plt.xlim(10**3, 10**5)
plt.yscale('log')
plt.grid()
plt.show()
