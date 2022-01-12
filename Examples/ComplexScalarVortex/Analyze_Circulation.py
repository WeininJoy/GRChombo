import numpy as np
import matplotlib.pyplot as plt

# read flash.dat to a list of lists
datContent = [i.strip().split() for i in open("./outputs_1000.dat").readlines()]
m_omega = 1.0

cir_list = []
large_arg = []

for i in range(len(datContent)):
    cir = 0.0
    for j in range(len(datContent[i])):
        datContent[i][j] = float(datContent[i][j]) * ( 2.0*3.1416/len(datContent[i]) )
        if abs(datContent[i][j])>1.0:
            large_arg.append(j*360.0/len(datContent[i]))
        cir+=datContent[i][j]
    cir_list.append(cir)

theta = [2*3.1416*i/len(datContent[48]) for i in range(len(datContent[48]))]
plt.plot(theta, datContent[48],'.')
plt.yscale('log')
plt.savefig("Circulation_t_48.png")
print(cir_list)
print(large_arg)
