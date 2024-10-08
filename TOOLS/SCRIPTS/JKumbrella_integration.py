import numpy as np
import matplotlib.pyplot as plt


# Load the 2nd and 3rd columns (0-based index: column 1 and column 2)
data = np.loadtxt('out', usecols=(0,1))
data = data[~np.isnan(data).any(axis=1)]

plt.figure()
plt.plot(data[:, 0], data[:, 1], marker='.', linestyle='-', color='b', label='Column 2 vs Column 3')
plt.xlabel('coordinate')
plt.ylabel('PMF (kJ/mol)')
plt.savefig('plot.png')

data = np.column_stack((data[:, 0],+2*0.509/0.238848*np.log(data[:, 0]+0.00000001)+0.509/0.238848*np.log(4*3.14)+data[:, 1]))
data = np.column_stack((data[:, 0],0.238846*0.0433641153087705*data[:, 1]))
data = np.column_stack((data[:, 0],data[:, 1]-np.min(data[:, 1])))

min_index = np.argmin(data[:, 1])
print(min_index)

def model(w, a, b, depth):
    return (1 - np.exp(-a * (w - b)))**2 * depth
initial_guess = [1.0, 3, 3.0]

from scipy.optimize import curve_fit

popt, pcov = curve_fit(model, data[min_index:, 0], data[min_index:, 1]/0.0433641153087705, p0=initial_guess, bounds=(0, np.inf))
a_opt, b_opt, depth_opt = popt
print(a_opt, b_opt, depth_opt)

w_fit = np.linspace(min(data[min_index:, 0]), max(data[:, 0]), 100)
y_fit = model(w_fit, *popt)
plt.figure()

plt.plot(data[:, 0], data[:, 1]/0.0433641153087705, marker='.', linestyle='-', color='b', label='Column 2 vs Column 3')
plt.plot(w_fit, y_fit, color='red', label='Fitted curve')
plt.xlabel('coordinate')
plt.ylabel('PMF (kcal/mol)')
plt.savefig('plot4PiR.png')

V0 = 8.314*298.15/101325*10**30/6.022/10**23
mult=1

oo = np.column_stack((data[:, 0], np.exp(-data[:, 1]/0.0433641153087705/0.509)*4*3.14*data[:, 0]**2*(data[1, 0]-data[0, 0])))
tab = np.array([np.array([oo[i,0], -0.509*np.log(0.0000001+np.sum(np.exp(depth_opt/0.509)/V0*oo[1:i,1]/mult))]) for i in range(1, len(oo))])

plt.figure()
plt.plot(tab[:, 0], tab[:, 1], marker='.', linestyle='-', color='b', label='Column 2 vs Column 3')
plt.xlabel('coordinate')
plt.ylabel('dG (kcal/mol)')
plt.savefig('plotdG.png')

print(tab[-1, 1])

