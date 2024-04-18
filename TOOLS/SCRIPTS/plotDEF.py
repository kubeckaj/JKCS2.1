import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import argparse

filename = "YOUR-DATA-FILE"
x_col = 1
y_col = 2
usecols = (x_col - 1, y_col - 1) if x_col and y_col else None
data = np.loadtxt(filename, usecols=usecols)

plt.figure()

### SCATTER/LINE/BAR/HIST- PLOT ### {:,--,-.}
plt.plot(data[:, 0], data[:, 1], linestyle=':')
#plt.scatter(data[:, 0], data[:, 1], label='Data')
#plt.bar(data[:, 0], data[:, 1])
#plt.hist(data)

### REGRESSION ###
#fit_coefficients = np.polyfit(data[:, 0], data[:, 1], 1)
#fit_function = np.poly1d(fit_coefficients)
#x_values_for_fit = np.linspace(min(data[:, 0]), max(data[:, 0]), num=100)  # Change 0 to 1 if you want to spah over y range
#plt.plot(x_values_for_fit, fit_function(x_values_for_fit), color='red', linestyle=':', label=f'Linear fit: y={fit_coefficients[0]:.2f}x+{fit_coefficients[1]:.2f}')

plt.legend()

#SCALE: {linear,log,symlog,logit}
plt.xscale('linear')
plt.yscale('linear')

#TITLE
plt.title('')

#AXES LABEL
plt.xlabel('X')
plt.ylabel('Y')

#GRID
#plt.grid()

#VISUALIZE
plt.show()
#plt.savefig("plot.svg", dpi=600)


