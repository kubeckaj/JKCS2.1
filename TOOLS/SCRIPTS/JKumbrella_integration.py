import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.linalg import eigh

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

print("")
print("Fitted paramteres of function: (1 - np.exp(-a * (w - b)))**2 * depth")
print("# a b depth")
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
mult=int(sys.argv[1])

oo = np.column_stack((data[:, 0], np.exp(-data[:, 1]/0.0433641153087705/0.509)*4*3.14*data[:, 0]**2*(data[1, 0]-data[0, 0])))
tab = np.array([np.array([oo[i,0], -0.509*np.log(0.0000001+np.sum(np.exp(depth_opt/0.509)/V0*oo[1:i,1]/mult))]) for i in range(1, len(oo))])

plt.figure()
plt.plot(tab[:, 0], tab[:, 1], marker='.', linestyle='-', color='b', label='Column 2 vs Column 3')
plt.xlabel('coordinate')
plt.ylabel('dG (kcal/mol)')
plt.savefig('plotdG.png')

print("")
print("dG =")
print(tab[-1, 1])

###########################################
#####  SOLVIND SCHRODINGER 1D #############
###########################################

try:
  M1=float(sys.argv[2])
  M2=float(sys.argv[3])
except:
  print("QC not possible to solve. I need commands_TODO.txt file for it")
  exit()

x_angstrom = data[:,0]
V_kjmol = data[:,1]/0.0433641153087705/0.238846

# Constants for unit conversion
angstrom_to_bohr = 1.88973           # 1 Å = 1.88973 Bohr radii
kjmol_to_hartree = 0.0003808798      # 1 kJ/mol = 0.0003808798 Hartree

# Convert units to atomic units (a.u.)
x = x_angstrom * angstrom_to_bohr  # Convert Å to Bohr radii
V = V_kjmol * kjmol_to_hartree     # Convert kJ/mol to Hartrees

# Discretization
dx = x[1] - x[0]  # Grid spacing (assumes uniform grid)
N = len(x)        # Number of grid points

# Potential minimum
V_min = np.min(V)

# Construct the kinetic energy operator using finite differences
#T = np.zeros((N, N))
#for i in range(N):
#    if i > 0:
#        T[i, i - 1] = -1
#    if i < N - 1:
#        T[i, i + 1] = -1
#    T[i, i] = 2
#T *= -1 / (2 * dx**2)  # Scale by -ħ²/2m (ħ=1, m=1 in atomic units)

# Construct the kinetic energy operator using the correct finite difference method
T = np.zeros((N, N))
for i in range(1, N-1):  # Avoid boundary points for simplicity
    T[i, i-1] = -0.5 / dx**2
    T[i, i+1] = -0.5 / dx**2
    T[i, i] = 1.0 / dx**2

# Boundary conditions: assume zero derivative at the boundaries
T[0, 0] = T[-1, -1] = 1.0 / dx**2

m = M1*M2/(M1+M2)/(1000*9.1093837*10**-31*6.022*10**23)
T *= 1 / (m)

# Construct the potential energy operator as a diagonal matrix
V_matrix = np.diag(V)

# Hamiltonian is the sum of kinetic and potential energy operators
H = T + V_matrix

# Solve the eigenvalue problem
eigenvalues, eigenvectors = eigh(H)

# Ground state energy (ZPE)
E_ground = eigenvalues[0]

# Difference between ZPE and potential minimum
delta_E = E_ground - V_min

# Print results
print("")
print(f"Zero-point energy (in Hartrees): {E_ground:.6f}")

# Convert to kJ/mol for clarity
E_ground_kjmol = E_ground / kjmol_to_hartree
E_ground_kcalmol = E_ground_kjmol * 0.239006
V_min_kjmol = V_min / kjmol_to_hartree
delta_E_kjmol = delta_E / kjmol_to_hartree

print(f"Zero-point energy (in kJ/mol): {E_ground_kjmol:.6f}")
print(f"Zero-point energy (in kcal/mol): {E_ground_kcalmol:.6f}")

print("dG_qc-corr: [kcal/mol]")
print(E_ground_kcalmol+tab[-1, 1])

#Aho = 0.509*np.log(2*E_ground_kcalmol/0.509)
#Aho = 0.509*np.log(1-np.exp(-2*E_ground_kcalmol/0.509))
#Aqo = 0.509*np.log(0.5*np.sinh(E_ground_kcalmol/0.509))
#Aqo = -0.509*np.log(0.5*np.sinh(E_ground_kcalmol/0.509)**-1)+E_ground_kcalmol
#print("")
#print("HO =",Aho)
#print("QO =",Aqo)
#print("QO-HO =",-Aho+Aqo)
#print(tab[-1, 1]+Aqo-Aho)
