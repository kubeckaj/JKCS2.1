import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.linalg import eigh

# 0.238846 kJ/mol -> kcal/mol
# 0.0433641153087705 kcal/mol -> eV

print(" ----- Post-Processing ----- ")

data = np.loadtxt('out', usecols=(0,1)) #Ang kJ/mol
data = data[~np.isnan(data).any(axis=1)]

plt.figure()
plt.plot(data[:, 0], data[:, 1], marker='.', linestyle='-', color='b', label='Column 2 vs Column 3')
plt.xlabel('coordinate')
plt.ylabel('PMF (kJ/mol)')
plt.savefig('PMF_3D_kJmol-1.png')

#from umbrellaintegration
print("PIC: Histograms.png created.")

### 3D PMF ###
data = np.loadtxt('out', usecols=(0,1)) #Ang kJ/mol
data = data[~np.isnan(data).any(axis=1)]
PMF_3D = np.column_stack((data[:, 0], 0.238846*0.0433641153087705*data[:, 1])) #Ang eV

plt.figure()
plt.plot(PMF_3D[:, 0], PMF_3D[:, 1], marker='.', linestyle='-', color='b', label='Column 2 vs Column 3')
plt.xlabel('coordinate')
plt.ylabel('PMF (kcal/mol)')
plt.savefig('PMF_3D_kcalmol-1.png')
print("PIC: PMF_3D_kcalmol-1.png created.")

### 1D PMF ###
PMF_1D = np.column_stack((PMF_3D[:, 0], +2*0.592*0.0433641153087705*np.log(PMF_3D[:, 0])+PMF_3D[:, 1]))
PMF_1D = np.column_stack((PMF_1D[:, 0], PMF_1D[:, 1]-np.min(PMF_1D[:, 1])))

plt.figure()
plt.plot(PMF_1D[:, 0], PMF_1D[:, 1]/0.0433641153087705, marker='.', linestyle='-', color='b', label='Column 2 vs Column 3')
plt.xlabel('coordinate')
plt.ylabel('PMF (kcal/mol)')
plt.savefig('PMF_1D_kcalmol-1.png')
print("PIC: PMF_1D_kcalmol-1.png created.")

min_index = np.argmin(PMF_1D[:, 1]) #index
print("MIN: Position of minimun has index: "+str(min_index))

if 0==int(sys.argv[2]):
  def model(w, a, b, depth):
      return (1 - np.exp(-a * (w - b)))**2 * depth
  initial_guess = [1.0, PMF_1D[min_index, 0], 3.0]
  
  from scipy.optimize import curve_fit
  popt, pcov = curve_fit(model, PMF_1D[min_index:, 0], PMF_1D[min_index:, 1]/0.0433641153087705, p0=initial_guess, bounds=(0, np.inf))
  a_opt, b_opt, depth_opt = popt #depth is in kcal/mol
  print("FIT:Fitting parameters of function: (1 - np.exp(-a * (w - b)))**2 * depth")
  print("FIT: a b depth = "+str(a_opt)+" "+str(b_opt)+" "+str(depth_opt))

  w_fit = np.linspace(min(PMF_1D[min_index:, 0]), max(PMF_1D[:, 0]), 100)
  y_fit = model(w_fit, *popt)
else:
  max_index = np.argmax(PMF_3D[min_index:, 1]) #index
  def model(w,a):
      return a
  initial_guess = PMF_1D[(min_index+max_index), 1]/0.0433641153087705
  from scipy.optimize import curve_fit
  popt, pcov = curve_fit(model, PMF_1D[min_index+max_index:, 0], PMF_1D[min_index+max_index:, 1]/0.0433641153087705, p0=initial_guess, bounds=(0, np.inf))
  depth_opt = popt[0] #depth is in kcal/mol
  print("FIT:Fitting parameters of function: depth")
  print("FIT: depth = "+str(depth_opt))
  w_fit = np.linspace(min(PMF_1D[(min_index+max_index):, 0]), max(PMF_1D[:, 0]), 100)
  y_fit = model(w_fit, np.array([depth_opt]*100))

plt.figure()
plt.plot(PMF_1D[:, 0], PMF_1D[:, 1]/0.0433641153087705, marker='.', linestyle='-', color='b', label='Column 2 vs Column 3')
plt.plot(w_fit, y_fit, color='red', label='Fitted curve')
plt.xlabel('coordinate')
plt.ylabel('PMF (kcal/mol)')
plt.savefig('FITTED.png')
print("PIC: FITTED.png created.")

V0 = 8.314*298.15/101325*10**30/6.022/10**23 #Ang^3
mult=int(sys.argv[1]) #multiplier, should be 2 for symmetric reactions, otherwise 1

ddx = PMF_1D[1, 0] - PMF_1D[0, 0]
oo = np.column_stack((PMF_1D[:, 0], np.exp(-PMF_1D[:, 1]/0.0433641153087705/0.592)*4*3.14*PMF_1D[:, 0]**2*ddx))
tab = np.array([np.array([oo[i,0], -0.592*np.log(0.000000001+np.sum(np.exp(depth_opt/0.592)/V0*oo[1:i,1]/mult))]) for i in range(1, len(oo))])

plt.figure()
plt.plot(tab[:, 0], tab[:, 1], marker='.', linestyle='-', color='b', label='Column 2 vs Column 3')
plt.xlabel('coordinate')
plt.ylabel('dG (kcal/mol)')
plt.savefig('dG_convergence.png')
print("PIC: dG_convergence.png created.")

print("\ndG = "+str(tab[-1, 1])+" kcal/mol\n")

###########################################
#####  SOLVIND SCHRODINGER 1D #############
###########################################

try:
  M1=float(sys.argv[3])
  M2=float(sys.argv[4])
except:
  print("QC not possible to solve. I need commands_TODO.txt file for it")
  exit()

x_angstrom = PMF_1D[:,0]
V_kjmol = PMF_1D[:,1]/0.0433641153087705/0.238846

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
