import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.linalg import eigh

# All is in Hartree and Agstrom if not mentioned otherwise in the variable

## UNITS ##
kB = 0.000003166811563 #Eh/K
patm = 101325 #Pa
NA = 6.022*10**23 #mol^-1
pi = 3.14159265359
## UNIT CONVERSIONS ##
Eh2kcalmol = 627.5094740631
Eh2kJmol = 2625.4996394799
Eh2J = 4.3597447222071*10**-18
Ang2Bohr = 1.88973
amu2au = 1822.888486209
print(" ----- Post-Processing ----- ")

### ARGUMENTS ###

Qfit = int(sys.argv[2])
Qflatten = float(sys.argv[3])
T = float(sys.argv[6]) #K
kBT = T*kB #Eh

### PREPARING DATA ###
data = np.loadtxt('out', usecols=(0,1)) #Ang kJ/mol
data = data[~np.isnan(data).any(axis=1)]
data[:, 1] /= Eh2kJmol

print("PIC: Histograms.png created.")

### 3D PMF ###
PMF_3D = np.column_stack((data[:, 0], data[:, 1])) #Ang Eh

plt.figure()
plt.plot(PMF_3D[:, 0], PMF_3D[:, 1]*Eh2kcalmol, marker='.', linestyle='-', color='b')
plt.xlabel('coordinate')
plt.ylabel('PMF (kcal/mol)')
plt.savefig('PMF_3D_kcalmol-1.png')
print("PIC: PMF_3D_kcalmol-1.png created.")

### 1D PMF ###
PMF_1D = np.column_stack((PMF_3D[:, 0], PMF_3D[:, 1]+ 2*kBT*np.log(PMF_3D[:, 0]))) #Ang Eh
PMF_1D = np.column_stack((PMF_1D[:, 0], PMF_1D[:, 1]-np.min(PMF_1D[:, 1]))) #Ang Eh

PMF_1D_orginal = np.copy(PMF_1D)
if Qflatten > 0:
  mask=PMF_1D[:, 0] <= Qflatten
  PMF_1D[mask, 1] = (PMF_1D[mask, 1])[-1] 
 
plt.figure()
plt.plot(PMF_1D[:, 0], PMF_1D_orginal[:, 1]*Eh2kcalmol, marker='.', linestyle='-', color='b')
#TODO remove 2 lines
#plt.plot(PMF_3D[:, 0], PMF_3D[:, 1]*Eh2kcalmol, marker='.', linestyle='-', color='g')
#plt.plot(PMF_3D[:, 0], 2*kBT*np.log(PMF_3D[:, 0])*Eh2kcalmol, marker='.', linestyle='-', color='r')
if Qflatten > 0:
  plt.plot(PMF_1D[mask, 0], PMF_1D[mask, 1]*Eh2kcalmol, marker='.', linestyle='-', color='g')
plt.xlabel('coordinate')
plt.ylabel('PMF (kcal/mol)')
plt.savefig('PMF_1D_kcalmol-1.png')
print("PIC: PMF_1D_kcalmol-1.png created.")

### MINIMUM INDEX ###
min_index = np.argmin(PMF_1D[:, 1]) #index
print("MIN: Position of minimun has index: "+str(min_index))

### FITTING ###
max_index = min_index + np.argmax(PMF_3D[min_index:, 1]) #index
print("MAX: Position of maximum has index: "+str(max_index))
print("Max = ", PMF_1D[max_index, 0]," Angstrom")
if Qfit==0:
  def model(w, a, xe, De):
      return (1 - np.exp(-a * (w - xe)))**2 * De
  initial_guess = [1.0, PMF_1D[min_index, 0], 0.001]
  
  from scipy.optimize import curve_fit
  popt, pcov = curve_fit(model, PMF_1D[min_index:, 0], PMF_1D[min_index:, 1], p0=initial_guess, bounds=(0, np.inf))
  a, xe, De = popt #depth is in Eh
  print("FIT:Fitting parameters of function: (1 - np.exp(-a * (w - b)))**2 * depth")
  print("FIT: a[1/Ang] xe[Ang] De[kcal/mol] = "+str(a)+" "+str(xe)+" "+str(De*Eh2kcalmol))

  w_fit = np.linspace(PMF_1D[min_index, 0], PMF_1D[-1, 0], 100)
  y_fit = model(w_fit, *popt)
else:
  def model(w,De):
      return De
  initial_guess = PMF_1D[(max_index), 1] #Eh
  from scipy.optimize import curve_fit
  popt, pcov = curve_fit(model, PMF_1D[max_index:, 0], PMF_1D[max_index:, 1], p0=initial_guess, bounds=(0, np.inf))
  De = popt[0] #Eh
  print("FIT:Fitting parameters of function: 1*De")
  print("FIT: De[kcal/mol] = "+str(De*Eh2kcalmol))
  w_fit = np.linspace(PMF_1D[max_index, 0], PMF_1D[-1, 0], 100)
  y_fit = model(w_fit, np.array([De]*100))

plt.figure()
plt.plot(PMF_1D[:, 0], PMF_1D_orginal[:, 1]*Eh2kcalmol, marker='.', linestyle='-', color='b')
if Qflatten > 0:
  plt.plot(PMF_1D[mask, 0], PMF_1D[mask, 1]*Eh2kcalmol, marker='.', linestyle='-', color='g')
plt.plot(w_fit, y_fit*Eh2kcalmol, color='red')
plt.xlabel('coordinate')
plt.ylabel('PMF (kcal/mol)')
plt.savefig('FITTED.png')
print("PIC: FITTED.png created.")

try:
  M1=float(sys.argv[4])
  M2=float(sys.argv[5])
except:
  print("QC not possible to solve. I need commands_TODO.txt file for it")
  exit()

### V0 ###
V0 = kBT*Eh2J/patm * 10**30 #Ang^3
print("V0 = "+str(V0)+" Ang^3")
#print("L0 = "+str((3/4/3.14*V0)**(1/3))+" Ang")

### MULTIPLIER ###
mult=float(sys.argv[1]) #multiplier, should be 2 for symmetric reactions, otherwise 1
print("Multiplier = "+str(mult)+" [1 or 2 for symmetric reactions]")

### INTEGRATION/SUM ###
dx = PMF_1D[1, 0] - PMF_1D[0, 0]
sum_element = np.column_stack((PMF_1D[:, 0], np.exp(-PMF_1D[:, 1]/kBT)*4*pi*PMF_1D[:, 0]**2*dx))
#sum_element_ivo = np.column_stack((PMF_1D[:, 0], np.exp(-PMF_1D[:, 1]/kBT)*PMF_1D[:, 0]**2))
#sum_table_ivo = np.array([np.array([sum_element[i,0], np.sum(sum_element_ivo[i,1])]) for i in range(0, len(sum_element))])
sum_table = np.array([np.array([sum_element[i,0], np.sum(sum_element[0:i,1])]) for i in range(0, len(sum_element))])
F_well = sum_table[-1, 1] #Eh
F_table = -kBT*np.log(np.exp(De/kBT)*sum_table[:, 1]/V0/mult) 
F = F_table[max_index] #Eh

plt.figure()
plt.plot(sum_table[:, 0], F_table[:]*Eh2kcalmol, marker='.', linestyle='-', color='b')
#plt.plot(sum_table[:, 0], sum_table_ivo[:,1], marker='.', linestyle='-', color='b')
#print((sum_element[:, 1]-sum_element[-1, 1]))
#print((np.max(sum_element[:, 1])-sum_element[-1, 1]))
#print(np.max(F_table[0])-F_table[-1])
plt.plot(sum_element[:, 0], ((sum_element[:, 1]-sum_element[-1, 1])/(np.max(sum_element[:, 1])-sum_element[-1, 1])*(F_table[1]-F_table[max_index])*Eh2kcalmol+F_table[max_index]*Eh2kcalmol), marker='.', linestyle='-', color='r')
plt.xlabel('coordinate')
plt.ylabel('dF (kcal/mol)')
plt.savefig('dF_convergence.png')
print("PIC: dF_convergence.png created.")

### PTINTING CLASSICAL RESAULT ###
print("RESUTLS (CLASSICAL):")
print("\ndF = "+str(F*Eh2kcalmol)+" kcal/mol")
print("\ndG = "+str((F-kBT)*Eh2kcalmol)+" kcal/mol\n")

###########################################
#####  SOLVIND SCHRODINGER 1D #############
###########################################

# Convert units to atomic units (a.u.)
x = PMF_1D[:,0]*Ang2Bohr # Bohr
dx = x[1] - x[0]
V = PMF_1D[:,1] #Eh

N = len(x)   # Number of grid points

# Construct the kinetic energy operator using the correct finite difference method
T = np.zeros((N, N))
for i in range(1, N-1):  # Avoid boundary points for simplicity
    T[i, i-1] = -0.5 / dx**2
    T[i, i+1] = -0.5 / dx**2
    T[i, i] = 1.0 / dx**2

# Boundary conditions: assume zero derivative at the boundaries
T[0, 0] = 1.0 / dx**2
T[-1, -1] = 1.0 / dx**2
#Dirichlet, WF vanishes in boundaries
#T[0, :] = 0
#T[-1, :] = 0

m = M1*M2/(M1+M2)*amu2au
T *= 1 / m

# Construct the potential energy operator as a diagonal matrix
V_matrix = np.diag(V)

# Hamiltonian is the sum of kinetic and potential energy operators
H = T + V_matrix

# Solve the eigenvalue problem
eigenvalues, eigenvectors = eigh(H)
eigenvalues = eigenvalues[eigenvalues < De]

# Ground state energy (ZPE)
ZPE = eigenvalues[0]
print("Energy levels: ", eigenvalues*Eh2kcalmol, " kcal/mol")
print(f"Zero-point energy = {ZPE:.6f} Eh = {ZPE*Eh2kcalmol:.6f} kcal/mol")
print("")

print(" --Approximation by quantum and classical Morse wells--")
### FREE ENERGY OF APPROXIMATED QUANTUM WELL ###
#Z_qmw = np.exp(ZPE/kBT)*(1-np.exp(-2*ZPE/kBT))
Z_qmw = 2*np.sinh(ZPE/kBT)
F_qmw = kBT * np.log(Z_qmw)
print("Quantum well dF = ", F_qmw, "Hartree =", F_qmw*Eh2kcalmol, " kcal/mol")

#### FREE ENERGY OF APPROXIMATED CLASSICAL WELL ###
#Z_cmw = np.exp(-De/kBT)*pi*kBT/a*np.sqrt(m/2/De)
Z_cmw = 2*ZPE/kBT
F_cmw = kBT * np.log(Z_cmw)
print("Classical well dF = ", F_cmw, "Hartree =", F_cmw*Eh2kcalmol, " kcal/mol")
print("ddG = ddF = ", (F_qmw-F_cmw)*Eh2kcalmol, "kcal/mol")

print("")
print("RESUTLS:")
print("dF_qc-corr = ", (F+F_qmw-F_cmw)*Eh2kcalmol, " kcal/mol")
print("dG_qc-corr = ", (F-kBT+F_qmw-F_cmw)*Eh2kcalmol, " kcal/mol")

