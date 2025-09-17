from schnetpack.md.data import HDF5Loader
import numpy as np
import matplotlib.pyplot as plt
from schnetpack import units as spk_units
import os
data = HDF5Loader("SIM/simulation.hdf5")
T=data.get_temperature() 
Ek=data.get_kinetic_energy() 
Ep=data.get_potential_energy()
#print(spk_units.convert_units("kJ/mol", "kcal/mol"))
Ep *= 0.2390057361376673 #spk_units.convert_units("kJ/mol", "kcal/mol")
Ek *= 0.2390057361376673 #spk_units.convert_units("kJ/mol", "kcal/mol")
Ep -= np.min(Ep)
time_axis = np.arange(data.entries) * data.time_step / spk_units.fs  # in fs

if not os.path.exists("movie.xyz"):
  atoms=data.convert_to_atoms()
  print(atoms)
  (i.wrap() for i in atoms)
  from ase.io import write
  f = open("movie.xyz","w")
  for i in range(0,len(atoms),int(np.round(len(atoms)/100))):
    write(".movie.xyz",atoms[i])
    with open(".movie.xyz", 'r') as f2:
      lines = f2.readlines()
      lines[1] = str(i)+"\n"
    f.writelines(lines)
    f2.close()
  f.close()

if not os.path.exists("step_T_Ek_Ep_Etot.txt"):
  with open("step_T_Ek_Ep_Etot.txt","a") as f:
    for i in range(len(T)):
      print(f"%i %.2f %.4f %.4f %.4f" % (i,T[i],Ek[i],Ep[i],Ek[i]+Ep[i]), file=f)

########################
########################
########################
from schnetpack.md.data import PowerSpectrum #, IRSpectrum, RamanSpectrum
# Initialize the spectrum
equilibrated_data = HDF5Loader("SIM/simulation.hdf5", skip_initial=100)
Pspectrum = PowerSpectrum(equilibrated_data, resolution=2048)
#IRspectrum = IRSpectrum(data, resolution=2048)
#Rspectrum = RamanSpectrum(data, resolution=2048, incident_frequency=4000,temperature=300)
# Compute the spectrum for the first molecule (default)
Pspectrum.compute_spectrum(molecule_idx=0)     
#IRspectrum.compute_spectrum(molecule_idx=0)     
#Rspectrum.compute_spectrum(molecule_idx=0)     
# Get frequencies and intensities
Pfrequencies, Pintensities = Pspectrum.get_spectrum()
#IRfrequencies, IRintensities = IRspectrum.get_spectrum()
#Rfrequencies, Rintensities = RSspectrum.get_spectrum()

def plot_temperature(mplt):
    temperature_mean = np.cumsum(T) / (np.arange(data.entries)+1)
    #plt.figure(figsize=(8,4))
    mplt.plot(time_axis, T, label='T')
    mplt.plot(time_axis, temperature_mean, label='T (avg.)')
    #plt.ylabel('T [K]')
    #plt.xlabel('t [fs]')
    mplt.legend()
    #plt.tight_layout()
    #plt.show()

def plot_energies(mplt):
    #plt.figure()
    mplt.plot(time_axis, Ep, label="potential", ls=":")
    mplt.plot(time_axis, Ek, label="kinetic", ls="--")
    mplt.plot(time_axis, Ek+Ep, label="total")
    #plt.ylabel("E [kcal/mol]")
    #plt.xlabel("t [fs]")
    #plt.xlim(9800,10000)
    mplt.legend()
    #plt.tight_layout()
    #plt.show()

def plot_spectrum(mplt, freq, ins):
    #plt.figure()
    mplt.plot(freq, ins)
    mplt.set_xlim(0,4500)
    mplt.set_ylim(0,100)
    #plt.ylabel('I [a.u.]')
    #plt.xlabel('$\omega$ [cm$^{-1}$]')
    #plt.show() 

fig = plt.figure(figsize=(12,6))
gspec = fig.add_gridspec(ncols=3, nrows=2, width_ratios = [1]*3, height_ratios = [1]*2)
plot_temperature(fig.add_subplot(gspec[0,0]))
plot_energies(fig.add_subplot(gspec[0,1:]))
plot_spectrum(fig.add_subplot(gspec[1,:]),Pfrequencies,Pintensities)
#plot_spectrum(fig.add_subplot(gspec[1,1]),IRfrequencies,IRintensities)
#plot_spectrum(fig.add_subplot(gspec[1,2]),Rfrequencies,Rintensities)
plt.tight_layout()
mng = plt.get_current_fig_manager()
mng.resize(*mng.window.maxsize())
plt.show()
#plt.tight_layout(pad=1.5, h_pad=2.5, w_pad=2)


