import meep as mp
import meep.materials as mt
import numpy as np
import matplotlib.pyplot as plt

Size_x=2
Size_y=2
Size_z=2 #setting up simulation box
cell_size=mp.Vector3(Size_x,Size_y,Size_z)

#adding symmetry:
symmetries = [ mp.Mirror(direction=mp.X, phase=+1), mp.Mirror(direction=mp.Y, phase=+1)]


center_wavelength = 0.7
min_wavelength=0.4
max_wavelength=1.0
center_freq=1.0/center_wavelength
freq_width=1.0/min_wavelength - 1.0/max_wavelength


source_position=Size_z*0.5-0.3
gaussian_pulse=mp.GaussianSource(frequency=center_freq,fwidth=2*freq_width)
source=[mp.Source(gaussian_pulse,component=mp.Ex,center=mp.Vector3(0,0,source_position),size=(Size_x,Size_y,0))]


Ag_top=mp.Block(material=mt.Ag,center=mp.Vector3(0,0,.075),size=mp.Vector3(1.2,1.2,0.035))
Ag_bottom=mp.Block(material=mt.Ag,center=mp.Vector3(0,0,-0.0525),size=mp.Vector3(1.2,1.2,0.1))
TiO2_index=2.4667
TiO2_material=mp.Medium(index=TiO2_index)
InSb_material=mp.Medium(epsilon=16.8,E_susceptibilities=[mp.DrudeSusceptibility(frequency=0.1167, gamma=0.0011)])
TiO2=mp.Block(material=TiO2_material,center=mp.Vector3(0,0,0.0175),size=mp.Vector3(1.2,1.2,0.04))
InSb=mp.Block(material=InSb_material,center=mp.Vector3(0,0,0.0475),size=mp.Vector3(1.2,1.2,0.02))
geometry=[Ag_top,InSb,TiO2,Ag_bottom]#Ag_top

resolution=100
pmls=[mp.PML(thickness=0.2)]

sim=mp.Simulation(resolution=resolution,cell_size=cell_size,sources=source,boundary_layers=pmls,geometry=[],symmetries=symmetries)

wavelengths = np.linspace(min_wavelength,max_wavelength,120)
frequencies=1.0/wavelengths

ref1=sim.add_flux(frequencies,mp.FluxRegion(center=mp.Vector3(0,0,0.2),size=mp.Vector3(1.2,1.2,0)))

sim.run(until=10) ##test run
#sim.run(until_after_sources=mp.stop_when_fields_decayed(5,mp.Ex,mp.Vector3(0,0,0),1e-3))

incident_flux = mp.get_fluxes(ref1)
refl_data = sim.get_flux_data(ref1)


incident_flux = np.array(incident_flux)

# Save to reuse later
#np.save('incident_flux.npy', incident_flux)
#sim.save_flux('refl_data', refl_data)

print("incident Flux Array:", incident_flux)
print("Shape:", np.shape(incident_flux))

sim.reset_meep()
sim=mp.Simulation(resolution=resolution,cell_size=cell_size,sources=source,boundary_layers=pmls,geometry=geometry,symmetries=symmetries)


ref2=sim.add_flux(frequencies,mp.FluxRegion(center=mp.Vector3(0,0,0.2),size=mp.Vector3(1.2,1.2,0)))
trans2=sim.add_flux(frequencies,mp.FluxRegion(center=mp.Vector3(0,0,-0.2),size=mp.Vector3(1.2,1.2,0)))

sim.load_minus_flux_data(ref2, refl_data)

sim.run(until=10) ##test run
#sim.run(until_after_sources=mp.stop_when_fields_decayed(5,mp.Ex,mp.Vector3(0,0,0),1e-3))

#refl_data = np.load('refl_data.npy', allow_pickle=True).item()


reflected_flux = np.array(mp.get_fluxes(ref2))
transmitted_flux=np.array(mp.get_fluxes(trans2))


print("Reflected Flux Array:", reflected_flux)
print("Shape:", np.shape(reflected_flux))

reflectance=np.absolute(reflected_flux)/np.absolute(incident_flux)
transmittance=np.absolute(transmitted_flux)/np.absolute(incident_flux)


absorption=1-reflectance

print("Transmitted Flux Array:", transmitted_flux)
print("Shape:", np.shape(transmitted_flux))

print("\n\nAbsorption fractional array:",absorption)


plt.plot(wavelengths[30:70], reflectance[30:70], label="REFLECTION")
plt.plot(wavelengths[30:70], absorption[30:70], '--', label="ABSORPTION")
plt.xlabel("Wavelenght(micron)")
plt.ylabel("Absorption and reflection")
plt.legend()
plt.savefig("InSb_0V_symmetrized_140res.jpeg")
plt.show()
plt.close()
#plt.close()
