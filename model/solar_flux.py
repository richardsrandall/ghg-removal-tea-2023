
from solcore.light_source import LightSource
import numpy as np

# Compute the total photon flux below a specified wavelength in the AM1.5g solar spectrum
n_steps = 500
lower_bound_nm = 50
upper_bound_nm = 365
wavelengths=np.linspace(lower_bound_nm,upper_bound_nm, n_steps)
am15g = LightSource(source_type='standard', x=wavelengths,
                    version='AM1.5g', output_units='photon_flux_per_nm')
integrated_photons = sum(am15g.spectrum()[1])*(upper_bound_nm-lower_bound_nm)/n_steps
solar_uv_photon_flux = integrated_photons #Measured in photons, not moles of photons.

# Convert to moles
avogadro = 6.022e23

print("")
print("Solar photon flux, in moles per (m2*second), between "+
      str(lower_bound_nm)+" and "+str(upper_bound_nm)+" nm")
print(solar_uv_photon_flux/avogadro)
