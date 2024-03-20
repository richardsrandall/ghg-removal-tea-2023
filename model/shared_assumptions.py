
# Contains assumptions shared between the models to make sure
# there aren't conflicting versions in different places

from solcore.light_source import LightSource
import numpy as np

# Mostly used in constructors for parameter objects

class shared_assumptions:

    def __init__(self,which_GHG):

        # Air properties
        self.air_kin_visc = 1.48e-5 #Kinematic viscosity
        self.rho_air = 1.225 #Density/kg/m3

        # Properties of specific GHGs
        self.CH4_molar_mass = 0.016042 #kg/mole
        self.CH4_moles_per_ton = 1e6/self.CH4_molar_mass #moles/ton
        self.CH4_diffusivity = 2.21e-5 #m2/s
        self.CH4_schmidt_number = self.air_kin_visc/self.CH4_diffusivity #dimensionless

        self.N2O_molar_mass = 0.044013 #kg/mole
        self.N2O_moles_per_ton = 1e6/self.N2O_molar_mass #moles/ton
        self.N2O_diffusivity = 1.43e-5 #m2/s
        self.N2O_schmidt_number = self.air_kin_visc/self.N2O_diffusivity #dimensionless

        # Values that depend on which GHG is specified
        if which_GHG == 'CH4':
            self.ghg_name = 'CH4'
            self.ghg_molar_mass = self.CH4_molar_mass
            self.ghg_moles_per_ton = self.CH4_moles_per_ton
            self.ghg_diffusivity = self.CH4_diffusivity
            self.schmidt_number = self.CH4_schmidt_number
        elif which_GHG == 'N2O':
            self.ghg_name = 'N2O'
            self.ghg_molar_mass = self.N2O_molar_mass
            self.ghg_moles_per_ton = self.N2O_moles_per_ton
            self.ghg_diffusivity = self.N2O_diffusivity
            self.schmidt_number = self.N2O_schmidt_number
        elif which_GHG == 'Test':
            pass #This exists for when we just want to print the contents
                    # of the shared_assumptions object for testing
        else:
            raise Exception("Invalid GHG selected.")

        # Unit conversions
        self.avogadro = 6.022e23
        self.ev_per_j = 6.24e18 #eV per Joule
        self.j_per_kwh = 3.6e6
        self.seconds_per_year = 3600*24*365 # Seconds per year
        self.seconds_per_day = 3600*24
        self.m_per_mm = 0.001 #meters/mm
        self.hours_per_year = 8760

        # Some other shared assumptions
        self.uv_photon_max_wavelength = 365 #nm; depends on the wavelength at which photocatalyst was studied
        self.CRF = 0.075
        self.temp = 300 #K
        self.air_molar_volume = 0.0246 #m3/mole, at 1 bar and 300K
        self.moles_per_m3_per_ppm = (1/self.air_molar_volume)*(1e-6) #(moles/m3) per ppm
        
        # Photocatalyst layer and support properties
        self.photocatalyst_thickness = 2 # microns
        self.photocatalyst_density = 5600 # kg/m3; ZnO
        self.support_density = 2650 # kg/m3; silica
        self.support_cost_per_ton = 500
        self.photocatalyst_cost_per_ton = 4500
        
        # Compute the total photon flux below a specified wavelength in the AM1.5g solar spectrum
        n_steps = 500
        lower_bound_nm = 50
        wavelengths=np.linspace(lower_bound_nm,self.uv_photon_max_wavelength, n_steps)
        am15g = LightSource(source_type='standard', x=wavelengths, version='AM1.5g', output_units='photon_flux_per_nm')
        integrated_photons = sum(am15g.spectrum()[1])*(self.uv_photon_max_wavelength-lower_bound_nm)/n_steps
        self.solar_uv_photon_flux = integrated_photons #Measured in photons, not moles of photons.

        # Physical constants
        self.k = 8.314 #J/K-mol; Boltzmann constant / universal gas constant
        
        

