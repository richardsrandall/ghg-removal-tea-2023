
from numpy import interp
from math import log, pi

import model.shared_assumptions as sa

# Script for modeling GHG removal with tropospheric aerosols

# Object to hold assumptions and inputs for single runs of the script
class aerosol_model_params(sa.shared_assumptions):

    # Constructor for the parameter object
    def __init__(self, which_GHG):
        super().__init__(which_GHG)

        # All shared assumptions are stored in the "shared assumptions" module.
        
        # Independent variables that must be explicitly set
        self.altitude = None
        self.particle_diameter = None
        self.AQY = None

        # Values that depend on which GHG is specified
        self.which_GHG = which_GHG
        if which_GHG == 'CH4':
            self.start_ppm = 1.9 #Always run this model on 1.9 ppm CH4 unless overridden
        elif which_GHG == 'N2O':
            self.start_ppm = 0.33 #Always run this model on 0.33 ppm N2O onless overridden
        elif which_GHG == 'Test':
            self.start_ppm = None #This is just for when you're printing a table of the model's assumptions
        else:
            raise Exception("Invalid GHG specified")

        # Assumptions unique to the aerosol model
        self.GHG_mean_free_path = 67e-9 # Particle mean free path of 68 nanometers in air
        self.capacity_factor = 0.35 #Sunlight for 1/4 of the day
        
        # Aviation costs and metrics
        self.cropduster_cost_per_ton_launched = 580
        
        # Constants used in Jaenicke's formula for aerosol residence times in the atmosphere
        self.residence_time_R = 0.3 #microns -- this is the radius of the particle size with the longest residence time
        self.residence_time_K = 2*6.4e7 #seconds -- 2x the duration 0.3um particles would last in the absence of wet deposition
        self.tau_wet_lower_troposphere = 8*self.seconds_per_day #Max residence time 8 days in lower troposphere due to wet deposition
        self.tau_wet_upper_troposphere = 21*self.seconds_per_day #Corresponding value is 21 days in the upper troposphere
        
        # Flag whether or not to print intermediate values for debugging
        self.print_tests = False
        

# Calculate GHG removed per particle per second
def get_ghg_removal_rate_per_particle(params,outputs):
    
    # Calculate Knudsen number, which is used in our mass transport formula
    kn = params.GHG_mean_free_path/(0.5*params.particle_diameter)

    # Calculate average particle velocity
    c_a = ((8*params.k*params.temp)/(3.14159*params.ghg_molar_mass))**0.5
    # Unit check: (J/K-mol)*(K)/(kg/mol) = J/kg = m2/s2... it checks out.

    # Get particle surface area
    radius_m = 0.5*params.particle_diameter*(1e-6)
    particle_sa = radius_m*radius_m*3.14159*4
        
    # Calculate the mass transfer rate to the particle
    # These equations and the cleaned-up ones in the SI are from: https://authors.library.caltech.edu/25069/7/AirPollution88-Ch5.pdf
    # It's assumed that the particle rotates quickly and that depletion/mass transfer happen at the same rate all around the particle
    j_ac = 4*3.14159*radius_m*params.ghg_diffusivity*params.start_ppm*params.moles_per_m3_per_ppm # Formula for a big particle relative to the MFP and mass transport boundary layer
    b_f = (1+kn)/(1+((4*params.ghg_diffusivity)/(c_a*params.GHG_mean_free_path))*kn*(1+kn)) # Adjustment for local depletion of particle near surface
    j_a = j_ac*b_f
    moles_transport_per_particle_second = j_a

    # Account for photon supply
    # It's assumed that the particle rotates quickly, so the whole surface gets an even photon flux on average 25% the incident light intensity
    # That's because (frontal area) / (surface area) = (pi*r*r)/(4*pi*r*r) = 1/4.
    AQY_limit_moles = params.AQY*params.solar_uv_photon_flux/params.avogadro
    moles_photons_per_particle_second = 0.25*AQY_limit_moles*particle_sa

    # Determine which of them is limiting
    moles_removed_per_particle_second = min(moles_transport_per_particle_second, moles_photons_per_particle_second)
    # Kludgily, I temporarily store some intermediate values (used by later functions) in the params object instead of returning them
    outputs["Limiting Factor"] = "Photon Flux" if moles_transport_per_particle_second>moles_photons_per_particle_second else "Mass Transport"
    outputs["Removal Rate Per Particle"] = moles_removed_per_particle_second
    # Return output
    return moles_removed_per_particle_second*params.capacity_factor


# Calculate GHG removed per ton particle per second
def get_ghg_removal_rate_per_ton(params,outputs):
    # Find the removal rate per particle per second
    rate_per_particle = get_ghg_removal_rate_per_particle(params,outputs)
    # Compute the mass of a particle consisting of a catalyst shell and an inert core (or pure catalyst for small particles)
    core_volume = pi*max(0,(4/3)*(((0.5*params.particle_diameter*1e-6)-params.photocatalyst_thickness*1e-6)**3))
    total_volume = (4/3)*pi*((0.5*params.particle_diameter*1e-6)**3)
    shell_volume = total_volume-core_volume
    core_mass = core_volume*params.support_density
    shell_mass = shell_volume*params.photocatalyst_density
    particle_mass_tons = 0.001*(core_mass+shell_mass)
    # Combine to get the removal rate per ton per second
    rate_per_ton = rate_per_particle/particle_mass_tons
    return rate_per_ton

# Calculate lifespan in seconds/days of a particle of a certain size at a certain altitude
def get_lifespan(params,outputs):
    
    # Source: https://www.sciencedirect.com/science/article/pii/S0074614208602107
    # Good explainer of what's at work here: https://www.climate-policy-watcher.org/hydrology/residence-times-of-particles-in-the-troposphere.html

    if params.altitude == "Lower Troposphere":
        tau_wet = params.tau_wet_lower_troposphere
    elif params.altitude == "Upper Troposphere":
        tau_wet = params.tau_wet_upper_troposphere
    
    # Implement Jaenicke's equation
    radius_ratio = params.particle_diameter*0.5/params.residence_time_R
    tau = ( (1/params.residence_time_K)*(radius_ratio**2) + (1/params.residence_time_K)*(radius_ratio**-2) + (1/tau_wet))**-1
    
    outputs["Particle Lifespan"] = tau
    # Return output
    return tau

# Calculate GHG removed per ton particle over its lifetime
def get_lifetime_moles_GHG_per_ton(params,outputs):
    lifespan = get_lifespan(params,outputs)
    lifetime_moles_per_ton = lifespan*get_ghg_removal_rate_per_ton(params,outputs)
    outputs["Lifetime Moles Removed Per Ton"] = lifetime_moles_per_ton
    return lifetime_moles_per_ton

# Calculate cost (and GHG) to launch a ton of particles into the atmosphere at some height
def get_launch_cost_per_ton(params,outputs):
    if params.altitude == "Lower Troposphere" or params.altitude=="Upper Troposphere":
        # Use cropdusting cost estimates, dividing a plane's hourly dispersal capacity by its hourly flight cost
        cost_per_ton = params.cropduster_cost_per_ton_launched
    elif params.altitude == "Tropopause":
        # Use sulfate geoengineering cost estimates ( https://iopscience.iop.org/article/10.1088/1748-9326/aba7e7?s=09 )
        cost_per_ton = 2300 # Not used in this paper
    else:
        raise Exception("Invalid altitude zone specified.")
    outputs["Launch Cost per Ton"] = cost_per_ton
    return cost_per_ton

# Calculate the material cost of the catalyst particles. Synthesis costs besides raw materials are not considered here.
def get_material_cost_per_ton(params,outputs):
    # Figure our how much of the particle mass is in the core and shell
    core_volume = max(0,(4/3)*((0.5*params.particle_diameter-params.photocatalyst_thickness)**3))
    total_volume = (4/3)*((0.5*params.particle_diameter)**3)
    shell_volume = total_volume-core_volume
    core_mass = core_volume*params.support_density
    shell_mass = shell_volume*params.photocatalyst_density
    total_mass = core_mass+shell_mass
    # Figure out the fraction of particle mass in each
    t_core_per_t_total = core_mass/total_mass
    t_shell_per_t_total = shell_mass/total_mass
    # Calculate core and photocatalyst shell costs per ton of particles
    core_cost = params.support_cost_per_ton*t_core_per_t_total
    shell_cost = params.photocatalyst_cost_per_ton*t_shell_per_t_total
    material_cost_per_ton = core_cost+shell_cost
    # Store data in the output object
    outputs["Volume Fraction Catalyst"] = shell_volume/total_volume
    outputs["Mass Fraction Catalyst"] = shell_mass/total_mass
    outputs["Catalyst Material Cost"] = shell_cost
    outputs["Support Material Cost"] = core_cost
    outputs["Material Cost"] = material_cost_per_ton
    return material_cost_per_ton

# Calculate the cost per mole of GHG removed via photocatalytic aerosols
def get_cost_per_mole(params):
    outputs = dict() #Wipe the params dictionary
    total_cost_per_ton = get_launch_cost_per_ton(params,outputs)+get_material_cost_per_ton(params,outputs)
    lifetime_moles_per_ton = get_lifetime_moles_GHG_per_ton(params,outputs)
    cost_per_mole = total_cost_per_ton/lifetime_moles_per_ton
    outputs["Removal Cost per Mole"] = cost_per_mole
    return outputs
