
# Script for modeling cost and performance of passive, solar, ground-based CH4/N2O removal
import model.shared_assumptions as sa
from math import sin, cos, radians

# Object to hold assumptions and inputs for single runs of the script
class passive_model_params(sa.shared_assumptions):
    def __init__(self,which_GHG):
        super().__init__(which_GHG)

        # All shared assumptions are stored in the "shared assumptions" module.

        # Independent variables that must be explicitly set
        self.AQY = None # Apparent quantum yield; dimensionless
        self.air_velocity = None # Free stream velocity in the atmosphere; m/s
        self.ghg_ppm = None # GHG concentration in ppm
        self.which_GHG = which_GHG

        # Assumptions and values unique to the ground-based passive model
        self.capacity_factor = 0.25 #Change to availability or utilization factor
        self.painting_capex_per_m2 = 10.8 #$/m2
        self.opex_per_m2 = 1.6 #$/m2-y
        self.catalyst_capex_per_m2 = 0.86 #$/m2

        # The correlation that will be used to find the sherwood number and the corresponding characteristic length
        self.sherwood_correlation = lambda reynolds, schmidt: 1.1*0.037*(reynolds**(4.0/5)-0)*(schmidt**(1.0/3))
        # Assume roughness factor of 1.1 for 'smooth plaster', turbulent flow as LBL paper found generally applied
        # Source: https://gundog.lbl.gov/dirpubs/47275.pdf
        
        # Set the characteristic length according to the effective length formula given
        self.l_char = 44.1
        
        
        # If you want to run the 'solar PV-like' model, uncomment the following Sherwood correlation:
        #self.sherwood_correlation = lambda reynolds, schmidt: 0.905*(reynolds**(1/2))*(schmidt**(1/3)) # For angled square plates
        # And the following l_char, costs
        #self.painting_capex_per_m2 = 76
        #self.opex_per_m2 = 2
        #self.l_char = 1.5 #Length of panel

def calculate_surface_activity(params,outputs): # Returns the average GHD drawdown flux in moles per m2s
    
    # From AQY and solar flux, calculate max reaction flux
    AQY_limit_moles = params.AQY*params.solar_uv_photon_flux/params.avogadro
    # From Sherwood number, calculate max diffusion flux
    # Retrieve Schmidt number
    schmidt = params.schmidt_number
    # Find Reynolds number
    reynolds = params.air_velocity*params.l_char/params.air_kin_visc
    # Plug into a convective heat/mass transfer relationship
    sherwood = params.sherwood_correlation(reynolds, schmidt)
    h_m = sherwood*params.ghg_diffusivity/params.l_char
    # Calculate flux based on concentration gradient
    surface_ppm = 0
    bulk_ppm = params.ghg_ppm
    mass_transfer_limit_moles = h_m*params.moles_per_m3_per_ppm*(bulk_ppm-surface_ppm)
    # Determine rate-limiting step
    activity = min([mass_transfer_limit_moles, AQY_limit_moles])
    # Adjust for fraction of time that is daylight
    adjusted_activity = params.capacity_factor*activity
    # Log helpful info to the output dict and return the desired value
    outputs["Mass Transfer Flux Limit"] = mass_transfer_limit_moles
    outputs["Reaction Flux Limit"] = AQY_limit_moles
    return adjusted_activity

def calculate_cost_per_m2_year(params,outputs): # Returns the annualized cost per m2 of mounted active phase

    # Get annualized capex
    annual_capex = params.CRF*(params.painting_capex_per_m2+params.catalyst_capex_per_m2)
    # Get opex
    annual_opex = params.opex_per_m2
    # Combine them
    annual_cost = annual_capex + annual_opex
    outputs["Capex Fraction"] = annual_capex/annual_cost
    outputs["Opex Fraction"] = annual_opex/annual_cost
    return annual_cost
    
def calculate_figures_of_merit(params):

    outputs = dict()

    # Calculate total GHG drawdown in moles/m2-year
    activity = calculate_surface_activity(params,outputs) #in moles/m2-sec

    # Get all-in cost per m2-year
    cost_per_m2_year = calculate_cost_per_m2_year(params,outputs)

    # Compute and organize the different output values
    outputs["Cost per Mole"] = cost_per_m2_year/(float(activity)*params.seconds_per_year)
    outputs["Surface Activity"] = activity
    
    # Return them
    return outputs


