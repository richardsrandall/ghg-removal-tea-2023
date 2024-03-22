

# Script for modeling cost / performance of active AMR
# System architecture based on DAC

from math import log, sqrt
from numpy import real
import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.optimize import brute
from scipy.special import lambertw
import scipy.integrate
import model.shared_assumptions as sa



# Use class "active AMR params" to provide inputs to model.
class active_model_params(sa.shared_assumptions):
    def __init__(self,which_GHG):
        super().__init__(which_GHG)

        # All shared assumptions are stored in the "shared assumptions" module.

        # Independent variables that must be explicitly set for the active contacting model
        self.AQY = None # Apparent quantum yield; dimensionless
        self.start_ppm = None # GHG concentration in ppm
        self.LCOE = None # Electricity cost in $/kWh
        self.end_ppm = None # Optional; leaving this set to "None" makes the model choose the least-cost GHG removal fraction

        # Assumptions unique to the ground-based active system
        
        # Values for evaluating transport correlations
        self.turbulent_cutoff = 2300
        self.geometry_factor = 57 # Friction factor = 57/Re for laminar flow in square channel
        self.get_ssa_from_d_h = lambda s: 4/s # Specific surface area as function of channel hydraulic diameter

        # Bounds and initial values for optimization
        self.velocity_limits = (0.1,20) # Flow velocity limits in m/s
        self.d_h_limits = (20,1000) # Hydraulic diameter limits in mm
        self.ppm_bounds = (0.01,0.99) #As fractions of the start ppm
        self.velocity_0 = 0.1 #m/s
        self.d_h_0 = 10 #mm
        self.ppm_0 = 0.1 # As fraction of the start ppm

        # Values for TEA
        self.removal_capacity = 1.0 #moles per second removal capacity of the contactor section analyzed; is normalized out and doesn't affect the outcome.
        self.LCOE = 0.02
        self.fan_and_shell_cost_per_m2 = 5000
        self.packing_cost_per_m3 = 340
        self.packing_lifespan_hours = 50000
        self.fan_efficiency = 0.56
        self.utilization_factor = 0.85 # Fraction uptime for the whole machine; Keith uses 85%
        self.non_energy_opex_capex_ratio = 0.05 # 5% capex = annualized opex
        self.contingency = 0.20 #20% contingency for an early plant

        # Functions to calculate the surface reaction rate & pressure drop
        # Surface activity should be f(velocity, hydraulic diameter, ppm of GHG, params)
        # Activity may be limited either by mass transfer or by the catalyst rate law, which is assumed to be first-order in CH4 concentration.
        self.get_surface_activity = lambda conc, vel, d, p: min(get_mass_transfer_flux_limit_to_surface(conc,vel,d,p),get_reaction_rate_limit(conc,p))
        # Pressure drop should be f(velocity, hydraulic diameter, depth, params)
        self.get_pressure_drop = monolith_pressure_drop
        
        # Values for LED
        self.LED_cost_per_kwh = 0.074 #Divide cost per kW by lifetime in hours
        self.LED_efficiency = 0.60
        self.uv_photon_energy_ev = 3.5 #Photon energy in eV

# Method that calls packing_objective_function many times to minimize contacting cost
def calculate_active_amr_cost(params):

    minimize_me = lambda x: contactor_objective_function(x[0],x[1],x[2],params,False) # Velocity, channel size, and removal fraction
    
    if params.end_ppm is None: #Iterate velocity, channel size, and removal fraction
        end_ppm_limits = (q*params.start_ppm for q in params.ppm_bounds)
        end_ppm_0 = params.ppm_0*params.start_ppm
    else: # Only iterate velocity and channel size for a user-specified removal fraction
        end_ppm_limits = (params.end_ppm,params.end_ppm)
        end_ppm_0 = params.end_ppm
    
    # Optimize using the scipy brute force algorithm with nelder-mead as a "polishing" algorithm
    dhl = tuple(params.d_h_limits)
    vl = tuple(params.velocity_limits)
    ppml = tuple(end_ppm_limits)
    rranges = (vl,dhl,ppml)
    
    # This is the number of grid points the brute-force solver tests for each variable before starting gradient descent
    n = 50 #Reduce this value to 10 or so to speed the solver up a lot; use at least 50 for final results.
    #I've previously tested the optimizer stability by randomizing the grid points the brute force optimizer tests
    
    # Use this to avoid the fact that brute force doesn't support bounds on the finishing function
    def helper(x):
        if x[0]<vl[0] or x[0]>vl[1] or x[1]<dhl[0] or x[1]>dhl[1] or x[2]<ppml[0] or x[2]>ppml[1]:
            return np.inf
        return minimize_me(x)
    
    res = brute(helper, ranges=rranges, Ns=n, full_output=True)
    r = res[0]
    #cost = contactor_objective_function(r[0],r[1],r[2],params,False)
    out = contactor_objective_function(r[0],r[1],r[2],params,True)
    
    # Compute LED costs, which are decoupled from the air contacting cost
    out["LED Opex per Mole"] = calculate_LED_cost_per_mole(params)

    # Combine into a total cost
    out["Total Cost per Mole"] = out["LED Opex per Mole"]+out["Contacting Cost per Mole"]
    
    print_everything = False #Use this to print diagnostics to see what the optimal contactor looks like
    if print_everything:
        print("Contactor data for "+params.ghg_name+" at "+str(params.start_ppm)+" ppm w/ initial rate "+str(params.reaction_rate_at_ambient_conc))
        for k in out.keys():
            print(str(k)+": "+str(out[k]))
        print("")

    # The commented part below is to help with convergence and stability testing
    #import pprint
    #pprint.pprint(out)
    #print("")
    
    # The commented part below tells us what subcomponents contribute the most contacting cost
    #print(out["Shell Cost per Mole"]/out["Total Cost per Mole"])
    #print(out["Fan Opex per Mole"]/out["Total Cost per Mole"])
    #print("")

    # Return the result    
    return out


# Objective function for optimization 

def contactor_objective_function(velocity, channel_hydraulic_diameter, end_ppm, params, return_all_data = False):

    # Find dimensionless numbers
    d_h = channel_hydraulic_diameter*params.m_per_mm # For normalization, value is fed in mm
    reynolds = velocity*d_h/params.air_kin_visc
    # Find moles GHG removed per m3 of air throughput; find frontal area for desired removal capacity
    moles_GHG_removed_per_m3_air = (params.start_ppm - end_ppm)*params.moles_per_m3_per_ppm
    frontal_area = (params.removal_capacity/moles_GHG_removed_per_m3_air)/velocity

    # Calculate depth of packing necessary to hit certain GHG removal fraction
    # Volume-specific surface area (m2/m3) is a useful intermediate
    specific_surface_area = params.get_ssa_from_d_h(d_h)

    #Integrate to find the packing depth required to achieve the specified outlet GHG concentration

    # Volume flow rate is frontal_area*velocity
    # Total moles of GHG flowing in/out of a control volume per second equals frontal_area*velocity*moles_per_m3_per_ppm*C_ghg
    # To remove d_ppm from the stream, you need to remove d_moles = frontal_area*velocity*moles_per_m3_per_ppm*d_ppm, requiring
    # an extra surface area of dA = d_moles/rate = (frontal_area*velocity*moles_per_m3_per_ppm*d_ppm)/rate
    # and an extra packing depth of d_Depth = dA/(frontal_area*SSA) =  (frontal_area*velocity*moles_per_m3_per_ppm*d_ppm)/(rate*frontal_area*SSA)

    d_depth_over_d_ppm = lambda ppm_ghg: (frontal_area*params.moles_per_m3_per_ppm*velocity)/(params.get_surface_activity(ppm_ghg,velocity,d_h,params)*frontal_area*specific_surface_area)
    packing_depth = scipy.integrate.quad(d_depth_over_d_ppm, end_ppm, params.start_ppm)[0]
    
    # Calculate pressure drop for given depth, assuming a square channel profile
    delta_P = params.get_pressure_drop(velocity, d_h, packing_depth, params)
    
    # Calculate frontal area and volume needed
    packing_volume = frontal_area*packing_depth
    
    # Calculate the contacting costs for whatever removal capacity (e.g. 1 mole/second) was specified; return value is a tuple with different cost components
    contacting_costs_total = calculate_annualized_contacting_costs(delta_P,velocity,packing_depth,frontal_area,moles_GHG_removed_per_m3_air,specific_surface_area,params)
    # Normalize by whatever removal capacity (e.g. 1 mole/second) was specified at the beginning
    contacting_costs = dict()
    for k in contacting_costs_total.keys():
        contacting_costs[k] = contacting_costs_total[k]/(params.removal_capacity*params.seconds_per_year)
    
    if return_all_data:
        out = contacting_costs
        out["Depth"] = packing_depth
        out["Velocity"] = velocity
        out["Channel Size"] = channel_hydraulic_diameter
        out["Start ppm"] = params.start_ppm
        out["LCOE"] = params.LCOE
        out["End ppm"] =end_ppm
        out["Delta-P"] = delta_P
        out["Reynolds"] = reynolds
        out["Frontal Area per Mole per Second"] = frontal_area/params.removal_capacity
        return out
    else:
        return sum(contacting_costs.values())

# Do the TEA and compute capex/opex for contactor and fans according to Keith and Holmes methodology
def calculate_annualized_contacting_costs(delta_P,flow_speed,packing_depth,frontal_area,moles_GHG_removed_per_m3_air,specific_surface_area,params):

    # All of these values are defined for whatever removal capacity (e.g. 1 mole/second) was input at the beginning

    # Get packing cost and shell cost for the given frontal area and depth
    catalyst_cost_per_m2 = (1e-6*params.photocatalyst_thickness)*(params.photocatalyst_density)*(params.photocatalyst_cost_per_ton*0.001)
    packing_cost = frontal_area*packing_depth*(params.packing_cost_per_m3+catalyst_cost_per_m2*specific_surface_area)
    fan_and_shell_cost = params.fan_and_shell_cost_per_m2*frontal_area

    # Account for CRF and utilization factor to get annual costs for each
    yearly_packing_capex = packing_cost*params.CRF/params.utilization_factor
    yearly_packing_opex = packing_cost*1/(params.packing_lifespan_hours*params.hours_per_year) #Replace the packing every n years, i.e. replace 1/n of the packing every year as opex
    yearly_shell_capex = fan_and_shell_cost*params.CRF/params.utilization_factor
    # Account for the cost of replacing the packing
    yearly_non_energy_opex = params.non_energy_opex_capex_ratio*(yearly_packing_capex+yearly_shell_capex)+yearly_packing_opex
    yearly_contingency = (packing_cost+fan_and_shell_cost)*params.contingency*params.CRF
    
    # Fan energy
    fan_joules_per_mole_removed = delta_P*(1/moles_GHG_removed_per_m3_air)*(1/params.fan_efficiency)# Make sure i'm doing this right!!!
    fan_watts = params.removal_capacity*fan_joules_per_mole_removed
    yearly_fan_opex = params.seconds_per_year*params.LCOE*fan_watts/params.j_per_kwh

    out = dict()
    out["Packing Cost per Mole"] = yearly_packing_capex
    out["Shell Cost per Mole"] = yearly_shell_capex
    out["Non-Energy Opex per Mole"] = yearly_non_energy_opex
    out["Fan Opex per Mole"] = yearly_fan_opex
    out["Contingency per Mole"] = yearly_contingency
    out["Capex per Mole"] = yearly_packing_capex+yearly_shell_capex+yearly_contingency
    out["Contacting Cost per Mole"] = yearly_packing_capex+yearly_shell_capex+yearly_contingency+yearly_fan_opex+yearly_non_energy_opex
    # Return results

    return out

# Figure out how much we're spending on LEDs -- calculated separately because it doesn't affect the optimization (and so doesn't belong in the objective function)
def calculate_LED_cost_per_mole(params):
    # If AQY has not been set, assume the user wants to ignore LEDs for now
    if params.AQY is None:
        return 0
    # Find LED cost from AQY and surface activity
    kwh_per_mole_GHG = params.avogadro*params.uv_photon_energy_ev/(params.AQY*params.ev_per_j*params.j_per_kwh*params.LED_efficiency)
    LED_cost_per_mole_GHG = (params.LCOE+params.LED_cost_per_kwh)*kwh_per_mole_GHG #LED cost is measured per watt of input, not output.
    return LED_cost_per_mole_GHG

# Calculate pressure drop in a square channel
# We sanity-checked this script using values from an automotive catalytic converter & also values from Keith et al. with good results
def monolith_pressure_drop(velocity, d_h, depth, params):
    c_f = params.geometry_factor
    reynolds = velocity*d_h/params.air_kin_visc
    # Plug into one of several relations from literature.
    if reynolds<2000: # An exact solution exists for laminar flow
        friction_factor = c_f/reynolds
    elif reynolds<3000:
        return 1e10
    else: # Explicit formula for friction factor
        friction_factor = (0.79*log(reynolds)-1.64)**-2
    delta_P = depth*friction_factor*(0.5*params.rho_air)*(velocity**2)/d_h
    return delta_P

# Calculate how fast GHGs can transport to the surface of the monolith
def get_mass_transfer_flux_limit_to_surface(concentration, velocity, d_h, params):

    # Calculate Reyonolds no.; Schmidt no. is already in the parameters object
    reynolds_number = velocity*d_h/params.air_kin_visc
    # Calculate mass transfer limit from Dittus-Boelter eqn or closed-form solution
    if reynolds_number>3000:#params.turbulent_cutoff:
        f = (0.79*log(reynolds_number)-1.64)**-2
        sh = ((f/8)*(reynolds_number-1000)*params.schmidt_number)/(1+12.7*((f/8)**0.5)*((params.schmidt_number**(2/3))-1))
        #Gnielinski equation
    elif reynolds_number>2000:
        return 1e-10
    else:
        sh = 3.66 #Sherwood no. is constant for laminar flow in a square pipe
    h_m = sh*params.ghg_diffusivity/d_h
    surface_ppm = 0 #Assumption for sake of model
    mass_transfer_flux_limit_moles = h_m*params.moles_per_m3_per_ppm*(concentration-surface_ppm)

    return mass_transfer_flux_limit_moles

# Calculate the reaction rate at the surface of the monolith in the absence of mass transfer limitations
def get_reaction_rate_limit(conc,params):
    rate_at_ambient_conc = params.reaction_rate_at_ambient_conc #Should be input in moles/m2s
    if params.ghg_name == 'CH4':
        ambient_conc = 1.8
    else:
        ambient_conc = 0.33
    return rate_at_ambient_conc*(conc/ambient_conc)