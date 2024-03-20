
# This module contains functions to shuffle between {$/mole, $/ton} of {CH4, N2O, CO2e}

# Fixed values
ch4_molar_mass = 16.04
n2o_molar_mass = 44.01
co2_molar_mass = 44.01
ch4_gwp = 81.2 #GWP20 -- an earlier version used GWP100
n2o_gwp = 273 #GWP20

# Convert something (probably $) per mole to per ton GHG
def per_mole_to_per_ton(cost, which_GHG):
    if which_GHG == 'CH4':
        out = cost*(1e6/ch4_molar_mass) # Don't adjust for GWP at this point
    elif which_GHG == 'N2O':
        out = cost*(1e6/n2o_molar_mass)
    elif which_GHG == 'CO2':
        out = cost*(1e6/co2_molar_mass)
    else:
        out = 0
    return out

# Convert something (probably $) per mole to per ton CO2e
def per_mole_to_per_ton_CO2e(cost, which_GHG):
    if which_GHG == 'CH4':
        out = cost*(1e6/ch4_molar_mass)/ch4_gwp # Adjust for global warming potentials
    elif which_GHG == 'N2O':
        out = cost*(1e6/n2o_molar_mass)/n2o_gwp
    elif which_GHG == 'CO2':
        out = cost*(1e6/co2_molar_mass)
    else:
        out = 0
    return out

# Convert from $/tCO2e from removed CH4 to $/tCO2e for removing the same # of moles of N2O
def secax_helper_forward(value):
    return value*(ch4_gwp/n2o_gwp)*(ch4_molar_mass/n2o_molar_mass)

# Convert from $/tCO2e from removed N2O to $/tCO2e for removing the same # of moles of CH4
def secax_helper_backward(value):
    return value*(n2o_gwp/ch4_gwp)*(n2o_molar_mass/ch4_molar_mass)

# Return a tuple of the previous two functions
def get_secax_helpers():
    return (secax_helper_forward, secax_helper_backward)
