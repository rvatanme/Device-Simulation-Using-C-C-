import numpy as np
import math

def abs(x):
    if x >= 0:
        return x
    else:
        return -x

    
def BER(x):
    if x > 1.0E-10:
        return x*math.exp(-x)/(1-math.exp(-x))
    elif x < 0.0 and abs(x) > 1.0E-10:
        return x/(math.exp(x)-1)
    elif (x == 0.0):
        return 1.0
    else:    
        temp_term = 1
        sum = temp_term
        i = 0.0
        while flag_sum == 0:
            i += 1
            temp_term = temp_term*x/(i+1)
            if abs(temp_term) < 1.0E-14:
                flag_sum = 1
            else:
                sum += temp_term
        return 1/sum

# Defining the Fundamental and Material Constants 
q = 1.602E-19;        # C or [J/eV]
kb = 1.38E-23;        # [J/K]
eps = 1.05E-12;       # This includes the eps  = 11.7 for Si [F/cm]
T = 300.0;            # [K]
ni = 1.5E10;          # Intrinsic carrier concentration [1/cm^3]
Vt = kb*T/q;          # [eV]
RNc = 2.8E19;         # This is 2.8e20 in the FORTRAN file
TAUN0 = 0.1E-6;       # Electron SRH life time
TAUP0 = 0.1E-6;       # Hole SRH life time
mun0   = 1500.0;      # Electron Mobility in cm2/V-s
mup0   = 1000.0;      # Hole Mobility in cm2/V-s
dEc = Vt*math.log(RNc/ni);

# Define Doping Values
Na = 1.0E16;                                      # [1/cm^3]
Nd = 1.0E17;                                      # [1/cm^3]

# Calculate relevant parameters for the simulation 
Vbi = Vt*math.log(Na*Nd/(ni*ni));                 # Built-in Potential Energy [eV] at zero bias (open-circuit)
W   = math.sqrt(2*eps*(Na+Nd)*Vbi/(q*Na*Nd));     # Depletion width [cm]
Wn  = W*math.sqrt(Na/(Na+Nd));                    # N-side depletion width [cm]
Wp  = W*math.sqrt(Nd/(Na+Nd));                    # P-side depletion width [cm]
E_p = q*Nd*Wn/eps;                                # Maximum electric field [V/cm] at the junction
Ldn = math.sqrt(eps*Vt/(q*Nd));                   # Debye length for N-side
Ldp = math.sqrt(eps*Vt/(q*Na));                   # Debye length for P-side
Ldi = math.sqrt(eps*Vt/(q*ni));                   # Debye length for Intrinsic Si

# Setting the size of the simulated pn diode (x_max) based on the analytical results from the width of the 
# depletion regions for a simple pn-diode
x_max = 0.0
if(x_max < Wn):
    x_max = Wn
if(x_max < Wp):
    x_max = Wp
x_max = 20*x_max

# Setting the uniform grid size (dx) based on the extrinsic Debye lengths
dx = Ldn
if(dx > Ldp):
    dx=Ldp
dx = dx/20

# Calculate the required number of grid points and normalize dx by Intrinsic Debye length (Ldi)
n_max = math.floor(x_max/dx) +1
dx = dx/Ldi                # Renormalize lengths with Ldi

# Setup the normalized doping for a pn diode with eqaul N-side and P-side length size
dop = [0.0]*n_max
for i in range(n_max):
    if i <= n_max/2:
        dop[i] = -Na/ni
    else:
        dop[i] = Nd/ni

               




               


