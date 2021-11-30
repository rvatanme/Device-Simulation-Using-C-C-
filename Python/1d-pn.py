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

# Setup the normalized doping by ni for a pn diode with eqaul N-side and P-side length size        
dop = [0]*n_max
for i in range(n_max):
    if i <= n_max/2:
        dop[i] = -Na/ni         # P-side from 0 to x_max/2
    else:
        dop[i] = Nd/ni          # N-side from x_max/2 to x_max
        
 # Initialize the normalized potential based on the requirement of charge neutrality throughout the whole structure
fi = [0]*n_max
n = [0]*n_max
p = [0]*n_max
for i in range(n_max):
    zz = 0.5*dop[i]
    if zz > 0:
        xx = zz*(1 + math.sqrt(1+1/(zz*zz)))
    if zz <  0:
        xx = zz*(1 - math.sqrt(1+1/(zz*zz)))       # Taylor series of sqrt(1+x) = 1 + x/2
    fi[i] = math.log(xx)                           # \phi = (kT/q)ln(Nd/ni)  or  \phi = (kT/q)ln(Na/ni) in kT/q unit
    n[i] = xx                                      # n = Nd in N-side  and  n = ni*2/Na in P-side in ni unit
    p[i] = 1/xx                                    # p = Na in P-side  and  n = ni*2/Nd in P-side in ni unit

delta_acc = 1E-5               # Set the Tolerance

##########################################################################
##                                                                      ##
##               EQUILIBRIUM  SOLUTION PART BEGINS                      ##
##                                                                      ##
##########################################################################

#(A) Define the elements of the coefficient matrix for the internal nodes and
#    initialize the forcing function
a = [0]*n_max
c = [0]*n_max
b = [0]*n_max
f = [0]*n_max
dx2 = dx*dx
for i in range(n_max):      # Setup the coefficients of discritized Poisson equation 
    a[i] = 1/dx2
    c[i] = 1/dx2
    b[i] = -(2/dx2+math.exp(fi[i])+math.exp(-fi[i]))
    f[i] = math.exp(fi[i]) - math.exp(-fi[i]) - dop[i] - fi[i]*(math.exp(fi[i])+math.exp(-fi[i]))

#(B) Define the elements of the coefficient matrix and initialize the forcing function at the ohmic contacts
a[0] = 0                      # Setup the Dirichlet boundary condition for the boundary nodes
c[0] = 0
b[0] = 1
f[0] = fi[0]
a[n_max-1] = 0
c[n_max-1] = 0
b[n_max-1] = 1
f[n_max-1] = fi[n_max-1]

#  Start the iterative procedure for the solution of the linearized Poisson equation using LU decomposition method:
flag_conv = 0           # convergence of the Poisson loop
k_iter= 0
alpha = [0]*n_max
beta = [0]*n_max
v = [0]*n_max
delta = [0]*n_max
while flag_conv == 0:          
    k_iter = k_iter + 1                      #Count iteration cycles
    alpha[0] = b[0]
    for i in range(1, n_max):                #Set up L and U matrix 
        beta[i]=a[i]/alpha[i-1]
        alpha[i]=b[i]-beta[i]*c[i-1]
    
    # Solution of Lv = f %    
    v[0] = f[0]
    for i in range(1,n_max):                      # Calculate v array using L matrix and fi array from the previous cycle
        v[i] = f[i] - beta[i]*v[i-1]
    
    # Solution of U*fi = v %
    temp = v[n_max - 1]/alpha[n_max - 1]          # Calculate the new fi using U matrix and v array.
    delta[n_max - 1] = temp - fi[n_max - 1]
    fi[n_max - 1]=temp
    for i in range(n_max-1):       #delta
        i = n_max - 2 - i
        temp = (v[i]-c[i]*fi[i+1])/alpha[i]
        delta[i] = temp - fi[i]
        fi[i] = temp
    
    delta_max = 0
    
    for i in range(n_max):                     # Finding the maximum element difference between new and old fi 
        xx = abs(delta[i])
        if xx > delta_max:
            delta_max=xx
            
    print (delta_max)
    if(delta_max < delta_acc):                # Check if the maximum difference is less or greater than threshold 
        flag_conv = 1
    else:
        for i in range(1,n_max-1):            # If difference is greater then the new b and f would be calculated and used in the new iter
            b[i] = -(2/dx2 + math.exp(fi[i]) + math.exp(-fi[i]));
            f[i] = math.exp(fi[i]) - math.exp(-fi[i]) - dop[i] - fi[i]*(math.exp(fi[i]) + math.exp(-fi[i]))

# End of Poisson-Boltzmann solver
            
# Printing out the results
xx1 = [0]*n_max
Ec = [0]*n_max
ro = [0]*n_max
el_field1 = [0]*n_max
el_field2 = [0]*n_max
xx1[i] = dx*1e4
for i in range(1,n_max-1):       # Calculate electric field, charge carrier densities and conduction band profile based on calculated fi
    Ec[i] = dEc - Vt*fi[i];    
    ro[i] = -ni*(math.exp(fi[i]) - math.exp(-fi[i]) - dop[i])
    el_field1[i] = -(fi[i+1] - fi[i])*Vt/(dx*Ldi)
    el_field2[i] = -(fi[i+1] - fi[i-1])*Vt/(2*dx*Ldi)
    n[i] = math.exp(fi[i])
    p[i] = math.exp(-fi[i])
    xx1[i] = xx1[i-1] + dx*Ldi*1e4



Ec[0] = Ec[1]
Ec[n_max-1] = Ec[n_max-2]
xx1[n_max-1] = xx1[n_max-2] + dx*Ldi*1e4
el_field1[0] = el_field1[1]
el_field2[0] = el_field2[1]
el_field1[n_max-1] = el_field1[n_max-2]
el_field2[n_max-1] = el_field2[n_max-2]
nf = [x*ni for x in n]
pf = [x*ni for x in p]
ro[1] = ro[2]
ro[n_max-1] = ro[n_max-2]


f1 = open('els_plot_py.txt', 'w')
f1.write("x [um]\tPotential [V]\tElectric Field1 [V/cm]\tElectric Field2 [V/cm]\tElectron Densities [1/cm^3]\tHole Densities [1/cm^3]\tTotal Charge Density [C/cm^3]\tConduction Band Energy (eV)\n")
for i in range(n_max):
    f1.write(str(xx1[i])+'\t'+str(Vt*fi[i])+'\t'+str(el_field1[i])+'\t'+str(el_field2[i])+'\t'+str(nf[i])+'\t'+str(pf[i])+'\t'+str(ro[i])+'\t'+str(Ec[i])+'\n')
f1.close()
