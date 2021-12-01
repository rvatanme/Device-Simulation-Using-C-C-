import numpy as np
import math

def abs(x):
    if x >= 0:
        return x
    else:
        return -x

    
def BER(x):
    flag_sum = 0
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
    for i in range(n_max-2,-1,-1):       #delta
        temp = (v[i]-c[i]*fi[i+1])/alpha[i]
        delta[i] = temp - fi[i]
        fi[i] = temp
    
    delta_max = 0
    
    for i in range(n_max):                     # Finding the maximum element difference between new and old fi 
        xx = abs(delta[i])
        if xx > delta_max:
            delta_max=xx
            
    #print (delta_max)
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

##########################################################################
##                 END OF EQUILIBRIUM  SOLUTION PART                    ##
##########################################################################

##########################################################################
##                                                                      ##
##               NON-EQUILIBRIUM  SOLUTION PART BEGINS                  ##
##                                                                      ##
##########################################################################

##########################################################################
##              1. Calculate Low filed mobility                         ## 
##########################################################################

##  Prameters for Low field mobility calculation 

TL = 300                    # Temp in Kelvin
N  = Na + Nd                # Local (total) impurity concentration

MU1N_CAUG   = 55.24         # cm2/(V.s)
MU2N_CAUG   = 1429.23       # cm2/(V.s)
ALPHAN_CAUG = 0.0           # unitless
BETAN_CAUG  = -2.3          # unitless
GAMMAN_CAUG = -3.8          # unitless
DELTAN_CAUG = 0.73          # unitless
NCRITN_CAUG = 1.072E17      # cm-3

MU1P_CAUG   = 49.7          # cm2/(V.s)
MU2P_CAUG   = 479.37        # cm2/(V.s)
ALPHAP_CAUG = 0.0           # unitless
BETAP_CAUG  = -2.2          # unitless
GAMMAP_CAUG = 13.7          # unitless
DELTAP_CAUG = 0.70          # unitless
NCRITP_CAUG = 1.606E17      # cm-3
BETAN = 2.0
BETAP = 1.0

VSATN = (2.4E7) / (1 + 0.8*math.exp(TL/600))  # Saturation Velocity of Electrons
VSATP = VSATN                                 # Saturation Velocity of Holes

#################### END of Low Field Mobility Calculation ###############


##########################################################################
##   2. Start the main Loop to increment the Anode voltage by Vt=KbT/q  ## 
##      till it reaches 0.625V.                                         ##
##########################################################################
vindex=0
Vplot = [0]*74
Ef = [0]*n_max
Ev = [0]*n_max
Ei = [0]*n_max
Efn = [0]*n_max
Efp = [0]*n_max
mup = [0]*n_max
mun = [0]*n_max
an = [0]*n_max
bn = [0]*n_max
cn = [0]*n_max
fn = [0]*n_max
ap = [0]*n_max
bp = [0]*n_max
cp = [0]*n_max
fp = [0]*n_max
alphan = [0]*n_max
betan = [0]*n_max
vn = [0]*n_max
alphap = [0]*n_max
betap = [0]*n_max
vp = [0]*n_max
Jnim1by2 = [[0]*n_max]*74
Jnip1by2 = [[0]*n_max]*74
Jelec = [[0]*n_max]*74
Jpim1by2 = [[0]*n_max]*74
Jpip1by2 = [[0]*n_max]*74
Jhole = [[0]*n_max]*74
Jtotal = [[0]*n_max]*74

for iV in range(74):                # Start VA increment loop

    VA = iV*0.33*Vt
    Each_Step   = 0.33*Vt 
    Total_Steps = math.ceil(0.625/(0.33*Vt))
    Vplot[iV] = VA
    
    fi[0] = fi[0] + VA            # Apply potential to Anode (1st node)
    
    flag_conv2 = 0               # Convergence of the Poisson loop
    k_itern= 0
    
    ## Initialize the First and Last Node for Poisson's eqn
    a[0] = 0
    c[0] = 0
    b[0] = 1
    f[0] = fi[0]
    a[n_max-1] = 0
    c[n_max-1] = 0
    b[n_max-1] = 1
    f[n_max-1] = fi[n_max-1]
    
    #######################################################################
    ## 3. Start the Poisson equation solver loop to calculate the        ##
    ##    potential for each Anode voltage increase                      ##
    #######################################################################
        
    #######################################################################
    ## 3. Start the Poisson equation solver loop to calculate the        ##
    ##    potential for each Anode voltage increase                      ##
    #######################################################################
    
    while flag_conv2 == 0:             # Start Poisson's eqn
        
        k_itern = k_itern + 1

        
        #######################################################################
        ## 3.1 . Calculate Field Dependant Mobility for each value of 'fi'   ##
        ##       at each node point of the PN diode.                         ##
        #######################################################################
                
### To test with Constant Mobility without field dependancy.
        ## Calculate the Electric Field at each Node
        for i in range(1,n_max-1):
            Ef[i] = abs(fi[i] - fi[i+1])*Vt/(dx*Ldi)
        
        Ef[0] = Ef[1]
        Ef[n_max-2] = Ef[n_max-1]
        
        for i in range(n_max):           
            pdeno  = pow((mup0 * Ef[i] / VSATP) , BETAP)
            mup[i]   = mup0 * (pow( (1/(1 + pdeno)) , (1/BETAP)))
                            
            ndeno  = pow((mun0 * Ef[i] / VSATN) , BETAN)
            mun[i] = mun0 * (pow( (1/(1 + ndeno)) , (1/BETAN)))
            
        mup[0]     = mup[1];
        mup[n_max-1] = mup[n_max-2];
        
        mun[0]     = mun[1];
        mun[n_max-1] = mun[n_max-2];
        ############# END of FIELD Dependant Mobility Calculation ############
        
                      
        #############################################################################
        ## Solve Continuity Equation for Electron and Holes using LU Decomposition ##                                    
        #############################################################################
        
        #(A) Define the elements of the coefficient matrix and initialize the forcing
        #    function at the ohmic contacts for ELECTRON and HOLE Continuity Eqns 


        an[0] = 0              #Co-ef for electron at Anode
        bn[0] = 1              #Co-ef for electron at Anode
        cn[0] = 0              #Co-ef for electron at Anode
        ap[0] = 0              #Co-ef for hole     at Anode
        bp[0] = 1              #Co-ef for hole     at Anode
        cp[0] = 0              #Co-ef for hole     at Anode
        fn[0] = n[0]
        fp[0] = p[0]
        
        an[n_max-1] = 0          #Co-ef for electron at Cathode
        bn[n_max-1] = 1          #Co-ef for electron at Cathode
        cn[n_max-1] = 0          #Co-ef for electron at Cathode
        ap[n_max-1] = 0          #Co-ef for hole     at Cathode
        bp[n_max-1] = 1          #Co-ef for hole     at Cathode
        cp[n_max-1] = 0          #Co-ef for hole     at Cathode
        fn[n_max-1] = n[n_max-1]
        fp[n_max-1] = p[n_max-1]


        #(B) Define the elements of the coefficient matrix for the internal nodes and
        #    initialize the forcing function
        for i in range(1,n_max-1):
            munim1by2 = (mun[i-1]+mun[i])/2         
            munip1by2 = (mun[i]+mun[i+1])/2         
            mupim1by2 = (mup[i-1]+mup[i])/2         
            mupip1by2 = (mup[i]+mup[i+1])/2
            
            ## Co-efficients for HOLE Continuity eqn
            cp[i] = mupip1by2 * BER(fi[i] - fi[i+1])
            ap[i] = mupim1by2 * BER(fi[i] - fi[i-1])
            bp[i] = -( mupim1by2 * BER(fi[i-1] - fi[i]) + mupip1by2 * BER(fi[i+1] - fi[i]))
            ## Co-efficients for ELECTRON Continuity eqn
            cn[i] = munip1by2 * BER(fi[i+1] - fi[i])
            an[i] = munim1by2 * BER(fi[i-1] - fi[i])
            bn[i] = -( munim1by2 * BER(fi[i] - fi[i-1]) + munip1by2 * BER(fi[i] - fi[i+1]))

            ## Forcing Function for ELECTRON and HOLE Continuity eqns
            fn[i] = (Ldi*Ldi*dx2/Vt) * ( p[i]*n[i] - 1 ) / ( TAUP0*(n[i] + 1) + TAUN0*(p[i]+1))
            fp[i] = (Ldi*Ldi*dx2/Vt) * ( p[i]*n[i] - 1 ) / ( TAUP0*(n[i] + 1) + TAUN0*(p[i]+1))
            
        #(C)  Start the iterative procedure for the solution of the linearized Continuity
        
        ###########  Start "ELECTRON" continuity solver using LU decomposition method  ###########
        alphan[0] = bn[0]
        for i in range(1,n_max):
            betan[i]=an[i]/alphan[i-1]
            alphan[i]=bn[i]-betan[i]*cn[i-1]
        
        ## Solution of Lv = f %    

        vn[0] = fn[0]
        for i in range(1,n_max):
            vn[i] = fn[i] - betan[i]*vn[i-1]

        # Solution of U*fi = v %
        tempn = vn[n_max - 1]/alphan[n_max - 1];    
        n[n_max - 1]=tempn
        for i in range(n_max-2, -1, -1):            
            tempn = (vn[i]-cn[i]*n[i+1])/alphan[i]
            n[i] = tempn       
       ################## END of ELECTRON Continuty Solver ##################
        
            
        #(D)  Start the iterative procedure for the solution of the linearized Continuity
        #########  Start of HOLE continuity solver using LU decomposition method  #######
        alphap[0] = bp[0];
        for i in range(1,n_max):
            betap[i]=ap[i]/alphap[i-1]
            alphap[i]=bp[i]-betap[i]*cp[i-1]

        # Solution of Lv = f %

        vp[0] = fp[0]
        for i in range(1,n_max):
            vp[i] = fp[i] - betap[i]*vp[i-1]

        # Solution of U*fi = v %    

        tempp = vp[n_max - 1]/alphap[n_max - 1];    
        p[n_max - 1]=tempp;
        for i in range(n_max - 2, -1, -1):
            tempp = (vp[i]-cp[i]*p[i+1])/alphap[i]  
            p[i] = tempp    
       ################## END of HOLE Continuty Solver ###############
    
       
       ################################################################# 
       ## Calculate potential fi again with new values of "n" and "p" ##
       ##     and check for convergence                               ##     
       #################################################################
       
       # Recalculate forcing function and central coefficient b for fi 

        for i in range(1,n_max-1):
            b[i] = -(2/dx2 + n[i] + p[i])
            f[i] = n[i] - p[i] - dop[i] - (fi[i]*(n[i] + p[i]))
       ## here values of n(i) and p(i) are used in place of exp(fi(i))

        
        # Solve for Updated potential given the new value of Forcing 
        # Function using LU decomposition 
        
        alpha[0] = b[0]
        for i in range(1,n_max):
            beta[i]=a[i]/alpha[i-1]
            alpha[i]=b[i]-beta[i]*c[i-1]

        # Solution of Lv = f 

        v[0] = f[0];
        for i in range(1,n_max):
            v[i] = f[i] - beta[i]*v[i-1]

        # Solution of U*fi = v
        temp = v[n_max - 1]/alpha[n_max - 1]
        delta[n_max - 1] = temp - fi[n_max - 1]
        fi[n_max - 1] = temp
        for i in range(n_max - 2, -1, -1):
            temp = (v[i]-c[i]*fi[i+1])/alpha[i]
            delta[i] = temp - fi[i]
            fi[i] = temp

        delta_max = 0
        for i in range(n_max):
            xx = abs(delta[i])
            if xx > delta_max:
                delta_max = xx        ## Calculate the max error

        # Test convergence and start the loop if necessary else increase
        # the applied potential
        #print (delta_max)
        if(delta_max < delta_acc):
            flag_conv2 = 1
    #####################################################################
    ##                        CALCULATE CURRENT                        ##
    #####################################################################

    # Electron Current       
    for i in range(1,n_max-1): 
        Jnim1by2[iV][i] = (q*mun[i]*Vt/(dx*Ldi)) * ni*( n[i]*BER(fi[i]-fi[i-1]) - n[i-1]*BER(fi[i-1]-fi[i]) )
        Jnip1by2[iV][i] = (q*mun[i]*Vt/(dx*Ldi)) * ni*( n[i+1]*BER(fi[i+1]-fi[i]) - n[i]*BER(fi[i]-fi[i+1]) )
        Jelec[iV][i] = (Jnip1by2[iV][i] + Jnim1by2[iV][i])/2
        
        Jpim1by2[iV][i] = (q*mup[i]*Vt/(dx*Ldi)) * ni*( p[i]*BER((fi[i-1]-fi[i])) - p[i-1]*BER((fi[i]-fi[i-1])) )
        Jpip1by2[iV][i] = (q*mup[i]*Vt/(dx*Ldi)) * ni*( p[i+1]*BER((fi[i]-fi[i+1])) - p[i]*BER((fi[i+1]-fi[i])) )
        Jhole[iV][i] = (Jpip1by2[iV][i] + Jpim1by2[iV][i])/2

        Jtotal[iV][i] = Jelec[iV][i] + Jhole[iV][i]
        
# End of main FOR loop for VA increment.
##########################################################################
##                 END OF NON-EQUILIBRIUM  SOLUTION PART                ##
##########################################################################


# Write the results of the simulation in files %

xx1[0] = dx*1e4
for i in range(1,n_max - 1): 
    Ec[i] = dEc - Vt*fi[i]     #Values from the second Node%
    ro[i] = -ni*(n[i] - p[i] - dop[i])
    el_field1[i] = -(fi[i+1] - fi[i])*Vt/(dx*Ldi)
    el_field2[i] = -(fi[i+1] - fi[i-1])*Vt/(2*dx*Ldi)
    xx1[i] = xx1[i-1] + dx*Ldi*1e4

for i in range(74):
    Jtotal[i][0] = Jtotal[i][1]
    Jelec[i][0] = Jelec[i][1]
    Jhole[i][0] = Jhole[i][1]
    Jtotal[i][n_max-1] = Jtotal[i][n_max-2]
    Jelec[i][n_max-1] = Jelec[i][n_max-2]
    Jhole[i][n_max-1] = Jhole[i][n_max-2]

Ec[0] = Ec[1]
Ec[n_max-1] = Ec[n_max-2]
xx1[n_max-1] = xx1[n_max-2] + dx*Ldi*1e4
el_field1[0] = el_field1[1]
el_field2[0] = el_field2[1]
el_field1[n_max - 1] = el_field1[n_max - 2]
el_field2[n_max - 1] = el_field2[n_max - 2]
nf = [x*ni for x in n]
pf = [x*ni for x in p]
ro[0] = ro[1]
ro[n_max - 1] = ro[n_max - 2]

# Calculate Quasi Fermi Level - Efn Efp
for i in range(n_max):
    Ei[i]   = Ec[i] - 0.56
    Efn[i]  = Ei[i] + Vt*math.log(nf[i]/ni)
    Efp[i]  = Ei[i] - Vt*math.log(pf[i]/ni)
Ev = [(x - 1.12) for x in Ec]

f1 = open("quasi_py.txt", "w")
f1.write("x [um]\tEc\tEv\tEi\tEfn\tEfp\n")
for i in range(n_max):
    f1.write(str(xx1[i])+'\t'+str(Ec[i])+'\t'+str(Ev[i])+'\t'+str(Ei[i])+'\t'+str(Efn[i])+'\t'+str(Efp[i])+'\n')
f1.close()

f1 = open("pote_py.txt", "w")
f1.write("x [um]\tPotential [eV]\n")
for i in range(n_max):
    f1.write(str(xx1[i])+'\t'+str(Vt*fi[i])+'\n')
f1.close();

f1 = open("elec_py.txt", "w")
f1.write("x [um]\tElectric Field [V/cm]\tElectric Field [V/cm]\n")
for i in range(n_max):
    f1.write(str(xx1[i])+'\t'+str(el_field1[i])+'\t'+str(el_field2[i])+'\n')
f1.close()

f1 = open("car_dens_py.txt", "w")
f1.write("x [um]\tElectron Densities [1/cm^3]\tHole Densities [1/cm^3]\n")
for i in range(n_max):
    f1.write(str(xx1[i])+'\t'+str(nf[i])+'\t'+str(pf[i])+'\n')
f1.close()

f1 = open("charge_py.txt", "w")
f1.write("x [um]\tTotal Charge Density [C/cm^3]\n")
for i in range(n_max):
    f1.write(str(xx1[i])+'\t'+str(q*ro[i])+'\n')
f1.close()

f1 = open("curr_tot_ele_hol_py.txt", "w");
f1.write("Vplot\tTotal Current Density [Amp/cm^2]\tElectron Current Density [Amp/cm^2]\tHole Current Density [Amp/cm^2]\n");
for i in range(74):
    f1.write(str(Vplot[i])+'\t'+str(Jtotal[i][1])+'\t'+str(Jelec[i][1])+'\t'+str(Jhole[i][1])+'\n')
f1.close()

f1 = open("curr_pos_py.txt", "w")
f1.write("x [um]\tTotal Current Density [Amp/cm^2]\n")
for i in range(n_max):
    f1.write(str(xx1[i])+'\t'+str(Jtotal[73][i])+'\n')
f1.close()
