#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

double abs1(double x);    // Return an absolute value
double abs1(double x) {
    if (x > 0) return x;
    if (x < 0) return -x;
}
double BER(double x);     // Bernoulli Functiton
double BER(double x) {
    int flag_sum = 0, i = 0;
    double temp_term = 0.0, sum = 0.0;
    if (x > 1.0E-10) return x*exp(-x)/(1-exp(-x));
    else if (x < 0.0 && abs1(x) > 1.0E-10) return x/(exp(x)-1);
    else if (x == 0.0) return 1.0;
    else {
        temp_term = 1;
        sum = temp_term;
        i = 0.0;
        while (flag_sum == 0) {
            i += 1;
            temp_term = temp_term*x/(i+1);
            if (abs1(temp_term) < 1.0e-14) flag_sum = 1;
            else sum += temp_term;
        }
        return 1/sum;
    }
}

typedef struct rnode {
    double dop0;      //doping concentration in each grid point
    double fi;        //potential of each grid point
    double n;         //electron concentration normalized with ni intrinsic electron concentration
    double p;         //hole concentration normalized with ni
    double a, b, c, f, alpha, beta, v;  //coefficients of Poisson equation
    double delta;                       // Difference between elements of two matrixs
    double xx1, Ei, Ec, Ev, ro, el_field1, el_field2;   // Position grid point, intr/cond/val band, char density, Elect field.
    double nf, pf;    //electron and hole density
    double Ef, Efn, Efp;  // Quasi electron and hole fermi level
    double mup, mun;      // electron and hole mobility at each grid
    double an, bn, cn, fn, ap, bp, cp, fp;   // coefficients of DD-continuity equation
    double alphan, betan, alphap, betap;     // LU decomposition of current equations
    double vn, vp;                           // electron and hole potential
} Rnode;

typedef struct cur {
    double Jnim1by2, Jnip1by2;
    double Jpim1by2, Jpip1by2;
    double Jelec, Jhole, Jtotal;
} Cur;

typedef struct vpro {
    float Vplot, Jelec, Jhole, Jtotal;
} Vpro;

int main() {
    /* Defining the Fundamental and Material Constants */
    double q = 1.602E-19;        // C or [J/eV]
    double kb = 1.38E-23;         // [J/K]
    double eps = 1.05E-12;         // This includes the eps  = 11.7 for Si [F/cm]
    double T = 300.0;              // [K]
    double ni = 1.5E10;           // Intrinsic carrier concentration [1/cm^3]
    double Vt = kb*T/q;           // Thermal Energy [eV] equals to 25 mV thermal voltage [V]
    double RNc = 2.8E19;           // This is 2.8e20 in the FORTRAN file
    double TAUN0 = 0.1E-6;           // Electron SRH life time
    double TAUP0 = 0.1E-6;           // Hole SRH life time
    double mun0   = 1500.0;            // Electron Mobility in cm2/V-s
    double mup0   = 1000.0;             // Hole Mobility in cm2/V-s
    double dEc = Vt*log(RNc/ni);


    /* Define Doping Values */
    double Na = 1.0E16;             // [1/cm^3]
    double Nd = 1.0E17;             // [1/cm^3]


    /* Calculate relevant parameters for the simulation */
    double Vbi = Vt*log(Na*Nd/(ni*ni));                 // Built-in Potential Energy [eV] at zero bias (open circuit)
    double W   = sqrt(2*eps*(Na+Nd)*Vbi/(q*Na*Nd));     // Depletion width [cm]
    double Wn  = W*sqrt(Na/(Na+Nd));                    // N-side depletion width [cm]
    double Wp  = W*sqrt(Nd/(Na+Nd));                    // P-side depletion width [cm]
    double E_p = q*Nd*Wn/eps;                           // Maximum electric field [V/cm] at the junction
    double Ldn = sqrt(eps*Vt/(q*Nd));                   // Debye length for N-side
    double Ldp = sqrt(eps*Vt/(q*Na));                   // Debye length for P-side
    double Ldi = sqrt(eps*Vt/(q*ni));                   // Debye Length for Intrisic Si


    /* Setting the size of the simulated pn diode based on the analytical results from the width of the
       depletion regions for a simple pn-diode %*/
    double x_max = 0.0;
    if(x_max < Wn) x_max = Wn;
    if(x_max < Wp) x_max = Wp;
    x_max = 20*x_max;


    /* Setting the uniform grid size based (dx) on the extrinsic Debye lengths */
    double dx = Ldn;
    if(dx > Ldp) dx=Ldp;
    dx = dx/20;


    // Calculate the required number of grid points and normalize dx by Intrinsic Debye length (Ldi)
    int n_max = x_max/dx +1;
    dx = dx/Ldi;    // Renormalize lengths with Ldi

    /* Set up the doping C(x)=Nd(x)-Na(x) that is normalized with ni */
    Rnode *rnd0 = NULL, *rnd = NULL, *rnd1 = NULL;
    Cur *Jnp0 = NULL, *Jnp = NULL, *Jnp1 = NULL;
    rnd0 = new Rnode[n_max];
    rnd1 = rnd0 + n_max -1;
    Jnp0 = new Cur[n_max];
    Jnp1 = Jnp0 + n_max -1;
    for (rnd = rnd0; rnd <= rnd1; rnd++) {
        if (rnd - rnd0 < n_max/2) rnd->dop0 = -Na/ni;   // P-side from 0 to x_max/2
        else rnd->dop0 = Nd/ni;   // N-side from x_max/2 to x_max
    }

    /* Initialize the normalized potential based on the requirement of charge neutrality throughout the whole structure */
    double zz = 0.0, xx = 0.0, delta_acc = 0.0;
    for (rnd = rnd0; rnd <= rnd1; rnd++) {
        zz = 0.5*(rnd->dop0);
        if (zz > 0.0) xx = zz*(1 + sqrt(1+1/(zz*zz)));
        if (zz < 0.0) xx = zz*(1 - sqrt(1+1/(zz*zz)));     // Taylor series: sqrt(1+x) = 1 + x/2
        rnd->fi = log(xx);           // \phi = (kT/q)ln(Nd/ni)  or  \phi = (kT/q)ln(Na/ni) in kT/q unit
        rnd->n = xx;                 //  n = Nd in N-side  and  n = ni*2/Na in P-side in ni unit
        rnd->p = 1/xx;               //  p = Na in P-side  and  n = ni*2/Nd in P-side in ni unit
    }
    delta_acc = 1.0E-5;               // Preset the Tolerance


    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%                                                                      %%
    %%               EQUILIBRIUM  SOLUTION PART BEGINS                      %%
    %%                                                                      %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %(A) Define the elements of the coefficient matrix for the internal nodes and
    %    initialize the forcing function*/
    double dx2 = dx*dx;
    for (rnd = rnd0; rnd <= rnd1; rnd++) {   // Setup the coefficients of discritized Poisson equation
        rnd->a = 1/dx2;
        rnd->c = 1/dx2;
        rnd->b = -(2/dx2+exp(rnd->fi)+exp(-rnd->fi));
        rnd->f = exp(rnd->fi) - exp(-rnd->fi) - rnd->dop0 - (rnd->fi)*(exp(rnd->fi)+exp(-rnd->fi));
    }


    /*%(B) Define the elements of the coefficient matrix and initialize the forcing
    %    function at the ohmic contacts*/
    rnd0->a = 0.0;                      // Setup the Dirichlet boundary condition for the boundary nodes
    rnd0->c = 0.0;
    rnd0->b = 1.0;
    rnd0->f = rnd0->fi;
    rnd1->a = 0.0;
    rnd1->c = 0.0;
    rnd1->b = 1.0;
    rnd1->f = rnd1->fi;


    /*%(C)  Start the iterative procedure for the solution of the linearized Poisson
      %     equation using LU decomposition method:*/

    double temp = 0.0;
    double delta_max = 0.0;
    int flag_conv = 0;		           //% convergence of the Poisson loop
    int k_iter= 0, K = 0;
    while(flag_conv == 0) {
        k_iter = k_iter + 1;             // count iteration cycles
        rnd0->alpha = rnd0->b;
        for (rnd =  rnd0 + 1; rnd <= rnd1; rnd++) {            // Set up L and U matrix
            rnd->beta = (rnd->a)/((rnd-1)->alpha);
            rnd->alpha = rnd->b - (rnd->beta)*((rnd-1)->c);
        }

        /*% Solution of Lv = f %*/
        rnd0->v = rnd0->f;
        for (rnd = rnd0 + 1; rnd <= rnd1; rnd++) {
                rnd->v = rnd->f - (rnd->beta)*((rnd-1)->v);    // Calculate v array using L matrix and fi array from the previous cycle
        }

        /*% Solution of U*fi = v %*/
        temp = (rnd1->v)/(rnd1->alpha);
        rnd1->delta = temp - rnd1->fi;
        rnd1->fi = temp;
        for (rnd = rnd1 -1; rnd >= rnd0; rnd--){              // Calculate the new fi using U matrix and v array.
            temp = (rnd->v - (rnd->c)*((rnd + 1)->fi))/(rnd->alpha);
            rnd->delta = temp - rnd->fi;
            rnd->fi = temp;
        }
        delta_max = 0.0;
        for (rnd = rnd0; rnd <= rnd1; rnd++) {     // Finding the maximum element difference between new and old fi
            if (rnd->delta > 0.0) xx = rnd->delta;
            else xx = -rnd->delta;
            if (xx > delta_max) {delta_max = xx; K = rnd - rnd0;}
        }

        /*% Test convergence and recalculate forcing function and
        % central coefficient b if necessary*/
        cout << delta_max << endl;
        if (delta_max < delta_acc) flag_conv = 1;    // Check if the maximum difference is less or greater than threshold
        else
        for (rnd = rnd0 + 1; rnd <= rnd1 -1; rnd++) {//If difference is greater then the new b and f would be calculated and used in new iter
            rnd->b = -(2/dx2 + exp(rnd->fi) + exp(-rnd->fi));
            rnd->f = exp(rnd->fi) - exp(-rnd->fi) - rnd->dop0 - rnd->fi*(exp(rnd->fi) + exp(-rnd->fi));
        }
    }
// End of Poisson-Boltzmann solver


    /*Calculating electrostatic properties of the junction for plotting*/
    rnd0->xx1 = dx*1e4;
    for (rnd = rnd0 + 1; rnd <= rnd1 - 1; rnd++) {
        rnd->Ec = dEc - Vt*(rnd->fi);// Calculate electric field, charge carrier densities and conduction band profile based on calculated fi
        rnd->ro = -ni*(exp(rnd->fi) - exp(-rnd->fi) - rnd->dop0);
        rnd->el_field1 = -((rnd + 1)->fi - rnd->fi)*Vt/(dx*Ldi);
        rnd->el_field2 = -((rnd + 1)->fi - (rnd - 1)->fi)*Vt/(2*dx*Ldi);
        rnd->n = exp(rnd->fi);
        rnd->p = exp(-rnd->fi);
        rnd->xx1 = (rnd - 1)->xx1 + dx*Ldi*1e4;
    }
    rnd0->Ec = (rnd0 + 1)->Ec;
    rnd0->ro = (rnd0 + 1)->ro;
    rnd0->el_field1 = (rnd0 + 1)->el_field1;
    rnd0->el_field2 = (rnd0 + 1)->el_field2;
    rnd1->Ec = (rnd1 - 1)->Ec;
    rnd1->ro = (rnd1 - 1)->ro;
    rnd1->el_field1 = (rnd1 - 1)->el_field1;
    rnd1->el_field2 = (rnd1 - 1)->el_field2;
    rnd1->xx1 = (rnd1 - 1)->xx1 + dx*Ldi*1e4;
    for (rnd = rnd0; rnd <= rnd1; rnd++) {
        rnd->nf = (rnd->n)*ni;
        rnd->pf = (rnd->p)*ni;
    }

    // Dumping out the results
    ofstream myfile;
    myfile.open("els_plot_cpp.txt");
    myfile << "x [um]\tPotential [V]\tElectric Field1 [V/cm]\tElectric Field2 [V/cm]\tElectron Densities [1/cm^3]\tHole Densities [1/cm^3]\tTotal Charge Density [C/cm^3]\tConduction Band Energy (eV)\n";
    for (rnd = rnd0; rnd <= rnd1; rnd++) {
         myfile << rnd->xx1 << "\t" << Vt*(rnd->fi) << "\t" << rnd->el_field1 << "\t" << rnd->el_field2 << "\t" << rnd->nf << "\t" << rnd->pf << "\t" << rnd->ro << "\t" << rnd->Ec << "\n";
    }
    myfile.close();
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%                 END OF EQUILIBRIUM  SOLUTION PART                    %%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%                                                                      %%
      %%               NON-EQUILIBRIUM  SOLUTION PART BEGINS                  %%
      %%                                                                      %%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%              1. Calculate Low filed mobility                         %%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      %%  Prameters for Low field mobility calculation %%*/
    double TL = 300;                    //% Temp in Kelvin
    double N  = Na + Nd;                //% Local (total) impurity concentration
    double BETAN = 2.0;                 // Field dependency parameter of electron mobility (unitless)
    double BETAP = 1.0;                 // Field Dependency parameter of hole mobility (unitless)

    double VSATN = (2.4E7) / (1 + 0.8*exp(TL/600));  // Saturation Velocity of Electrons
    double VSATP = VSATN;                               // Saturation Velocity of Holes

    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%   2. Start the main Loop to increment the Anode voltage by Vt=KbT/q  %%
    %%      till it reaches 0.625V.                                         %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

    double VA = 0.0, Each_Step = 0.33*Vt, Vf = 0.625, pdeno = 0.0, ndeno = 0.0;
    double munim1by2 = 0.0, munip1by2 = 0.0, mupim1by2 = 0.0, mupip1by2 = 0.0;
    double tempn = 0.0, tempp = 0.0, deltan = 0.0, deltap = 0.0;
    int vindex= -1, Total_Steps = Vf/(0.33*Vt), flag_conv2 = 0, k_itern = 0, kk = 0;
    Vpro **IV0 = NULL, **IV = NULL, *IVptr = NULL;
    remove("Current.txt");

    Vpro  **dptrc, *sptrc, *sptrc0;
    IV0 = new Vpro*[74];
    for (IV = IV0; IV < IV0 + 74; IV++) *IV = new Vpro[5];


    for (VA = 0.0, IV = IV0; VA < Vf; VA+=0.33*Vt, IV++) {                //% Start VA increment loop

        printf("VA = %f\n",VA);
        vindex = vindex + 1;
        (*IV)->Vplot = VA;

        rnd0->fi = rnd0->fi + VA;            //% Apply potential to Anode (1st node)

        flag_conv2 = 0;		           //% Convergence of the Poisson loop
        k_itern= 0;

        //%% Initialize the First and Last Node for Poisson's eqn
        rnd0->a = 0;
        rnd0->c = 0;
        rnd0->b = 1;
        rnd0->f = rnd0->fi;
        rnd1->a = 0;
        rnd1->c = 0;
        rnd1->b = 1;
        rnd1->f = rnd1->fi;


    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 3. Start the Poisson equation solver loop to calculate the        %%
    %%    potential for each Anode voltage increase                      %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
        while(flag_conv2 == 0) {     //Start Poisson's eqn
            k_itern += 1;


        /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% 3.1 . Calculate Field Dependant Mobility for each value of 'fi'   %%
        %%       at each node point of the PN diode.                         %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

        //%%% To test with Constant Mobility without field dependancy.
        //% %         for i = 1:n_max           % Start Loop for Field Dep Mobility
        //% %             mup(i) = mup0;
        //% %             mun(i) = mun0;
        //% %         end
        //%% Calculate the Electric Field at each Node
            for (rnd = rnd0 + 1; rnd <= rnd1 -1; rnd++) {
                rnd->Ef = abs1(rnd->fi - (rnd + 1)->fi)*Vt/(dx*Ldi);
            }

            rnd0->Ef = (rnd0 + 1)->Ef;
            rnd1->Ef = (rnd1 -1 )->Ef;

        //%% Calculate the Field Dependant Mobility at each Node
            for (rnd = rnd0; rnd <= rnd1; rnd++) {
                pdeno = pow((mup0 * rnd->Ef / VSATP) , BETAP);
                rnd->mup = mup0 * pow((1/(1 + pdeno)) , (1/BETAP));
                ndeno = pow((mun0 * rnd->Ef / VSATN) , BETAN);
                rnd->mun = mun0 * pow((1/(1 + ndeno)) , (1/BETAN));
                //if (rnd - rnd0 < 10) printf("mup = %e\tmun = %e\n", rnd->mup, rnd->mun);
            }

            rnd0->mup = (rnd0 + 1)->mup;
            rnd1->mup = (rnd1 - 1)->mup;

            rnd0->mun = (rnd0 + 1)->mun;
            rnd1->mun = (rnd1 - 1)->mun;

        /*%%%%%%%%%%% END of FIELD Dependant Mobility Calculation %%%%%%%%%%%*/
        /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% 3.2 Solve Continuity Equation for Electron and Holes using LU Decomposition %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

        //%(A) Define the elements of the coefficient matrix and initialize the forcing
        //%    function at the ohmic contacts for ELECTRON and HOLE Continuity Eqns


            rnd0->an = 0;              //%Co-ef for electron at Anode
            rnd0->bn = 1;              //%Co-ef for electron at Anode
            rnd0->cn = 0;              //%Co-ef for electron at Anode
            rnd0->ap = 0;              //%Co-ef for hole     at Anode
            rnd0->bp = 1;              //%Co-ef for hole     at Anode
            rnd0->cp = 0;              //%Co-ef for hole     at Anode
            rnd0->fn = rnd0->n;
            rnd0->fp = rnd0->p;


            rnd1->an = 0;          //%Co-ef for electron at Cathode
            rnd1->bn = 1;          //%Co-ef for electron at Cathode
            rnd1->cn = 0;          //%Co-ef for electron at Cathode
            rnd1->ap = 0;          //%Co-ef for hole     at Cathode
            rnd1->bp = 1;          //%Co-ef for hole     at Cathode
            rnd1->cp = 0;          //%Co-ef for hole     at Cathode
            rnd1->fn = rnd1->n;
            rnd1->fp = rnd1->p;



        //%(B) Define the elements of the coefficient matrix for the internal nodes and
        //%    initialize the forcing function




            for (rnd = rnd0 + 1; rnd <= rnd1 - 1; rnd++) {
                munim1by2 = ((rnd - 1)->mun + rnd->mun)/2;      // Electron mobility between node (i-1) and i; D_(i-1/2)^n
                munip1by2 = (rnd->mun + (rnd + 1)->mun)/2;      // Electron mobility between node (i) and (i+1); D_(i+1/2)^n
                mupim1by2 = ((rnd - 1)->mup + rnd->mup)/2;      // Hole mobility between node (i-1) and i; D_(i-1/2)^n
                mupip1by2 = (rnd->mup + (rnd + 1)->mup)/2;      // Hole mobility between node (i) and (i+1); D_(i+1/2)^n

            //%% Co-efficients for HOLE Continuity eqn
                rnd->cp = mupip1by2 * BER(rnd->fi - (rnd + 1)->fi);
                rnd->ap = mupim1by2 * BER(rnd->fi - (rnd - 1)->fi);
                rnd->bp = -( mupim1by2 * BER((rnd - 1)->fi - rnd->fi) + mupip1by2 * BER((rnd + 1)->fi - rnd->fi));
            //%% Co-efficients for ELECTRON Continuity eqn
                rnd->cn = munip1by2 * BER((rnd + 1)->fi - rnd->fi);
                rnd->an = munim1by2 * BER((rnd - 1)->fi - rnd->fi);
                rnd->bn = -( munim1by2 * BER(rnd->fi - (rnd - 1)->fi) + munip1by2 * BER(rnd->fi - (rnd + 1)->fi));
            //%% Forcing Function for ELECTRON and HOLE Continuity eqns
                rnd->fn = (Ldi*Ldi*dx2/Vt) * ( (rnd->p)*(rnd->n) - 1 ) / ( TAUP0*(rnd->n + 1) + TAUN0*(rnd->p + 1)); // Electron Recom
                rnd->fp = (Ldi*Ldi*dx2/Vt) * ( (rnd->p)*(rnd->n) - 1 ) / ( TAUP0*(rnd->n + 1) + TAUN0*(rnd->p + 1)); // Hole Recom
            }


        /*(C)  Start the iterative procedure for the solution of the linearized Continuity
        %     equation for "ELECTRONS" using LU decomposition method:*/


            rnd0->alphan = rnd0->bn;
            for (rnd = rnd0 + 1; rnd <= rnd1; rnd++) {
                rnd->betan = rnd->an/(rnd - 1)->alphan;
                rnd->alphan = rnd->bn - rnd->betan*(rnd - 1)->cn;
            }

        //% Solution of Lv = f %

            rnd0->vn = rnd0->fn;
            for (rnd = rnd0 + 1; rnd <= rnd1; rnd++) {
                rnd->vn = rnd->fn - rnd->betan*(rnd - 1)->vn;
            }

        //% Solution of U*fi = v %
            tempn = rnd1->vn/rnd1->alphan;
            rnd1->n=tempn;
            for (rnd = rnd1 - 1; rnd >= rnd0; rnd--) {       //%delta%
                tempn = (rnd->vn - rnd->cn*(rnd + 1)->n)/rnd->alphan;
                rnd->n = tempn;
            }

       //%%%%%%%%%%%%%%%%%%%%%%% END of ELECTRON Continuty Solver %%%%%%%%%%%


       // %(D)  Start the iterative procedure for the solution of the linearized Continuity
        //%     equation for "HOLES" using LU decomposition method:

            rnd0->alphap = rnd0->bp;
            for (rnd = rnd0 + 1; rnd <= rnd1; rnd++) {
                rnd->betap = rnd->ap/(rnd-1)->alphap;
                rnd->alphap = rnd->bp - rnd->betap*(rnd - 1)->cp;
            }

        //% Solution of Lv = f %

            rnd0->vp = rnd0->fp;
            for (rnd = rnd0 + 1; rnd <= rnd1; rnd++) {
                rnd->vp = rnd->fp - rnd->betap*(rnd - 1)->vp;
            }

        //% Solution of U*fi = v %

            tempp = rnd1->vp/rnd1->alphap;
            rnd1->p=tempp;
            for (rnd = rnd1 - 1; rnd >= rnd0; rnd--) {       //%delta%
                tempp = (rnd->vp - rnd->cp*(rnd + 1)->p)/rnd->alphap;
                rnd->p = tempp;
            }

       //%%%%%%%%%%%%%%%%%%%%%%% END of HOLE Continuty Solver %%%%%%%%%%%


       /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %% 3.3 Calculate potential fi again with new values of "n" and "p"%%
       %%     and check for convergence                                  %%
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       % Recalculate forcing function and central coefficient b for fi*/

        for (rnd = rnd0 + 1; rnd <= rnd1 -1; rnd++) {
            rnd->b = -(2/dx2 + rnd->n + rnd->p);
            rnd->f = rnd->n - rnd->p - rnd->dop0 - (rnd->fi*(rnd->n + rnd->p));
        }

            /*% Solve for Updated potential given the new value of Forcing
            % Function using LU decomposition*/

        rnd0->alpha = rnd0->b;
        for (rnd = rnd0 + 1; rnd <= rnd1; rnd++) {
            rnd->beta = rnd->a/(rnd - 1)->alpha;
            rnd->alpha = rnd->b - rnd->beta*(rnd - 1)->c;
        }

        //% Solution of Lv = f %*/

        rnd0->v = rnd0->f;
        for (rnd = rnd0 + 1; rnd <= rnd1; rnd++) {
            rnd->v = rnd->f - rnd->beta*(rnd - 1)->v;
        }

        //% Solution of U*fi = v %

        temp = rnd1->v/rnd1->alpha;
        rnd1->delta = temp - rnd1->fi;
        rnd1->fi=temp;
        for (rnd = rnd1 - 1; rnd >= rnd0; rnd--) {      //%delta%
            temp = (rnd->v - rnd->c*(rnd + 1)->fi)/rnd->alpha;
            rnd->delta = temp - rnd->fi;
            rnd->fi = temp;
        }

        delta_max = 0.0;
        for (rnd = rnd0; rnd <= rnd1; rnd++) {
            xx = abs1(rnd->delta);
            if(xx > delta_max) delta_max=xx;
        }

        cout << delta_max << endl;
        if(delta_max < delta_acc) flag_conv2 = 1;
        }   // end of while loop

    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%                        CALCULATE CURRENT                             %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

        myfile.open("Current_cpp.txt");
        for (rnd = rnd0 + 1, Jnp = Jnp0 + 1; rnd <= rnd1 -1; rnd++, Jnp++) {
            Jnp0->Jnim1by2 = (q*(rnd->mun)*Vt/(dx*Ldi)) * ni*( (rnd->n)*BER(rnd->fi - (rnd - 1)->fi) - ((rnd - 1)->n)*BER((rnd - 1)->fi - rnd->fi) );   // Electron current from node (i-1) to (i) node at each applied voltage
            Jnp->Jnim1by2 = (q*(rnd->mun)*Vt/(dx*Ldi)) * ni*( ((rnd + 1)->n)*BER((rnd + 1)->fi - rnd->fi) - (rnd->n)*BER(rnd->fi - (rnd + 1)->fi) );       // Electron current from node (i) to (i+1) node at each applied voltage
            Jnp->Jelec = (Jnp0->Jnim1by2 + Jnp->Jnim1by2)/2;   // Electron current density on each node for each applied voltage

            Jnp->Jpim1by2 = (q*(rnd->mup)*Vt/(dx*Ldi)) * ni*( (rnd->p)*BER((rnd - 1)->fi - rnd->fi) - ((rnd - 1)->p)*BER(rnd->fi - (rnd - 1)->fi) );       // Hole current from node (i-1) to (i) node at each applied voltage
            Jnp->Jpip1by2 = (q*(rnd->mup)*Vt/(dx*Ldi)) * ni*( ((rnd + 1)->p)*BER(rnd->fi - (rnd + 1)->fi) - (rnd->p)*BER((rnd + 1)->fi - rnd->fi) );       // Hole current from node (i) to (i+1) node at each applied voltage
            Jnp->Jhole = (Jnp->Jpim1by2 + Jnp->Jpip1by2)/2;  // Hole current density on each node for each applied voltage

            Jnp->Jtotal = Jnp->Jelec + Jnp->Jhole;  // Total current density on each node for each applied voltage

            if (rnd - rnd0 == 1) {
                (*IV + 1)->Jelec = (*IV + 2)->Jelec = Jnp->Jelec;
                (*IV + 1)->Jhole = (*IV + 2)->Jhole = Jnp->Jhole;
                (*IV + 1)->Jtotal = (*IV + 2)->Jtotal = Jnp->Jtotal;
            }

            if (rnd1 - rnd == 1) {
                (*IV + 3)->Jelec = (*IV + 4)->Jelec = Jnp->Jelec;
                (*IV + 3)->Jhole = (*IV + 4)->Jelec = Jnp->Jhole;
                (*IV + 3)->Jtotal = (*IV + 4)->Jelec = Jnp->Jtotal;
            }

            myfile << Jnp->Jelec << "\t" << Jnp->Jhole << "\t" << Jnp->Jtotal;
            if ((rnd - rnd0)%3 == 0) myfile << "\n";
            else myfile << "\t";
        }
        myfile.close();
    }

        rnd0->xx1 = dx*1e4;
    for (rnd = rnd0 + 1; rnd <= rnd1 - 1; rnd++) {
        rnd->Ec = dEc - Vt*rnd->fi;     //%Values from the second Node%
        rnd->ro = -ni*(rnd->n - rnd->p - rnd->dop0);
        rnd->el_field1 = -((rnd + 1)->fi - rnd->fi)*Vt/(dx*Ldi);
        rnd->el_field2 = -((rnd + 1)->fi - (rnd - 1)->fi)*Vt/(2*dx*Ldi);
        rnd->xx1 = (rnd - 1)->xx1 + dx*Ldi*1e4;
    }

    rnd0->Ec = (rnd0 + 1)->Ec;
    rnd1->Ec = (rnd1 - 1)->Ec;
    rnd1->xx1 = (rnd1 - 1)->xx1 + dx*Ldi*1e4;
    rnd0->el_field1 = (rnd0 + 1)->el_field1;
    rnd0->el_field2 = (rnd0 + 1)->el_field2;
    rnd1->el_field1 = (rnd - 1)->el_field1;
    rnd1->el_field2 = (rnd1 - 1)->el_field2;
    for (rnd = rnd0; rnd <= rnd1; rnd++) {
        rnd->nf = (rnd->n)*ni;
        rnd->pf = (rnd->p)*ni;
    }
    rnd0->ro = (rnd0 + 1)->ro;
    rnd1->ro = (rnd1 - 1)->ro;

    //%% Calculate Quasi Fermi Level - Efn Efp
    for (rnd = rnd0; rnd <= rnd1; rnd++) {
        rnd->Ei   = rnd->Ec - 0.56;
        rnd->Efn  = rnd->Ei + Vt*log(rnd->nf/ni);
        rnd->Efp  = rnd->Ei - Vt*log(rnd->pf/ni);
        rnd->Ev = rnd->Ec - 1.12;
    }

    myfile.open("quasi_cpp.txt");
    myfile << "x [um]\tEc\tEv\tEi\tEfn\tEfp\n";
    for (rnd = rnd0; rnd <= rnd1; rnd++)
        myfile << rnd->xx1 << "\t" << rnd->Ec << "\t" << rnd->Ev << "\t" << rnd->Ei << "\t" << rnd->Efn << "\t" << rnd->Efp << "\n";
    myfile.close();

    myfile.open("pote_cpp.txt");
    myfile << "x [um]\tPotential [eV]\n";
    for (rnd = rnd0; rnd <= rnd1; rnd++) myfile << rnd->xx1 << "\t" << Vt*(rnd->fi) << "\n";
    myfile.close();

    myfile.open("elec_cpp.txt");
    myfile << "x [um]\tElectric Field [V/cm]\tElectric Field [V/cm]\n";
    for (rnd = rnd0; rnd <= rnd1; rnd++)
        myfile << rnd->xx1 << "\t" << rnd->el_field1 << "\t" << rnd->el_field2 << "\n";
    myfile.close();

    myfile.open("car_dens_cpp.txt");
    myfile << "x [um]\tElectron Densities [1/cm^3]\tHole Densities [1/cm^3]\n";
    for (rnd = rnd0; rnd <= rnd1; rnd++)
        myfile << rnd->xx1 << "\t" << rnd->nf << "\t" << rnd->pf << "\n";
    myfile.close();

    myfile.open("charge_cpp.txt");
    myfile << "x [um]\tTotal Charge Density [C/cm^3]\n";
    for (rnd = rnd0; rnd <= rnd1; rnd++)
        myfile <<  rnd->xx1 << "\t" << q*rnd->ro << "\n";
    myfile.close();

    myfile.open("curr_tot_ele_hol_cpp.txt");
    myfile << "Vplot\tTotal Current Density [Amp/cm^2]\tElectron Current Density [Amp/cm^2]\tHole Current Density [Amp/cm^2]\n";
    for (IV = IV0; IV - IV0 < 73; IV++)
        myfile << (*IV)->Vplot << "\t" << (*IV + 3)->Jtotal << (*IV + 3)->Jelec << "\t" << (*IV + 3)->Jhole << "\n";
    myfile.close();

    myfile.open("curr_pos_cpp.txt");
    myfile << "x [um]\tTotal Current Density [Amp/cm^2]\n";
    for (rnd = rnd0, Jnp = Jnp0; rnd <= rnd1; rnd++, Jnp++)
        myfile <<  rnd->xx1 << "\t" << Jnp->Jtotal << "\n";
    myfile.close();

    for(IV = IV0; IV < IV0 + 74; IV++) {
        delete [] *IV;
    }
    delete [] IV0;
    delete[] rnd0;
    delete[] Jnp0;
    return 0;
}
