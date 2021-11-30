#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

double abs1(double x);
double abs1(double x) {
    if (x > 0) return x;
    if (x < 0) return -x;
}
double BER(double x);
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
    double delta;                       // Difference of a matrix elements
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
    double Vt = kb*T/q;           // Thermal Energy [eV] equals to 25 mV termal voltage [V]
    double RNc = 2.8E19;           // This is 2.8e20 in the FORTRAN file
    double TAUN0 = 0.1E-6;           // Electron SRH life time
    double TAUP0 = 0.1E-6;           // Hole SRH life time
    double mun0   = 1500.0;            // Electron Mobility in cm2/V-s
    double mup0   = 1000.0;             // Hole Mobility in cm2/V-s
    double dEc = Vt*log(RNc/ni);
    //printf("q = %e\tkb = %e\teps = %e\tT = %f\tni = %e\tVt = %f\tRNc = %e\tTAUN0 = %e\tTAUP0 = %e\tmun0 =%f\tmunp0 =%f\tdEc = %f\n", q, kb, eps, T, ni, Vt, RNc, TAUN0, TAUP0, mun0, mup0, dEc);


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
    //printf("Vbi = %f\tW = %e\tWn = %e\tWp = %e\tWone = %e\tE_p = %e\tLdn = %e\tLdp = %e\tLdi = %e\n",Vbi,W,Wn,Wp,Wone,E_p,Ldn,Ldp,Ldi);



    /* Setting the size of the simulated pn diode based on the analytical results from the width of the
       depletion regions for a simple pn-diode %*/
    double x_max = 0.0;
    if(x_max < Wn) x_max = Wn;
    if(x_max < Wp) x_max = Wp;
    x_max = 20*x_max;


    /* Setting the grid size based (dx) on the extrinsic Debye lengths */
    double dx = Ldn;
    if(dx > Ldp) dx=Ldp;
    dx = dx/20;


    /* Calculate the required number of grid points and renormalize dx */
    int n_max = x_max/dx +1;
    dx = dx/Ldi;    // Renormalize lengths with Ldi
    //printf("n_max = %d\tx_max = %e\tdx = %e\n", n_max, x_max, dx);
        /* Set up the doping C(x)=Nd(x)-Na(x) that is normalized with ni */
    Rnode *rnd0 = NULL, *rnd = NULL, *rnd1 = NULL;
    Cur *Jnp0 = NULL, *Jnp = NULL, *Jnp1 = NULL;
    rnd0 = new Rnode[n_max];
    rnd1 = rnd0 + n_max -1;
    Jnp0 = new Cur[n_max];
    Jnp1 = Jnp0 + n_max -1;
    for (rnd = rnd0; rnd <= rnd1; rnd++) {
        if (rnd - rnd0 < n_max/2) rnd->dop0 = -Na/ni;
        else rnd->dop0 = Nd/ni;
        //if (rnd - rnd0 ==0 || rnd - rnd0 == 9653) printf("dop0 = %e\n", rnd->dop0);
    }

        /* Initialize the potential based on the requirement of charge neutrality throughout the whole structure */
    double zz = 0.0, xx = 0.0, delta_acc = 0.0;
    for (rnd = rnd0; rnd <= rnd1; rnd++) {
        zz = 0.5*(rnd->dop0);
        if (zz > 0.0) xx = zz*(1 + sqrt(1+1/(zz*zz)));
        if (zz < 0.0) xx = zz*(1 - sqrt(1+1/(zz*zz)));
        rnd->fi = log(xx);
        rnd->n = xx;
        rnd->p = 1/xx;
        //if (rnd - rnd0 == 0) printf("dop0 = %e\tzz = %e\txx = %e\tfi = %f\tn = %e\tp =%e\n",
                                    //rnd->dop0, zz, xx, rnd->fi, rnd->n, rnd->p);
        //if (rnd - rnd0 == 9653) printf("dop0 = %e\tzz = %e\txx = %e\tfi = %f\tn = %e\tp =%e\n",
                                       //rnd->dop0, zz, xx, rnd->fi, rnd->n, rnd->p);
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
    for (rnd = rnd0; rnd <= rnd1; rnd++) {
        rnd->a = 1/dx2;
        rnd->c = 1/dx2;
        rnd->b = -(2/dx2+exp(rnd->fi)+exp(-rnd->fi));
        rnd->f = exp(rnd->fi) - exp(-rnd->fi) - rnd->dop0 - (rnd->fi)*(exp(rnd->fi)+exp(-rnd->fi));
        //if (rnd - rnd0 == 0) printf("dx2 = %e\ta = %e\tc = %e\tb = %e\tf = %e\n", dx2, rnd->a, rnd->c, rnd->b, rnd->f);
        //if (rnd - rnd0 == 9653) printf("dx2 = %e\ta = %e\tc = %e\tb = %e\tf = %e\n", dx2, rnd->a, rnd->c, rnd->b, rnd->f);
    }


    /*%(B) Define the elements of the coefficient matrix and initialize the forcing
    %    function at the ohmic contacts*/
    rnd0->a = 0.0;
    rnd0->c = 0.0;
    rnd0->b = 1.0;
    rnd0->f = rnd0->fi;
    rnd1->a = 0.0;
    rnd1->c = 0.0;
    rnd1->b = 1.0;
    rnd1->f = rnd1->fi;
    //printf("a= %f\tc= %f\tb= %f\tf= %e\n",rnd0->a,rnd0->c,rnd0->b,rnd0->f);
    //printf("a= %f\tc= %f\tb= %f\tf= %e\n",rnd1->a,rnd1->c,rnd1->b,rnd1->f);


    /*%(C)  Start the iterative procedure for the solution of the linearized Poisson
      %     equation using LU decomposition method:*/

    double temp = 0.0;
    double delta_max = 0.0;
    int flag_conv = 0;		           //% convergence of the Poisson loo
    int k_iter= 0, K = 0;
    while(flag_conv == 0) {
        k_iter = k_iter + 1;
        rnd0->alpha = rnd0->b;
        for (rnd =  rnd0 + 1; rnd <= rnd1; rnd++) {
            rnd->beta = (rnd->a)/((rnd-1)->alpha);
            rnd->alpha = rnd->b - (rnd->beta)*((rnd-1)->c);
            //if (k_iter == 1 && rnd - rnd0 ==1) printf("alpha= %e\tbeta =%e\n",rnd->alpha,rnd->beta);
            //if (k_iter == 1 && rnd - rnd0 ==9653) printf("alpha= %e\tbeta =%e\n",rnd->alpha,rnd->beta);
        }

        /*% Solution of Lv = f %*/
        rnd0->v = rnd0->f;
        for (rnd = rnd0 + 1; rnd <= rnd1; rnd++) {
                rnd->v = rnd->f - (rnd->beta)*((rnd-1)->v);
                //if (k_iter == 1 && rnd - rnd0 ==1) printf("v= %e\n",rnd->v);
                //if (k_iter == 1 && rnd - rnd0 ==9653) printf("v= %e\n",rnd->v);
        }

        /*% Solution of U*fi = v %*/
        temp = (rnd1->v)/(rnd1->alpha);
        rnd1->delta = temp - rnd1->fi;
        rnd1->fi = temp;
        for (rnd = rnd1 -1; rnd >= rnd0; rnd--){
            temp = (rnd->v - (rnd->c)*((rnd + 1)->fi))/(rnd->alpha);
            rnd->delta = temp - rnd->fi;
            rnd->fi = temp;
            //if (k_iter == 1 && rnd - rnd0 ==0) printf("temp= %e\tdelta[%d]= %e\tfi[i]= %e\n",temp,rnd-rnd0,rnd->delta,rnd-rnd0,rnd->fi);
            //if (k_iter == 1 && rnd - rnd0 ==9652) printf("temp= %e\tdelta[%d]= %e\tfi[i]= %e\n",temp,rnd-rnd0,rnd->delta,rnd-rnd0,rnd->fi);
        }
        delta_max = 0.0;
        for (rnd = rnd0; rnd <= rnd1; rnd++) {
            if (rnd->delta > 0.0) xx = rnd->delta;
            else xx = -rnd->delta;
            if (xx > delta_max) {delta_max = xx; K = rnd - rnd0;}
            //if (rnd - rnd0 == 9653) printf("delta_max = %e\tK = %d\n", delta_max, K);
            //if (k_iter == 1 && rnd - rnd0 == 4822) printf("delta = %e\txx = %e\tK = %d\n",rnd->delta,xx,rnd-rnd0);
            //if (k_iter == 1 && rnd - rnd0 == 4824) printf("delta = %e\txx = %e\tK = %d\n",rnd->delta,xx,rnd-rnd0);
            //if (k_iter == 1 && rnd - rnd0 == 4826) printf("delta = %e\txx = %e\tK = %d\n",rnd->delta,xx,rnd-rnd0);
        }


        /*%delta_max=max(abs(delta));
        % Test convergence and recalculate forcing function and
        % central coefficient b if necessary*/
        cout << delta_max << endl;
        if (delta_max < delta_acc) flag_conv = 1;
        else
        for (rnd = rnd0 + 1; rnd <= rnd1 -1; rnd++) {
            rnd->b = -(2/dx2 + exp(rnd->fi) + exp(-rnd->fi));
            rnd->f = exp(rnd->fi) - exp(-rnd->fi) - rnd->dop0 - rnd->fi*(exp(rnd->fi) + exp(-rnd->fi));
            //if (k_iter == 1 && rnd - rnd0 == 4822) printf("i = %d\tb = %e\tf = %e\n", rnd - rnd0, rnd->b, rnd->f);
            //if (k_iter == 1 && rnd - rnd0 == 4824) printf("i = %d\tb = %e\tf = %e\n", rnd - rnd0, rnd->b, rnd->f);
            //if (k_iter == 1 && rnd - rnd0 == 4826) printf("i = %d\tb = %e\tf = %e\n", rnd - rnd0, rnd->b, rnd->f);
        }
    }


    /*Calculating electrostatic properties of the junction for plotting*/
    rnd0->xx1 = dx*1e4;
    for (rnd = rnd0 + 1; rnd <= rnd1 - 1; rnd++) {
        rnd->Ec = dEc - Vt*(rnd->fi);
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

    delete[] rnd0;
    delete[] Jnp0;
    return 0;
}
