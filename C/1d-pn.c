#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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
    double n;         //electron concentration
    double p;         //hole concentration
    double a, b, c, f, alpha, beta, v;  //coefficients of Poisson equation
    double delta;
    double xx1, Ei, Ec, Ev, ro, el_field1, el_field2;
    double nf, pf;
    double Ef, Efn, Efp;
    double mup, mun;
    double an, bn, cn, fn, ap, bp, cp, fp;
    double alphan, betan, alphap, betap;
    double vn, vp;
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
    double Vt = kb*T/q;           // [eV]
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
    double Vbi = Vt*log(Na*Nd/(ni*ni));
    double W   = sqrt(2*eps*(Na+Nd)*Vbi/(q*Na*Nd));     // [cm]
    double Wn  = W*sqrt(Na/(Na+Nd));                    // [cm]
    double Wp  = W*sqrt(Nd/(Na+Nd));                    // [cm]
    double Wone = sqrt(2*eps*Vbi/(q*Na));               // [cm]
    double E_p = q*Nd*Wn/eps;                           // [V/cm]
    double Ldn = sqrt(eps*Vt/(q*Nd));
    double Ldp = sqrt(eps*Vt/(q*Na));
    double Ldi = sqrt(eps*Vt/(q*ni));
    //printf("Vbi = %f\tW = %e\tWn = %e\tWp = %e\tWone = %e\tE_p = %e\tLdn = %e\tLdp = %e\tLdi = %e\n",Vbi,W,Wn,Wp,Wone,E_p,Ldn,Ldp,Ldi);


    /* Calculate relevant parameters in an input file
    % Write to a file
    save input_params.txt Na Nd Vbi W Wn Wp E_p Ldn Ldp
    %Material_Constants    %Define some material constants
    % Setting the size of the simulation domain based
    % on the analytical results for the width of the depletion regions for a simple pn-diode %*/
    double x_max = 0.0;
    if(x_max < Wn) x_max = Wn;
    if(x_max < Wp) x_max = Wp;
    x_max = 20*x_max;


    /* Setting the grid size based on the extrinsic Debye lengths */
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
    rnd0 = (Rnode*)calloc(n_max,sizeof(Rnode));
    rnd1 = rnd0 + n_max -1;
    Jnp0 = (Cur *)calloc(n_max, sizeof(Cur));
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


    /*dumping out the results*/
    FILE *ofp = fopen("/home/reza/python/els_plot.txt", "w");
    fprintf(ofp, "x [um]\tPotential [V]\tElectric Field1 [V/cm]\tElectric Field2 [V/cm]\tElectron Densities [1/cm^3]\tHole Densities [1/cm^3]\tTotal Charge Density [C/cm^3]\tConduction Band Energy (eV)\n");
    fclose(ofp);
    ofp = fopen("/home/reza/python/els_plot.txt", "a");
    for (rnd = rnd0; rnd <= rnd1; rnd++) {
        fprintf(ofp,"%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",rnd->xx1,Vt*(rnd->fi),rnd->el_field1,rnd->el_field2,rnd->nf,rnd->pf,rnd->ro,rnd->Ec);
    }
    fclose(ofp);
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

    double MU1N_CAUG   = 55.24;         //% cm2/(V.s)
    double MU2N_CAUG   = 1429.23;       //% cm2/(V.s)
    double ALPHAN_CAUG = 0.0;           //% unitless
    double BETAN_CAUG  = -2.3;          //% unitless
    double GAMMAN_CAUG = -3.8;          //% unitless
    double DELTAN_CAUG = 0.73;          //% unitless
    double NCRITN_CAUG = 1.072E17;   //% cm-3

    double MU1P_CAUG   = 49.7;          //% cm2/(V.s)
    double MU2P_CAUG   = 479.37;        //% cm2/(V.s)
    double ALPHAP_CAUG = 0.0;           //% unitless
    double BETAP_CAUG  = -2.2;          //% unitless
    double GAMMAP_CAUG = 13.7;          //% unitless
    double DELTAP_CAUG = 0.70;          //% unitless
    double NCRITP_CAUG = 1.606E17;   //% cm-3
    double BETAN = 2.0;
    double BETAP = 1.0;
    /*% %     mun0 = ( MU1N_CAUG*((TL/300)^ALPHAN_CAUG) ) ...
      % %       + (( (MU2N_CAUG*((TL/300)^BETAN_CAUG)) - (MU1N_CAUG*((TL/300)^ALPHAN_CAUG)) ) ...
      % %            / ( 1 + ((TL/300)^GAMMAN_CAUG) * ((N/NCRITN_CAUG)^DELTAN_CAUG) ))
      % %
      % %
      % %     mup0 = ( MU1P_CAUG*((TL/300)^ALPHAP_CAUG) ) ...
      % %       + (( (MU2P_CAUG*((TL/300)^BETAP_CAUG)) - (MU1P_CAUG*((TL/300)^ALPHAP_CAUG)) ) ...
      % %            / ( 1 + ((TL/300)^GAMMAP_CAUG) * ((N/NCRITP_CAUG)^DELTAP_CAUG) ))*/

    double VSATN = (2.4E7) / (1 + 0.8*exp(TL/600));  // Saturation Velocity of Electrons
    double VSATP = VSATN;                               // Saturation Velocity of Holes

    //printf("BETAN = %f\tBETAP = %f\tVSATN = %e\tVSATP = %e\n", BETAN, BETAP, VSATN, VSATP);


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
    IV0 = (Vpro **)calloc(74, sizeof(Vpro *));
    for (IV = IV0; IV < IV0 + 74; IV++) *IV = (Vpro *)calloc(5, sizeof(Vpro));
    //printf("%d ", (IV0[0] + 0)->Jelec);

    for (VA = 0.0, IV = IV0; VA < Vf; VA+=0.33*Vt, IV++) {                //% Start VA increment loop

        printf("VA = %f\n",VA);
        vindex = vindex + 1;
        (*IV)->Vplot = VA;
        //Vplot[vindex] = VA;
        //printf("total_steps = %d,  n_max = %d,  vindex = %d\n", Total_Steps, n_max, vindex);

        rnd0->fi = rnd0->fi + VA;            //% Apply potential to Anode (1st node)
        //%fi(1)

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

        //printf("fi0 = %e\tfi1 = %e\n", rnd0->fi, rnd1->fi);

    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 3. Start the Poisson equation solver loop to calculate the        %%
    %%    potential for each Anode voltage increase                      %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
        while(flag_conv2 == 0) {     //Start Poisson's eqn
            k_itern += 1;
            //printf("k_itern = %d\t", k_itern);


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
                //if (rnd - rnd0 < 10) printf("Ef = %e\tfi = %e\tfi+1 = %e\n", rnd->Ef, rnd->fi, (rnd + 1)->fi);
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
            //%fnp(1) = (Ldi*Ldi*dx2/Vt) * ( p(1)*n(1) - 1 ) / ( TAUP0*(n(1) + 1 ) + TAUN0*(p(1) + 1 ) );
            rnd0->fn = rnd0->n;
            rnd0->fp = rnd0->p;


            rnd1->an = 0;          //%Co-ef for electron at Cathode
            rnd1->bn = 1;          //%Co-ef for electron at Cathode
            rnd1->cn = 0;          //%Co-ef for electron at Cathode
            rnd1->ap = 0;          //%Co-ef for hole     at Cathode
            rnd1->bp = 1;          //%Co-ef for hole     at Cathode
            rnd1->cp = 0;          //%Co-ef for hole     at Cathode
            //%fnp(n_max) = (Ldi*Ldi*dx2/Vt) * ( p(n_max)*n(n_max) - 1 ) / ( TAUP0*(n(n_max) + 1) + TAUN0*(p(n_max) + 1) );
            rnd1->fn = rnd1->n;
            rnd1->fp = rnd1->p;

            //printf("n0 = %e\tp0 = %e\tn1 = %e\tp1 = %e\n", rnd0->n, rnd0->p, rnd1->n, rnd1->p);


        //%(B) Define the elements of the coefficient matrix for the internal nodes and
        //%    initialize the forcing function




            for (rnd = rnd0 + 1; rnd <= rnd1 - 1; rnd++) {
                munim1by2 = ((rnd - 1)->mun + rnd->mun)/2;
                munip1by2 = (rnd->mun + (rnd + 1)->mun)/2;
                mupim1by2 = ((rnd - 1)->mup + rnd->mup)/2;
                mupip1by2 = (rnd->mup + (rnd + 1)->mup)/2;
                //if (rnd1 - rnd < 10) printf("%e\t%e\t%e\t%e\n",munim1by2,munip1by2,mupim1by2,mupip1by2);

            //%% Co-efficients for HOLE Continuity eqn
                rnd->cp = mupip1by2 * BER(rnd->fi - (rnd + 1)->fi);
                rnd->ap = mupim1by2 * BER(rnd->fi - (rnd - 1)->fi);
                rnd->bp = -( mupim1by2 * BER((rnd - 1)->fi - rnd->fi) + mupip1by2 * BER((rnd + 1)->fi - rnd->fi));
            //%% Co-efficients for ELECTRON Continuity eqn
                rnd->cn = munip1by2 * BER((rnd + 1)->fi - rnd->fi);
                rnd->an = munim1by2 * BER((rnd - 1)->fi - rnd->fi);
                rnd->bn = -( munim1by2 * BER(rnd->fi - (rnd - 1)->fi) + munip1by2 * BER(rnd->fi - (rnd + 1)->fi));
            //%% Forcing Function for ELECTRON and HOLE Continuity eqns
                rnd->fn = (Ldi*Ldi*dx2/Vt) * ( (rnd->p)*(rnd->n) - 1 ) / ( TAUP0*(rnd->n + 1) + TAUN0*(rnd->p + 1));
                rnd->fp = (Ldi*Ldi*dx2/Vt) * ( (rnd->p)*(rnd->n) - 1 ) / ( TAUP0*(rnd->n + 1) + TAUN0*(rnd->p + 1));
                //if (rnd1 - rnd < 200) printf("i = %d\tfi = %e\tBER = %e\n",rnd1-rnd,rnd->fi,BER(rnd->fi - (rnd + 1)->fi));
            }


        /*%(C)  Start the iterative procedure for the solution of the linearized Continuity
        %     equation for "ELECTRONS" using LU decomposition method:*/


            rnd0->alphan = rnd0->bn;
            //printf("alphan0 = %e\n", rnd0->alphan);
            for (rnd = rnd0 + 1; rnd <= rnd1; rnd++) {
                rnd->betan = rnd->an/(rnd - 1)->alphan;
                rnd->alphan = rnd->bn - rnd->betan*(rnd - 1)->cn;
                //if (rnd1 - rnd < 10)
                    //printf("i = %d\tbetan(i) = %e\talphan(i-1) = %e\talphan(i) = %e\n", rnd1-rnd, rnd->betan, (rnd-1)->alphan, rnd->alphan);
                    //printf("i = %d\tan(i) = %e\tbn(i) = %e\tcn(i-1) = %e\n",rnd1-rnd,rnd->an,rnd->bn,(rnd-1)->cn);
            }

        //% Solution of Lv = f %

            rnd0->vn = rnd0->fn;
            for (rnd = rnd0 + 1; rnd <= rnd1; rnd++) {
                rnd->vn = rnd->fn - rnd->betan*(rnd - 1)->vn;
                //if (rnd1 - rnd < 10) printf("i = %d\tfn = %e\tbetan = %e\n", rnd1-rnd, rnd->fn, rnd->betan);
            }

        //% Solution of U*fi = v %
            tempn = rnd1->vn/rnd1->alphan;
            //printf("vn1 = %e\talphan1 = %e\n", rnd1->vn, rnd1->alphan);
            //%deltan(n_max) = tempn - n(n_max);
            rnd1->n=tempn;
            for (rnd = rnd1 - 1; rnd >= rnd0; rnd--) {       //%delta%
                tempn = (rnd->vn - rnd->cn*(rnd + 1)->n)/rnd->alphan;
              //%  deltan(i) = tempn - n(i);
                rnd->n = tempn;
                //if (rnd1-1-rnd < 10) printf("i = %d\tn = %e\n",rnd1-1-rnd,rnd->n);
            }

       //%%%%%%%%%%%%%%%%%%%%%%% END of ELECTRON Continuty Solver %%%%%%%%%%%


       // %(D)  Start the iterative procedure for the solution of the linearized Continuity
        //%     equation for "HOLES" using LU decomposition method:

            rnd0->alphap = rnd0->bp;
            for (rnd = rnd0 + 1; rnd <= rnd1; rnd++) {
                rnd->betap = rnd->ap/(rnd-1)->alphap;
                rnd->alphap = rnd->bp - rnd->betap*(rnd - 1)->cp;
                //if (rnd - rnd0 < 10) printf("betap = %e\talphap = %e\n", rnd->betap, rnd->alphap);
            }

        //% Solution of Lv = f %

            rnd0->vp = rnd0->fp;
            for (rnd = rnd0 + 1; rnd <= rnd1; rnd++) {
                rnd->vp = rnd->fp - rnd->betap*(rnd - 1)->vp;
                //if (rnd - rnd0 < 10) printf("fp = %e\tbetap = %e\n", rnd->fp, rnd->betap);
            }

        //% Solution of U*fi = v %

            tempp = rnd1->vp/rnd1->alphap;
            //%deltap(n_max) = tempp - p(n_max);
            rnd1->p=tempp;
            for (rnd = rnd1 - 1; rnd >= rnd0; rnd--) {       //%delta%
                tempp = (rnd->vp - rnd->cp*(rnd + 1)->p)/rnd->alphap;
             //%   deltap(i) = tempp - p(i);
                rnd->p = tempp;
                //if (rnd1 - rnd < 10) printf("vp = %e\tcp = %e\talphap = %e\n", rnd->vp, rnd->cp, rnd->alpha);
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
            //if (rnd - rnd0 < 10) printf("n = %e\tp = %e\tfi = %e\n", rnd->n, rnd->p, rnd->fi);
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
            //if (rnd - rnd0 < 10)
                //printf("f = %e\tbeta = %e\n", rnd->f, rnd->beta);
        }

        //% Solution of U*fi = v %

        temp = rnd1->v/rnd1->alpha;
        rnd1->delta = temp - rnd1->fi;
        rnd1->fi=temp;
        for (rnd = rnd1 - 1; rnd >= rnd0; rnd--) {      //%delta%
            temp = (rnd->v - rnd->c*(rnd + 1)->fi)/rnd->alpha;
            rnd->delta = temp - rnd->fi;
            //if (rnd1 - 1 - rnd < 10)
                //printf("i = %d\tv = %e\tc = %e\talpha = %e\tfi = %e\n", rnd1 - 1 - rnd, rnd->v, rnd->c, rnd->alpha, (rnd+1)->fi);
            rnd->fi = temp;
        }

        delta_max = 0.0;
        for (rnd = rnd0; rnd <= rnd1; rnd++) {
            //if (rnd - rnd0 < 10) printf("delta(%d) = %e\n", rnd - rnd0 +1, (rnd->delta));
            xx = abs1(rnd->delta);
            if(xx > delta_max) delta_max=xx;
        }

        printf("delta_max = %e\n", delta_max);
        if(delta_max < delta_acc) flag_conv2 = 1;
        }   // end of while loop

    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%                        CALCULATE CURRENT                             %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

        ofp = fopen("Current.txt", "a");
        for (rnd = rnd0 + 1, Jnp = Jnp0 + 1; rnd <= rnd1 -1; rnd++, Jnp++) {
            Jnp0->Jnim1by2 = (q*(rnd->mun)*Vt/(dx*Ldi)) * ni*( (rnd->n)*BER(rnd->fi - (rnd - 1)->fi) - ((rnd - 1)->n)*BER((rnd - 1)->fi - rnd->fi) );
            Jnp->Jnim1by2 = (q*(rnd->mun)*Vt/(dx*Ldi)) * ni*( ((rnd + 1)->n)*BER((rnd + 1)->fi - rnd->fi) - (rnd->n)*BER(rnd->fi - (rnd + 1)->fi) );
            Jnp->Jelec = (Jnp0->Jnim1by2 + Jnp->Jnim1by2)/2;

            Jnp->Jpim1by2 = (q*(rnd->mup)*Vt/(dx*Ldi)) * ni*( (rnd->p)*BER((rnd - 1)->fi - rnd->fi) - ((rnd - 1)->p)*BER(rnd->fi - (rnd - 1)->fi) );
            Jnp->Jpip1by2 = (q*(rnd->mup)*Vt/(dx*Ldi)) * ni*( ((rnd + 1)->p)*BER(rnd->fi - (rnd + 1)->fi) - (rnd->p)*BER((rnd + 1)->fi - rnd->fi) );
            Jnp->Jhole = (Jnp->Jpim1by2 + Jnp->Jpip1by2)/2;

            Jnp->Jtotal = Jnp->Jelec + Jnp->Jhole;

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

            fprintf(ofp,"%e\t%e\t%e",Jnp->Jelec,Jnp->Jhole,Jnp->Jtotal);
            if ((rnd - rnd0)%3 == 0) fprintf(ofp,"\n");
            else fprintf(ofp,"\t");
        }
        fclose(ofp);
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

    ofp = fopen("quasi.txt", "w");
    fprintf(ofp,"x [um]\tEc\tEv\tEi\tEfn\tEfp\n");
    for (rnd = rnd0; rnd <= rnd1; rnd++) fprintf(ofp, "%e\t%e\t%e\t%e\t%e\t%e\n",rnd->xx1,rnd->Ec,rnd->Ev,rnd->Ei,rnd->Efn,rnd->Efp);
    fclose(ofp);

    ofp = fopen("pote.txt", "w");
    fprintf(ofp,"x [um]\tPotential [eV]\n");
    for (rnd = rnd0; rnd <= rnd1; rnd++) fprintf(ofp, "%e\t%e\n", rnd->xx1, Vt*(rnd->fi));
    fclose(ofp);

    ofp = fopen("elec.txt", "w");
    fprintf(ofp,"x [um]\tElectric Field [V/cm]\tElectric Field [V/cm]\n");
    for (rnd = rnd0; rnd <= rnd1; rnd++) fprintf(ofp, "%e\t%e\t%e\n", rnd->xx1, rnd->el_field1, rnd->el_field2);
    fclose(ofp);

    ofp = fopen("car_dens.txt", "w");
    fprintf(ofp,"x [um]\tElectron Densities [1/cm^3]\tHole Densities [1/cm^3]\n");
    for (rnd = rnd0; rnd <= rnd1; rnd++) fprintf(ofp, "%e\t%e\t%e\n", rnd->xx1, rnd->nf, rnd->pf);
    fclose(ofp);

    ofp = fopen("charge.txt", "w");
    fprintf(ofp,"x [um]\tTotal Charge Density [C/cm^3]\n");
    for (rnd = rnd0; rnd <= rnd1; rnd++) fprintf(ofp, "%e\t%e\n", rnd->xx1, q*rnd->ro);
    fclose(ofp);

    ofp = fopen("curr_tot_ele_hol.txt", "w");
    fprintf(ofp,"Vplot\tTotal Current Density [Amp/cm^2]\tElectron Current Density [Amp/cm^2]\tHole Current Density [Amp/cm^2]\n");
    for (IV = IV0; IV - IV0 < 73; IV++) fprintf(ofp, "%f\t%e\t%e\t%e\n", (*IV)->Vplot, (*IV + 3)->Jtotal, (*IV + 3)->Jelec, (*IV + 3)->Jhole);
    fclose(ofp);

    ofp = fopen("curr_pos.txt", "w");
    fprintf(ofp,"x [um]\tTotal Current Density [Amp/cm^2]\n");
    for (rnd = rnd0, Jnp = Jnp0; rnd <= rnd1; rnd++, Jnp++) fprintf(ofp, "%e\t%e\n", rnd->xx1, Jnp->Jtotal);
    fclose(ofp);

    free(rnd0);
    free(Jnp0);
    free(IV0);
    return 0;
}
