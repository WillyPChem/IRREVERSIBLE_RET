#include <fftw3.h>
#include<math.h>
#include<stdio.h>
#include<malloc.h>
#include<complex.h>

#define REAL 0
#define IMAG 1
//void Commutator (int dim, double H[dim][dim], double D[dim][dim], double P[dim][dim]);
void RK3(int Nlevel, double time, double *bas, double *E, double *Hint, double *Mu, double *Dis, double complex *D, double dt);
void RK3_MG(int Nlevel, double time, double *bas, double *E, double *Hint, double *Mu, double *Dis, double complex *D, double dt);
void Liouville(int dim, double complex *H, double complex *D, double complex *P);
void AntiCommutator(int dim, double *H, double complex *D, double complex *P);
void PrintComplexMatrix(int dim, double complex *M);
void PrintRealMatrix(int dim, double *M);
void FormDM(int dim, double complex *C, double complex *Ct, double complex *DM);
void L_Diss_MG(int Nlevel, double gamma, double complex *D, double *bas, double complex *P, int ex, int dex);
void L_Diss(int Nlevel, double *gamma, double complex *D, double *bas, double complex *P);
void Fourier (double complex *dm, int n, double dt, double complex *ftout, double *freqvec);
double complex TrMuD(int Nlevel, double *Mu, double complex *D);
double E_Field(double time);
void FillDFTArray(int fftw_iter, double real_val, double imag_val, fftw_complex* ar);
double D_Error(int dim, double complex *D);
// Function prototype for H_interaction
void H_interaction(int dim, double *Hint, double *mu, double dpm, double *R);
void DipoleAcceleration(int dim, double dt, fftw_complex* dp, fftw_complex* dpa);

// NOTE!!!  You need three global variables for the rates associated with 
// The Lindblad operators - gamma, beta, and alpha should be defined here according
// to their values in the journal of physical chemistry letters paper
//int Nlevel = 3;
//int dim = Nlevel*Nlevel;

double pi = 4*atan(1.0);
double wn_to_au = 4.55633528e-6; 
double mu_au_to_si = 8.47835326e-30; // 1 a.u. = 8.47835326e-30 C m
double E_au_to_si = 5.14220652e11;  // 1 a.u. = 5.14220652e11 V/m
double omega_au = 4.134147e+16;;
int main() {

  //Nanoparticles Variables here
  int numTime = 20000000;
  int zeropad = 2000000;
  double *E, *Mu, *Dis, *bas, *Hint;
  double complex *H, *D, *P;
  double dt = 0.001;
  int Nlevel, dim;
  
  // NP levels can be variable in principle
  //printf("  How many states are in your NP system? \n");
  //scanf("%i",&Nlevel);
  Nlevel = 2;
  dim = Nlevel*Nlevel;
  // MG variables here!
  int NlevelMG = 3;
  int dimMG = NlevelMG*NlevelMG;
  double *EMG, *MuMG, *MuZERO, *DisMG, *basMG, *HintMG;
  double complex *HMG, *DMG, *PMG;

  int dft_dim = numTime+zeropad;
  // FFTW variables here -> inputs to fft
  fftw_complex *dipole;
  fftw_complex *dipoleA;
  fftw_complex *efield;
  fftw_complex *dipoleMG;
  fftw_complex *dipoleMGA;
  fftw_complex *nps;
  fftw_complex *mgs;
  fftw_complex *efs;
  fftw_complex *emiss;
  fftw_complex *emissMG;

  // Allocate memory for FFTW arrays
/*  dipole = (fftw_complex*)malloc(dft_dim*sizeof(fftw_complex));
  dipoleA= (fftw_complex*)malloc(dft_dim*sizeof(fftw_complex));
  efield = (fftw_complex*)malloc(dft_dim*sizeof(fftw_complex));
  dipoleMG = (fftw_complex*)malloc(dft_dim*sizeof(fftw_complex));
  dipoleMGA= (fftw_complex*)malloc(dft_dim*sizeof(fftw_complex));
  nps = (fftw_complex*)malloc(dft_dim*sizeof(fftw_complex));
  mgs = (fftw_complex*)malloc(dft_dim*sizeof(fftw_complex));
  efs = (fftw_complex*)malloc(dft_dim*sizeof(fftw_complex));
  emiss = (fftw_complex*)malloc(dft_dim*sizeof(fftw_complex));
  emissMG = (fftw_complex*)malloc(dft_dim*sizeof(fftw_complex));

  printf("  just declared bunch of arrays! \n");
  fflush(stdout);
  fftw_plan npp = fftw_plan_dft_1d(dft_dim,
                                      dipole,
                                      nps,
                                      FFTW_BACKWARD,
                                      FFTW_ESTIMATE);

  fftw_plan mgp = fftw_plan_dft_1d(dft_dim,
                                      dipoleMG,
                                      mgs,
                                      FFTW_BACKWARD,
                                      FFTW_ESTIMATE);


  fftw_plan efp = fftw_plan_dft_1d(dft_dim,
                                      efield,
                                      efs,
                                      FFTW_BACKWARD,
                                      FFTW_ESTIMATE);

  fftw_plan emp = fftw_plan_dft_1d(dft_dim, 
 					dipoleA,
					emiss,
					FFTW_BACKWARD,
					FFTW_ESTIMATE);

  fftw_plan emgp = fftw_plan_dft_1d(dft_dim,
					dipoleMGA,
					emissMG,
					FFTW_BACKWARD,
					FFTW_ESTIMATE); 

  printf("  STILL GOOD AFTER FFTW PLANS\n");
  fflush(stdout);
*/
  // Allocate memory for all other arrays
  // NP
  H = (double complex*)malloc(dim*sizeof(double complex));
  D = (double complex *)malloc(dim*sizeof(double complex));
  P = (double complex *)malloc(dim*sizeof(double complex));
  E  = (double *)malloc(dim*sizeof(double));
  Mu = (double *)malloc(dim*sizeof(double));
  Dis = (double *)malloc(dim*sizeof(double));
  bas = (double *)malloc(dim*sizeof(double));
  Hint = (double *)malloc(dim*sizeof(double));
  // MG
  HMG = (double complex*)malloc(dimMG*sizeof(double complex));
  DMG = (double complex *)malloc(dimMG*sizeof(double complex));
  PMG = (double complex *)malloc(dimMG*sizeof(double complex));
  EMG  = (double *)malloc(dimMG*sizeof(double));
  MuMG = (double *)malloc(dimMG*sizeof(double));
  MuZERO = (double *)malloc(dimMG*sizeof(double));
  DisMG = (double *)malloc(dimMG*sizeof(double));
  basMG = (double *)malloc(dimMG*sizeof(double));
  HintMG = (double *)malloc(dimMG*sizeof(double));

  // Variables for instantaneous quantities  
  double tr, trMG;
  double complex dipole_moment, dipole_momentMG;
  FILE *dfp, *dfpMG;
  FILE *popfp, *popMGfp, *ecumfp;

  // Separation vector
  double *r;
  r = (double *)malloc(3*sizeof(double));
  r[0] = 50.;
  r[1] = 0.;
  r[2] = 0.;



  FILE *Efp, *Mufp, *Disfp, *EfpMG, *MufpMG, *DisfpMG, *absfp, *emsfp;

  // INPUT FILES!
  Efp = fopen("MATRICES/Energy_Au.txt","r");
  Mufp = fopen("MATRICES/Dipole_Au.txt","r");
  Disfp = fopen("MATRICES/Dissipation_Au.txt","r");

  EfpMG = fopen("MATRICES/Energy3LDye.txt","r");
  MufpMG = fopen("MATRICES/Dipole3LDye.txt","r");
  DisfpMG = fopen("MATRICES/Dissipation3LDye.txt","r");


  //printf("  STILL GOOD AFTER MATRICES READS\n");
  //fflush(stdout);
  // Density matrix element D(i,j) is accessed as D[i*Nlevel+j];
  D[0] = 1. + 0.*I;
  DMG[0] = 1. + 0.*I;
  // NP
  for (int i=1; i<dim; i++){
    D[i] = 0. + 0.*I;
  }
  // MG
  for (int i=1; i<dimMG; i++){
    DMG[i] = 0. + 0.*I;
  }

  // BUILD DM BASIS - this comes into play in Lindblad operator
  // NP
  for (int i=0; i<Nlevel; i++) {
    for (int j=0; j<Nlevel; j++) {
      if (i==j){
        bas[i*Nlevel+j] = 1.0;
      }     
      else{
        bas[i*Nlevel+j] = 0.;
      }
    }
  }
  // MG
  for (int i=0; i<NlevelMG; i++) {
    for (int j=0; j<NlevelMG; j++) {
      if (i==j){
        basMG[i*NlevelMG+j] = 1.0;
      }
      else{
        basMG[i*NlevelMG+j] = 0.;
      }
    }
  }
  
  // Get parameters for NP and MG from files
  double val;
  // NP
  for (int i=0; i<dim; i++) {


       // Read from energy file and store to the energy matrix
       fscanf(Efp,"%lf",&val);
       E[i] = val;

       fscanf(Mufp,"%lf",&val);
       Mu[i] = val;

       fscanf(Disfp,"%lf",&val);
       Dis[i] = val;
  }
  // MG
  for (int i=0; i<dimMG; i++) {

       fscanf(EfpMG,"%lf",&val);
       EMG[i] = val;

       fscanf(MufpMG,"%lf",&val);
       MuMG[i] = val;
       MuZERO[i] = 0.;
       fscanf(DisfpMG,"%lf",&val);
       DisMG[i] = val;


  }
  double E_donor    = E[2*1+1];
  double E_acceptor = EMG[3*2+2];
  double g_acceptor = DisMG[3*2+2];
  EMG[3*2+2] = E_donor;

  // Loop over donor energy
  int lmax, halflmax;
  lmax = 500;
  halflmax = lmax/2;

  double dE = (4.*E_acceptor-E_acceptor/4.)/(lmax-1.);
  double dG = (5.*g_acceptor-g_acceptor/500.)/(lmax-1.);

  //for (int i=1; i<=lmax; i++) {
  for (int i=0; i<1; i++) {
  printf("\n");
  //E[2*1+1] = E_acceptor/4. + i*dE;
  E[2*1+1] = E_donor;
  for (int j=0; j<lmax; j++) {

  DisMG[3*1+1] = g_acceptor/500. + j*dG;
  DisMG[3*2+2] = g_acceptor/500. + j*dG;
  // Get initial dipole moments

  // Re-initialize density matrix
  D[0] = 1. + 0.*I;
  DMG[0] = 1. + 0.*I;
  // NP
  for (int i=1; i<dim; i++){
    D[i] = 0. + 0.*I;
  }
  // MG
  for (int i=1; i<dimMG; i++){
    DMG[i] = 0. + 0.*I;
  }

  dipole_moment = TrMuD(Nlevel, Mu, D)*mu_au_to_si;
  dipole_momentMG = TrMuD(NlevelMG, MuMG, DMG)*mu_au_to_si;

  //FillDFTArray(0, creal(dipole_moment), cimag(dipole_moment), dipole);
  //FillDFTArray(0, creal(dipole_momentMG), cimag(dipole_momentMG), dipoleMG);
  //FillDFTArray(0, 0., 0., efield);


  //void H_interaction(int dim, double *Hint, double *mu, double dpm, double R) 
  H_interaction(Nlevel, Hint, Mu, creal(dipole_momentMG), r); 
  H_interaction(NlevelMG, HintMG, MuMG, creal(dipole_moment), r);

  for (int i=1; i<numTime; i++) {

    // Update the NP 
    RK3(Nlevel, dt*i, bas, E, Hint, Mu, Dis, D, dt);

    // Calculate dipole moment of NP
    dipole_moment = TrMuD(Nlevel, Mu, D); 

    // Calculate interaction matrix between NP and MG: H_int^{np->mg}    
    H_interaction(NlevelMG, HintMG, MuMG, creal(dipole_moment), r);
    
    // Update the MG
    RK3_MG(NlevelMG, dt*i, basMG, EMG, HintMG, MuMG, DisMG, DMG, dt);
    //RK3(NlevelMG, dt*i, basMG, EMG, HintMG, MuZERO, DisMG, DMG, dt);
    
    // Calculate the dipole moment of MG
    dipole_momentMG = TrMuD(NlevelMG, MuMG, DMG);
    
    // Calculate interaction matrix between MG and NP: H_int^{mg->np} 
    H_interaction(Nlevel, Hint, Mu, creal(dipole_momentMG), r); 
   
  }
  // Energy transferred is proportional to E_3 * pop_2  
  double E_transfer = EMG[3*2+2]*creal(DMG[3*1+1])*27.211;
  double gam_D = Dis[2*1+1]*27.211;
  double gam_A = DisMG[3*1+1]*27.211;
  double E_D = E[2*1+1]*27.211;
  double E_A = EMG[3*2+2]*27.211;
  
  //printf("  E_D          E_A         g_D        g_A        E_tot(eV)\n");
  printf("  %8.6e    %8.6e   %8.6e  %8.6e  %8.6e\n",E_D, E_A, gam_D, gam_A, E_transfer);

  }
  }
  
  return 0;
}


void PrintRealMatrix(int dim, double *M) {

  printf("\n");
  for (int i=0; i<dim; i++) {

    for (int j=0; j<dim; j++) {

      printf(" %f ",M[i*dim+j]);

    }
    printf("\n");
  }

  printf("\n");
}

void PrintComplexMatrix(int dim, double complex *M) {
 
  printf("\n");
  for (int i=0; i<dim; i++) {

    for (int j=0; j<dim; j++) {

      printf(" (%12.10e,%12.10e) ",creal(M[i*dim+j]),cimag(M[i*dim+j]));

    }
    printf("\n");
  }
  printf("\n");
}

void RK3(int Nlevel, double time, double *bas, double *E, double *Hint, double *Mu, double *Dis, double complex *D, double dt) {

  int i, j;
  double complex *D_dot, *D2, *D3, *D_np1, *k1, *k2, *k3;
  double complex *H, *LD;  // Contribution to Ddot from Lindblad dissipation
  double *gamma;
  int dim = Nlevel*Nlevel; 
  double Efield;

  D_dot = (double complex *)malloc(dim*sizeof(double complex));
  D2    = (double complex *)malloc(dim*sizeof(double complex));
  D3    = (double complex *)malloc(dim*sizeof(double complex));
  D_np1 = (double complex *)malloc(dim*sizeof(double complex));
  k1    = (double complex *)malloc(dim*sizeof(double complex));
  k2    = (double complex *)malloc(dim*sizeof(double complex));
  k3    = (double complex *)malloc(dim*sizeof(double complex));
  H     = (double complex *)malloc(dim*sizeof(double complex));
  LD    = (double complex *)malloc(dim*sizeof(double complex));
  gamma = (double *)malloc(Nlevel*sizeof(double));
  
  // Must zero out all elements of these arrays
  for (i=0; i<dim; i++) {
    D_dot[i] = 0. + 0.*I;
    D2[i] = 0. + 0.*I;
    D3[i] = 0. + 0.*I;
    D_np1[i] = 0. + 0.*I;
    k1[i] = 0. + 0.*I;
    k2[i] = 0. + 0.*I;
    k3[i] = 0. + 0.*I;
    H[i] = 0. + 0.*I;
   
  }

  for (i=0; i<Nlevel; i++) {
    gamma[i] = Dis[i*Nlevel+i];
  }
  Efield = E_Field(time);

  // Compute full Hamiltonian at current time t 
  for (i=0; i<dim; i++) {


      H[i] = E[i] + Hint[i] - Efield*Mu[i]; // - I*Dis[i];

  } 

  //PrintComplexMatrix(Nlevel, H);

  // Get dPsi(n)/dt at initial time!
  // Two main changes needed to couple the molecule and nanoparticle:
  // (1) Liouville function needs to include H_interaction
  // (2) We need to use Liouville/L_Diss to update both the molecule and the nanoparticle density matrix
  Liouville(Nlevel, H, D, D_dot);
  L_Diss(Nlevel, gamma, D, bas, LD);
  //PrintComplexMatrix(Nlevel, D);
  //PrintComplexMatrix(Nlevel, D_dot);


  // Compute approximate wfn update with Euler step
  for (i=0; i<dim; i++) {
    k1[i] = dt*(D_dot[i]+LD[i]);
    D2[i] = D[i] + k1[i]/2.;
  }

  // Update Field!
  Efield = E_Field(time+dt/2.);

  // Compute full Hamiltonian at partially updated time t 
  for (i=0; i<dim; i++) {

      H[i] = E[i] + Hint[i] - Efield*Mu[i]; // - I*Dis[i];

  }

  //PrintComplexMatrix(Nlevel, H);
  // Get dPsi(n+k1/2)/dt
  Liouville (Nlevel, H, D2, D_dot);
  L_Diss(Nlevel, gamma, D2, bas, LD);
  
  // Compute approximate wfn update with Euler step
  for (i=0; i<dim; i++) {
    k2[i] = dt*(D_dot[i] + LD[i]);
    D3[i] = D[i] + k2[i]/2.;
  }

  // Get dPsi(n+k2/2)/dt
  Liouville (Nlevel, H, D3, D_dot);
  L_Diss(Nlevel, gamma, D3, bas, LD);

  // Compute approximate update with Euler step
  for (i=0; i<dim; i++) {
    k3[i] = dt*(D_dot[i] + LD[i]);
    D_np1[i] = D[i] + k1[i]/6. + 2.*k2[i]/3. + k3[i]/6.;
    D[i] = D_np1[i];
}


  free(D_dot);
  free(D2);
  free(D3);
  free(D_np1);
  free(k1);
  free(k2);
  free(k3);
  free(H);
  free(LD);
  free(gamma);
}


void RK3_MG(int Nlevel, double time, double *bas, double *E, double *Hint, double *Mu, double *Dis, double complex *D, double dt) {

  int i, j;
  double complex *D_dot, *D2, *D3, *D_np1, *k1, *k2, *k3;
  double complex *H, *LD_rad, *LD_nonrad;  // Contribution to Ddot from Lindblad dissipation
  double *gamma;
  int dim = Nlevel*Nlevel; 
  double Efield;

  D_dot = (double complex *)malloc(dim*sizeof(double complex));
  D2    = (double complex *)malloc(dim*sizeof(double complex));
  D3    = (double complex *)malloc(dim*sizeof(double complex));
  D_np1 = (double complex *)malloc(dim*sizeof(double complex));
  k1    = (double complex *)malloc(dim*sizeof(double complex));
  k2    = (double complex *)malloc(dim*sizeof(double complex));
  k3    = (double complex *)malloc(dim*sizeof(double complex));
  H     = (double complex *)malloc(dim*sizeof(double complex));
  LD_rad    = (double complex *)malloc(dim*sizeof(double complex));
  LD_nonrad = (double complex *)malloc(dim*sizeof(double complex));
  gamma = (double *)malloc(Nlevel*sizeof(double));
  
  // Must zero out all elements of these arrays
  for (i=0; i<dim; i++) {
    D_dot[i] = 0. + 0.*I;
    D2[i] = 0. + 0.*I;
    D3[i] = 0. + 0.*I;
    D_np1[i] = 0. + 0.*I;
    k1[i] = 0. + 0.*I;
    k2[i] = 0. + 0.*I;
    k3[i] = 0. + 0.*I;
    H[i] = 0. + 0.*I;
   
  }

  for (i=0; i<Nlevel; i++) {
    gamma[i] = Dis[i*Nlevel+i];
  }
  double g_val1 = gamma[1];
  double g_val2 = gamma[2];
  Efield = E_Field(time);

  // Compute full Hamiltonian at current time t 
  for (i=0; i<dim; i++) {


      H[i] = E[i] + Hint[i] - Efield*Mu[i]; // - I*Dis[i];

  } 

  //PrintComplexMatrix(Nlevel, H);

  // Get dPsi(n)/dt at initial time!
  // Two main changes needed to couple the molecule and nanoparticle:
  // (1) Liouville function needs to include H_interaction
  // (2) We need to use Liouville/L_Diss to update both the molecule and the nanoparticle density matrix
  Liouville(Nlevel, H, D, D_dot);
  L_Diss_MG(Nlevel, g_val1, D, bas, LD_rad, 2, 0);
  L_Diss_MG(Nlevel, g_val2, D, bas, LD_nonrad, 2, 1);
  //PrintComplexMatrix(Nlevel, D);
  //PrintComplexMatrix(Nlevel, D_dot);


  // Compute approximate wfn update with Euler step
  for (i=0; i<dim; i++) {
    k1[i] = dt*(D_dot[i]+LD_rad[i]+LD_nonrad[i]);
    D2[i] = D[i] + k1[i]/2.;
  }

  // Update Field!
  Efield = E_Field(time+dt/2.);

  // Compute full Hamiltonian at partially updated time t 
  for (i=0; i<dim; i++) {

      H[i] = E[i] + Hint[i] - Efield*Mu[i]; // - I*Dis[i];

  }

  //PrintComplexMatrix(Nlevel, H);
  // Get dPsi(n+k1/2)/dt
  Liouville (Nlevel, H, D2, D_dot);
  L_Diss_MG(Nlevel, g_val1, D2, bas, LD_rad, 2, 0);
  L_Diss_MG(Nlevel, g_val2, D2, bas, LD_nonrad, 2, 1);
  
  // Compute approximate wfn update with Euler step
  for (i=0; i<dim; i++) {
    k2[i] = dt*(D_dot[i] + LD_rad[i]+LD_nonrad[i]);
    D3[i] = D[i] + k2[i]/2.;
  }

  // Get dPsi(n+k2/2)/dt
  Liouville (Nlevel, H, D3, D_dot);
  L_Diss_MG(Nlevel, g_val1, D3, bas, LD_rad, 2, 0);
  L_Diss_MG(Nlevel, g_val2, D3, bas, LD_nonrad, 2, 1);

  // Compute approximate update with Euler step
  for (i=0; i<dim; i++) {
    k3[i] = dt*(D_dot[i] + LD_rad[i]+LD_nonrad[i]);
    D_np1[i] = D[i] + k1[i]/6. + 2.*k2[i]/3. + k3[i]/6.;
    D[i] = D_np1[i];
}


  free(D_dot);
  free(D2);
  free(D3);
  free(D_np1);
  free(k1);
  free(k2);
  free(k3);
  free(H);
  free(LD_rad);
  free(LD_nonrad);
  free(gamma);
}




void FormDM(int dim, double complex *C, double complex *Ct, double complex *DM) {

  for (int i=0; i<dim; i++) {

    for (int j=0; j<dim; j++) {

      DM[i*dim+j] = C[i]*Ct[j];

    }

  }

}

void Liouville(int dim, double complex *H, double complex *D, double complex *P) {


  // write code here!
  for (int i=0; i<dim; i++) {

    for (int j=0; j<dim; j++) {


      double complex sum2 = 0.+0.*I;
      double complex sum1 = 0.+0.*I;
      for (int k=0; k<dim; k++) {

        sum1 -= H[i*dim+k]*D[k*dim+j]*I;
        sum2 += D[i*dim+k]*H[k*dim+j]*I;
      }
      P[i*dim+j] = sum1 + sum2;
    }
  }
}



void AntiCommutator(int dim, double *H, double complex *D, double complex *P) {

// write code here!
for (int i=0; i<dim; i++) {

  for (int j=0; j<dim; j++) {


    double complex sum2 = 0.+0.*I;
    double complex sum1 = 0.+0.*I;
    for (int k=0; k<dim; k++) {

      sum1 += H[i*dim+k]*D[k*dim+j];
      sum2 += D[i*dim+k]*H[k*dim+j];
    }
    P[i*dim+j] = sum1 + sum2;
    //printf(" Pb[%i][%i] is %f %f\n",i,j,creal(sum1-sum2),cimag(sum1-sum2));
}
}
}


double E_Field(double time) {

  double Ef;
  double tau = 75.;

  //Ef = 0.01*sin(pi*time/tau)*sin(pi*time/tau)*exp(-0.005*time)*(sin(0.07423*time)+sin(0.1*time)+sin(0.5*time));
  if (time<tau) {

    Ef = 0.001*sin(time*pi/tau)*sin(time*pi/tau)*sin(0.07423*time);

  }
  else Ef = 0.;
    
  return Ef;


}

void L_Diss(int Nlevel, double *gamma, double complex *D, double *bas, double complex *P) {

  int i, j, k;
  double *temp_bas, *g_bas;
  double complex *temp_t1, *temp_t2, *LD;
  temp_bas = (double *)malloc(Nlevel*Nlevel*sizeof(double));
  g_bas    = (double *)malloc(Nlevel*Nlevel*sizeof(double));
  temp_t1  = (double complex *)malloc(Nlevel*Nlevel*sizeof(double complex));
  temp_t2  = (double complex *)malloc(Nlevel*Nlevel*sizeof(double complex));
  LD       = (double complex *)malloc(Nlevel*Nlevel*sizeof(double complex));
 
  double gk;
  // Form |g><g| matrix
  for (i=0; i<Nlevel; i++) {
    for (j=0; j<Nlevel; j++) {
      g_bas[i*Nlevel+j] = bas[0*Nlevel+i]*bas[j*Nlevel+0];
      LD[i*Nlevel+j] = 0. + 0.*I;
    }
  }
  //printf("  |g><g|  \n");
  //PrintRealMatrix(Nlevel, g_bas);

  for (k=1; k<Nlevel; k++) {

    gk = gamma[k];
    for (i=0; i<Nlevel; i++) {

      for (j=0; j<Nlevel; j++) {

        temp_bas[i*Nlevel+j] = bas[k*Nlevel+i]*bas[j*Nlevel+k];
        temp_t1[i*Nlevel+j] = 2*D[k*Nlevel+k]*g_bas[i*Nlevel+j];
      }
    }
   // printf("   |%i><%i| \n",k,k);
   // PrintRealMatrix(Nlevel, temp_bas);

    AntiCommutator(Nlevel, temp_bas, D, temp_t2);
    for (i=0; i<Nlevel; i++) {
      for (j=0; j<Nlevel; j++) {
        LD[i*Nlevel+j] += gk*temp_t1[i*Nlevel+j] - gk*temp_t2[i*Nlevel+j];
      }
    }
 }
 for (i=0; i<Nlevel; i++) {
   for (j=0; j<Nlevel; j++) {
     P[i*Nlevel+j] = LD[i*Nlevel+j];
   }
 }

 //PrintComplexMatrix(Nlevel, P);

 free(temp_bas);
 free(g_bas);
 free(temp_t1);
 free(temp_t2);
 free(LD);
}


void L_Diss_MG(int Nlevel, double gamma, double complex *D, double *bas, double complex *P, int ex, int dex) {

  int i, j, k;
  double *temp_bas, *g_bas;
  double complex *temp_t1, *temp_t2, *LD;
  temp_bas = (double *)malloc(Nlevel*Nlevel*sizeof(double));
  g_bas    = (double *)malloc(Nlevel*Nlevel*sizeof(double));
  temp_t1  = (double complex *)malloc(Nlevel*Nlevel*sizeof(double complex));
  temp_t2  = (double complex *)malloc(Nlevel*Nlevel*sizeof(double complex));
  LD       = (double complex *)malloc(Nlevel*Nlevel*sizeof(double complex));
 
  double gk;
  // Form |g><g| matrix
  for (i=0; i<Nlevel; i++) {
    for (j=0; j<Nlevel; j++) {
      g_bas[i*Nlevel+j] = bas[dex*Nlevel+i]*bas[j*Nlevel+dex];
      LD[i*Nlevel+j] = 0. + 0.*I;
    }
  }
  //printf("  |g><g|  \n");
  //PrintRealMatrix(Nlevel, g_bas);


    gk = gamma;
    for (i=0; i<Nlevel; i++) {

      for (j=0; j<Nlevel; j++) {

        temp_bas[i*Nlevel+j] = bas[ex*Nlevel+i]*bas[j*Nlevel+ex];
        temp_t1[i*Nlevel+j] = 2*D[ex*Nlevel+ex]*g_bas[i*Nlevel+j];
      }
    }
   // printf("   |%i><%i| \n",k,k);
   // PrintRealMatrix(Nlevel, temp_bas);

    AntiCommutator(Nlevel, temp_bas, D, temp_t2);
    for (i=0; i<Nlevel; i++) {
      for (j=0; j<Nlevel; j++) {
        LD[i*Nlevel+j] += gk*temp_t1[i*Nlevel+j] - gk*temp_t2[i*Nlevel+j];
      }
    }
 for (i=0; i<Nlevel; i++) {
   for (j=0; j<Nlevel; j++) {
     P[i*Nlevel+j] = LD[i*Nlevel+j];
   }
 }

 //PrintComplexMatrix(Nlevel, P);

 free(temp_bas);
 free(g_bas);
 free(temp_t1);
 free(temp_t2);
 free(LD);
}

void Fourier(double complex *dm, int n, double dt, double complex *ftout, double *freqvec){
  //FILE *fp;
  //fp = fopen("Absorption_SpectrumAu.txt","w");
  double wmin=0.5*0.07;
  double wmax=2*0.07;
  int maxk = 500;
  double dw = (wmax-wmin)/maxk;
  double time;
 
  for (int k = 0; k <=maxk; k++) {
    double sumreal = 0;
    double sumimag = 0;
    double  w = wmin+k*dw;

    for (int t = 0; t < n; t++){
      time = dt*t;
      double angle = time*w;
      sumreal += creal(cexp(-I*angle)*dm[t])*dt;
      sumimag += cimag(cexp(-I*angle)*dm[t])*dt;
      //sumreal += creal(dm[t]) * cos(angle) + cimag(dm[t]) * sin(angle);
      //sumimag  += creal(dm[t]) * sin(angle) + cimag(dm[t]) * cos(angle);
    }
    ftout[k] = sumreal + sumimag*I;
    // Energy in eV
    freqvec[k] = w*27.2114;
    //fprintf(fp," %12.10e  %12.10e\n",w*(27.2114),sumreal*sumreal+sumimag*sumimag);
  }
  //fclose(fp);
}

double complex TrMuD(int Nlevel, double *Mu, double complex *D) {
  double complex tr = 0. + 0.*I;
  for (int i=0; i<Nlevel; i++) {

    double complex sum = 0. + 0.*I;
    for (int k=0; k<Nlevel; k++) {

      sum += Mu[i*Nlevel+k]*D[k*Nlevel+i];

    }
    tr += sum;

  } 

  return tr;

}


void H_interaction(int dim, double *Hint, double *mu, double dpm, double *R) {
  
  int i; 
 // double *tmp1, *tmp2;
  double oer2, oer3;
  double scal_R = sqrt(R[0]*R[0] + R[1]*R[1] + R[2]*R[2]);
 
  oer2 = pow(scal_R,-2.);
  oer3 = pow(scal_R,-3.);
 
  // Write code between here!
 
 for (i=0; i<dim*dim; i++){

   // Very important!  Currently assuming z-polarized light, so only the <mu>_z and mu_z terms
   //                  are non-zero, hence we consider only R[2]*mu and <mu>*R[2] 
   Hint[i] = oer3*(dpm*mu[i] -3*R[2]*mu[i]*R[2]*dpm*oer2);
 } 

}

void FillDFTArray(int fftw_iter, double real_val, double imag_val, fftw_complex* ar) {

  ar[fftw_iter][REAL] = real_val;
  ar[fftw_iter][IMAG] = imag_val;

}

// Computes the Frobenius norm of the difference between the current density matrix and the initial density matrix
double D_Error(int dim, double complex *D) {

  double complex e;
  double complex one = 1. + 0.*I;
  double complex zero = 0. + 0.*I;
  double FN;

  e = 0.+0*I;
 
  e = (one - D[0])*conj(one-D[0]);

  for (int i=1; i<dim; i++) {

    e += (zero - D[i])*conj(zero-D[i]);  

  }

  FN = creal(csqrt(e));

  return FN;
}

void DipoleAcceleration(int dim, double dt, fftw_complex* dp, fftw_complex* dpa) {

  double complex st1, st2, st3, st4, st5, accel;

  // initialize stencil variables to zero
  st1 = 0.+0.*I;
  st2 = 0.+0.*I;
  st3 = 0.+0.*I;
  st4 = 0.+0.*I;
  st5 = 0.+0.*I;  
  for (int i=0; i<dim; i++) {

    st1 = st2;
    st2 = st3;
    st3 = st4;
    st5 = dp[i][REAL] + dp[i][IMAG]*I;
    accel = (-1./12.)*(st1+st5) + (4./3.)*(st2+st4)-(5./2.)*st3;
    accel /= (dt*dt);
    dpa[i][REAL] = creal(accel);
    dpa[i][IMAG] = cimag(accel);

    //printf("  %12.10e  %12.10e  %12.10e  %12.10e  %12.10e\n",dt*i,creal(st5),cimag(st5),creal(accel),cimag(accel));    

  }

}
