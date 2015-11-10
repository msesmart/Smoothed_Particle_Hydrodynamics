
#define npar 10000   //Maximum Number of particles
#define nb_max 802  //Maximum Number of boundary particles.
#define nct_max 3000  //Maximum Number of Cells or Meshes
#define nplink_max 20  //

double PI=3.1415926535898;
double g=9.81;

double visc_dt=0.0; FILE *fp;
time_t time_begin, time_end;

double tdetail,grab_P,grab_E,grx,gry,grz; int ngrab,ngrabdet;

int ipoute,itime; char supp[10]="00000.dat", detsupp[10]="00000.dat";

int i_kernel,i_algorithm,i_densityFilter,i_viscos,iBC,i_periodicOBs[3],lattice,i_EoS;
int np,nb,nbf,ivar_dt,i_restartRun,i_kernelcorrection,iRiemannSolver,i_geometry;
int iTVD,i_vort,ndt_VerletPerform,ndt_FilterPerform,ndt_DBCPerform,nstep_verlet,nstep_DBC;
int nbfp1,nbp1,i_correct_pb,ncx,ncz,nct,ncz_ini,nct_ini,nbfm,index_tensile;

double h_SWL,B,Gamma,expont,coef,coef2,eps,rho0,coefficient,Co,viscos_val,visc_wall,vlx,vly,vlz,dx,dy,dz;
double h,two_h,four_h,six_h,h_over_100,h_over_1000,h2,fourh2,Time,dt,tmax,out,trec_ini,dtrec_det,t_sta_det,t_end_det,vnorm_mass;
double cs0,ddt_c,ddt_p,CFL_number,TE0,beta_lim,beta_coef2,beta_coef,one_over_beta_coef2,one_over_beta_coef;
double two_alpha_h,rNum_h,one_over_rNum_h,one_over_r_ij_dot_n,adh,eta,eta2,xmax_container,xmin_container,zmax_container,zmin_container;
double xtrash,ztrash,zcontrol,xmin,xmax,zmin,zmax,one_over_h,one_over_2h,a2,gauss,gauss2,f1,deltap2,Wdeltap,od_Wdeltap,Wap,Wab;
double frx,frz,frxi,frxj,frzi,frzj;

double xp[npar],zp[npar],up[npar],wp[npar],p[npar],rhop[npar],pm[npar],cs[npar],pVol[npar],sum_wab[npar],rhop_sum[npar];
double xcor[npar],zcor[npar],udot[npar],wdot[npar],xdot[npar],zdot[npar],rdot[npar];
double ax[npar],az[npar],pr[npar],ar[npar],ux[npar],wx[npar];
double xo[npar],zo[npar],uo[npar],wo[npar],po[npar],rhoo[npar];
double TEo[npar],TEp[npar],aTE[npar],TEdot[npar];
int iBP_Pointer_Info[npar][4],BP_xz_Data[npar][2],iflag[npar];

//BP Normals
int nstart,nstart_minus1; double xnb[nb_max],znb[nb_max],xtb[nb_max],ztb[nb_max],deltaptb[nb_max][2];
// Cells
int nc[nct_max][2],ibox[nct_max][2][nplink_max],keep_nc[nct_max];

