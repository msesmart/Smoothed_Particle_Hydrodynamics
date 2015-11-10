

void pos_veloc(int nn,double xpos,double ypos,double zpos,double uveloc,double vveloc,double wveloc)
{
	xp[nn]=xpos; zp[nn]=zpos;
	up[nn]=uveloc; wp[nn]=wveloc;
}

void pressure(int nn,double ZZmax,double dx,double dy,double dz,double expont,double g)
{
    rhop[nn]=rho0*pow(1.0+rho0*g*(ZZmax-zp[nn])/B,expont);
    p[nn]=B*(pow(rhop[nn]/rho0,Gamma)-1.0);
	pm[nn]=vnorm_mass*rhop[nn]*dx*dz;
}

void boundaries(int *nn,int N0,int N,int M,int L,double dx,double dy,double dz,double beta,double theta)
{
    int N_start,N_finish,L_start,L_finish,i,j,k,i_minus1,i_plus1; double mdx,mdy,mdz;
    //Lateral Walls: Left  x = x1; Left  wall (x = x1) starts at ',nn+1
    N_start=0; N_finish=2*N-1; L_start=0; L_finish=2*L-1; j=0;
    mdx=0.5*dx; mdy=0.5*dy; mdz=0.5*dz;

    if(iBC==1) //Repulsive BC - Repulsive Force Particles with Normals
    {
        for(i=N_start;i<=N_finish;i++)
        {
            *nn=*nn+1;
            pos_veloc(*nn,i*mdx,j*mdy,0*mdz,0.0,0.0,0.0);
            if(i==N_start){ i_minus1=*nn; i_plus1=*nn+1;}
            else if(i==N_finish){ i_minus1=*nn-1; i_plus1==*nn; }
            else { i_minus1=*nn-1; i_plus1=*nn+1; }
                // -- Normal information - neighbour data for Repulsive Boundary Particles (BPs)
            iBP_Pointer_Info[*nn][2]= i_minus1;      //i-1 neighbour
            iBP_Pointer_Info[*nn][3]= i_plus1;        //i+1 neighbour
        }
        for(k=L_start;k<=L_finish;k++)
        {
            *nn=*nn+1;
            pos_veloc(*nn,N_finish*mdx,j*mdy,k*mdz,0.0,0.0,0.0);
            if(k==L_start){ i_minus1=*nn; i_plus1=*nn+1; }
            else if(k==L_finish){ i_minus1=*nn-1; i_plus1=*nn; }
            else { i_minus1=*nn-1; i_plus1=*nn+1; }
            // -- Normal information - neighbour data for Repulsive Boundary Particles (BPs)
            iBP_Pointer_Info[*nn][2]= i_minus1;      //i-1 neighbour
            iBP_Pointer_Info[*nn][3]= i_plus1;        //i+1 neighbour
        }
        for(i=N_finish;i>=N_start;i--)
        {
            *nn=*nn+1;
            pos_veloc(*nn,i*mdx,j*mdy,L_finish*mdz,0.0,0.0,0.0);
            if(i==N_start){ i_minus1=*nn-1; i_plus1=*nn; }
            else if(i==N_finish){ i_minus1=*nn; i_plus1=*nn+1; }
            else { i_minus1=*nn-1; i_plus1=*nn+1; }
            iBP_Pointer_Info[*nn][2]= i_minus1;      //i-1 neighbour
            iBP_Pointer_Info[*nn][3]= i_plus1;        //i+1 neighbour
        }
        for(k=L_finish;k>=L_start;k--)
        {
            *nn=*nn+1;
            pos_veloc(*nn,0*mdx,j*mdy,k*mdz,0.0,0.0,0.0);
            if(k==L_start){ i_minus1=*nn-1; i_plus1=*nn; }
            else if(k==L_finish){ i_minus1=*nn; i_plus1=*nn+1; }
            else { i_minus1=*nn-1; i_plus1=*nn+1; }
            iBP_Pointer_Info[*nn][2]= i_minus1;      //i-1 neighbour
            iBP_Pointer_Info[*nn][3]= i_plus1;        //i+1 neighbour
        }
    }
}

void set(int *nn,double g,double expont,double ZZmin,double ZZmax,double dx,double dy,double dz)
{
    printf("Have not finished this function for Initial_Condition=1 !\n"); exit(0);
}

void fluid_particles(int *nn,int N,int M,int L,double dx,double dy,double dz,double expont,double g,int i_correct_pb,double XXmin,double XXmax,double YYmin,double YYmax,double ZZmin,double ZZmax,double beta,double theta)
{
	//Initial conditions'
	//1=Set of particles without grid; 2=Fill a square shaped region; 3=Round shaped drop  '
	int initial=2,N_end,N_ini,L_end,L_ini,i,k; double x1,y1,z1;

    //I_CORRECT_PB =1 B. pressure will  be corrected;0 B.pressure will not be corrected

	if(initial==1)set(nn,g,expont,ZZmin,ZZmax,dx,dy,dz);
	else if(initial==2)
    {
        N_end=(int)(XXmax/dx); N_ini=(int)(XXmin/dx); XXmin = N_ini*dx; XXmax = N_end*dx;
        L_ini=(int)(ZZmin/dz); ZZmin=L_ini*dz;
        if(ZZmax>h_SWL)ZZmax=h_SWL; L_end=(int)(ZZmax/dz); ZZmax=L_end*dz;

        for(i=N_ini;i<=N_end;i++)
        {
            for(k=L_ini;k<=L_end;k++)
            {
                x1=XXmin+(i-N_ini)*dx; y1=0.0; z1=ZZmin+(k-L_ini)*dz;
                *nn=*nn+1;
                pos_veloc(*nn,x1,y1,z1,0.0,0.0,0.0);
                pressure(*nn,ZZmax,dx,dy,dz,expont,g);
            }
        }
    }
	else { printf("is not a valid I. Condition \n"); exit(0); }
}

void p_boundaries(double dx,double dy,double dz)
{
    double dr; int nn; int nbfm=nb;//floating bodies, No use in here
	dr=dx*dz;
    for(nn=1;nn<=nb;nn++)
    {
        rhop[nn]=rho0; p[nn]=0.0;
        if(iBC==1)pm[nn]=rhop[nn]*dr;
        else if(iBC==2)
        {
            if(nn>nbfm)pm[nn]=rhop[nn]*dr;  //only one layer of particles for Floating Bodies faces
            else pm[nn]=0.5*rhop[nn]*dr;
        }
        else ;
    }
}

void box(double g,double expont,double dx,double dy,double dz)     // ... Generate a 2D BOX
{
    int N,M,L,nn; double beta_deg,beta,theta,XXmin,XXmax,YYmax,YYmin,ZZmin,ZZmax;
    beta_deg=0.0; beta=beta_deg*PI/180.0; theta=0.0;i_correct_pb=0;

    N=(int)(vlx/dx)+1; M=0; L=(int)(vlz/dz)+1; nn=0;	//Initial number of particles
    if(i_periodicOBs[0]!=1&&i_periodicOBs[2]!=1)boundaries(&nn,0,N,M,L,dx,dy,dz,beta,theta);
    //Number of B. Fixed Particles, B. particles can be fixed or move according to external dependence(F.E.gates, wavemakers).
    nbf=nn; nbfm=nn;/* No moving objects*/ nb=nn; /*Number of B. Particles*/
    XXmin=0.02; XXmax=1.0; YYmax=0.0;YYmin=0.0; ZZmin=0.02; ZZmax=2.0; //Cube containing particles'
    fluid_particles(&nn,N,M,L,dx,dy,dz,expont,g,i_correct_pb,XXmin,XXmax,YYmin,YYmax,ZZmin,ZZmax,beta,theta);
    np=nn; printf("nb=%d ,   np= %d \n",nb,np);
    p_boundaries(dx,dy,dz);
}
void Init_para_SPH(void)
{
    i_kernel=1; //1=Gaussian; 2=Quadratic; 3=Cubic Spline; 5=Quintic
    i_kernelcorrection=0; //Kernel correction:(0=none, 1=Kernel correct., 2=Kernel Gradient Correction.)
    i_algorithm=1; //1=Predictor1 Corrector; 2=Verlet; 3=Symplectic; 4=Beeman
    ndt_VerletPerform=0; //Number of time steps to apply the Eulerian Equa. with Verlet algorithm,if not Verlet,=0

    i_densityFilter=1;  //0=no filter; 1=Shepard; 2=MLS
    ndt_FilterPerform=30; // Number of time steps to apply the density filter

    i_viscos=2; //1=Artificial; 2=Laminar; 3=Laminar + SPS
    viscos_val=1.0E-6;    //Viscosity parameter,it corresponds to α(Monaghan and Koss 1999) if i_viscos=1
                    //and to ν (kinematical viscosity) if i_viscos =2 or 3.
    i_EoS=1;    //Equation of State:1=Tait equation(Weakly Compressible); 2=Ideal Gas;3=Poisson equation(ImCompressible)
    h_SWL=1.0;  //Still water level (m)
    coef=12.0;  //Coefficient of Speed of sound ( 10 , 40 )
    coef2=coef*coef;
    rho0=1000.0;  //Reference density (Default value 1000 kg/m3)
    Gamma=7.0; //Parameter in Equation of State (Monaghan and Koss, 1999) (Default value 7)
    expont=1.0/Gamma;
    TE0=0.0;   //Initial value for the thermal energy simulating an Ideal Gas
    B=coef2*rho0*g*(h_SWL)/Gamma; //Parameter in Equation of State (Monaghan and Koss, 1999).
    Co=sqrt((B*Gamma)/rho0); // the speed of sound;

    iBC=1;  // Boundary conditions. 1= Monaghan repulsive forces; 2= Dynamic boundaries.
    visc_wall=8.0E-4;     //Wall viscosity value for Repulsive Force BC
    ndt_DBCPerform=0; //Number of time steps to apply the Hughes & Graham(2010)correction for dynamic B.C.
    i_periodicOBs[0]=0; //Periodic Lateral boundaries in x direction? (1=yes)
    i_periodicOBs[1]=0; // Y direction
    i_periodicOBs[2]=0; // Z direction

    i_geometry=1; //1= BOX, 2= BEACH, 3= EXTERNALLY FIXED GEOMETRY
    lattice=1;  //Initial Fluid Particle Structure: 1=simple cubic; 2=BCC(body-centered cubic)
    vnorm_mass=1.0; // when lattice=1; =0.5 when lattice =2
    //Dimensions of container box, medium extent in X,Y,Z direction, vly=0 for 2D
    vlx=2.0; vly=0.0; vlz=2.0;
    //Initial interparticle spacing in x, y, z direction. dy=0 for 2D
    dx=0.02; dy=0.0; dz=0.02;

    box(g,expont,dx,dy,dz);

    eps=0.5;  //Epsilon parameter in XSPH approach (Default value 0.5)
    ivar_dt=0;   //Variable time step calculated when ivar_dt=1
    dt=0.0001;    //Initial time step. It is kept throughout the run when ivar_dt=0
    tmax=3.0;  //RUN duration (in seconds), End of the run(seconds)
    out=0.02;   // output every out seconds = Recording time step (in seconds).
            //The position, vel, den, p & mass of every particle is recorded in PART file every out seconds.
    trec_ini=0.0;      // Initial output = Initial recording time
    //Detailed recording time step & Start time in detailed recording & End time in detailed recording.
    dtrec_det=0.1; t_sta_det=0.0; t_end_det=3.0;
    CFL_number=0.2; //0.1-0.5,Constant to calculate the time step following CFL condition
    coefficient=0.92;  //Coeffficient to calculate the smoothing length (h) in terms of dx,dy,dz;
    h=coefficient*sqrt(dx*dx+dy*dy+dz*dz); //Smoothing length
    printf("h= %lf \n",h);
    //0=Start a new RUN; 1=Restart an old RUN; 2=New with CheckPointg; 3=Restart with CheckPointing
    i_restartRun=0;
    iRiemannSolver=0;    //Use of Riemann Solver: (0=None, 1=Conservative, 2=NonConservative)
    iTVD=0;      //Use TVD, slope limiter (beta_lim)? (1=yes)
    beta_lim=0.0;      //slope limiter using Riemann Solver
    i_vort=0;        //vorticity printing ? (1=yes)
    index_tensile=0;
}

void equation_of_state(double rho_EoS,double TE_EoS,double *press_EoS,double *cs_EoS,int i_EoS)
{
    double cs0_2,backgroundPressure;
    if(i_EoS==1)
    {   *press_EoS=B*(pow(rho_EoS/rho0,Gamma)-1.0); *cs_EoS=cs0*pow(rho_EoS/rho0,3);}
    else if(i_EoS==2)
    {   *press_EoS=rho_EoS*(Gamma-1.0)*TE_EoS; *cs_EoS=sqrt(Gamma*(*press_EoS/rho_EoS));}
    else if(i_EoS==3)
    {   *cs_EoS=cs0; cs0_2=pow(*cs_EoS,2); backgroundPressure=0.0; *press_EoS = backgroundPressure+cs0_2*(rho_EoS-rho0);}
    else {printf("ERROR: Invalid i_EoS \n");exit(0);}
}


void pout_BoundaryNormals(char *filename)
{
    FILE *fp; int i;
    fp=fopen(filename,"w");
	fprintf(fp," VARIABLES=\"X \",\"Z \",\"xnb \",\"znb \",\"xtb \",\"ztb \" \n");
	fprintf(fp," ZONE F=POINT \n");
	for(i=1;i<=nb;i++)
    fprintf(fp,"%lf %lf %lf %lf %lf %lf \n",xp[i],zp[i],xnb[i],znb[i],xtb[i],ztb[i]);
    fclose(fp);
}

void updateNormals_2D(int i_start,int i_stop)
{
    int i,iminus1_neighbour,iplus1_neighbour;
    double xp_i,zp_i,xp_BP_iminus1,zp_BP_iminus1,xp_BP_iplus1,zp_BP_iplus1;
    double xp_previous,zp_previous,xp_next,zp_next,dxn_forward,dzn_forward;
    double dxn_back,dzn_back,rbn_back,rbn_forward,dxn,dzn,rb;
    //-- Serial Version --
    for(i=i_start;i<=i_stop;i++)
    {
        xp_i=xp[i]; zp_i=zp[i];

        iminus1_neighbour=iBP_Pointer_Info[i][2];
        iplus1_neighbour =iBP_Pointer_Info[i][3];

        xp_BP_iminus1 = xp[iminus1_neighbour];
        zp_BP_iminus1 = zp[iminus1_neighbour];
        xp_BP_iplus1 = xp[iplus1_neighbour];
        zp_BP_iplus1 = zp[iplus1_neighbour];

        if(i==iminus1_neighbour)
        {
            xp_previous = xp_i-(xp_BP_iplus1-xp_i);
            zp_previous = zp_i-(zp_BP_iplus1-zp_i);
            xp_next = xp_BP_iplus1;
            zp_next = zp_BP_iplus1;
        }
        else if(i==iplus1_neighbour)
        {
            xp_previous = xp_BP_iminus1;
            zp_previous = zp_BP_iminus1;
            xp_next = xp_i+(xp_i-xp_BP_iminus1);
            zp_next = zp_i+(zp_i-zp_BP_iminus1);
        }
        else
        {
            xp_previous = xp_BP_iminus1;
            zp_previous = zp_BP_iminus1;
            xp_next = xp_BP_iplus1;
            zp_next = zp_BP_iplus1;
        }
        dxn_forward=xp_next-xp_i; dzn_forward=zp_next-zp_i;
        dxn_back=xp_i-xp_previous; dzn_back=zp_i-zp_previous;

        rbn_back=sqrt(dxn_back*dxn_back+dzn_back*dzn_back); //该边界粒子与后一个边界粒子的距离
        rbn_forward=sqrt(dxn_forward*dxn_forward+dzn_forward*dzn_forward);//该边界粒子与前一个边界粒子的距离
        dzn=dzn_forward+dzn_back;
        dxn=dxn_forward+dxn_back;
        rb=dxn*dxn+dzn*dzn;
        if(rb==0.0)
        {
            printf("xp_BP_iplus1 =%lf, i=%d, iplus1_neighbour=%d, iminus1_neighbour=%d \n",xp_BP_iplus1,i,iplus1_neighbour,iminus1_neighbour);
            printf("xp_i=%lf, zp_i=%lf \n",xp_i,zp_i);
            printf("xp_previous=%lf, zp_previous=%lf, xp_next=%lf, zp_next=%lf \n",xp_previous,zp_previous,xp_next,zp_next);
        }
        xnb[i]=-dzn/sqrt(rb);  //xnb,znb: Components of the unitary vector normal to the boundary at point i
        znb[i]=dxn/sqrt(rb);
        xtb[i]=znb[i];   //xtb,ztb: Components of the unitary vector tangential to the boundary at point
        ztb[i]=-xnb[i];

        if(fabs(rbn_back)>0.0)deltaptb[i][0]=rbn_back; else deltaptb[i][0]=h*1.0E-3;
        if(fabs(rbn_forward)>0.0)deltaptb[i][1]=rbn_forward; else deltaptb[i][1]=h*1.0E-3;
    }
    pout_BoundaryNormals("BoundaryNormals.dat");
}

void pout(char *filename)
{
    FILE *fp; int i;
    fp=fopen(filename,"w");
	fprintf(fp," VARIABLES=\"X \",\"Z  \",\"U \",\"W \",\"P \",\"den \",\"m \" \n");
	fprintf(fp," ZONE F=POINT \n");
	for(i=1;i<=np;i++)fprintf(fp,"%lf %lf %lf %lf %lf %lf %lf \n",xp[i],zp[i],up[i],wp[i],p[i],rhop[i],pm[i]);
    fclose(fp);
    printf("output %s \n",filename);
    time(&time_end); printf("Current Used Time : %d seconds\n",time_end-time_begin);
}


void getdata(void)
{
    FILE *fp; int i; double dp,twopi,a1,aa,b1,c1,c2,d1,e1,a24;
    double deltap,Awen,Bwen;

    Init_para_SPH();

    if(np>=npar){printf("Number of particles exceeds npar in common!\n"); exit(0);}
    if(nb>=nb_max){printf("Number of Boundary particles exceeds nb_max!\n"); exit(0);}

    nbfp1=nbf+1; nbp1=nb+1;
    if(i_EoS==1){cs0=sqrt(Gamma*B/rho0);/*i_gamma=(int)(gamma);*/}
    ddt_c=dt; ddt_p=dt; //To be used with variable time step & predictor corrector
    if(i_algorithm==2)nstep_verlet=0; //Initialisation of verlet counter

    if(lattice==1) vnorm_mass=1.0; else if(lattice==2)vnorm_mass=0.5; else ;

    if(iBC==1)
    {
        beta_coef2 =(coef/10.0)*(coef/10.0); beta_coef=coef/10.0;
        one_over_beta_coef2=1.0/beta_coef2; one_over_beta_coef=1.0/beta_coef;
    }
    else{ one_over_beta_coef2=0.0; one_over_beta_coef=0.0;}

    itime=0; ngrab=0; grab_E=0.0;
    for(i=1;i<=np;i++)
    {
        equation_of_state(rhop[i],TEp[i],&p[i],&cs[i],i_EoS);
        iflag[i]=1; //To detect particles inside the container
        pVol[i]=pm[i]/rhop[i];
    }

    if(iBC==1)
    {
        for(i=1;i<=nb;i++)
        {
            xnb[i]=0;znb[i]=0;
            iBP_Pointer_Info[i][0]=0;    //Absolute index BP
            iBP_Pointer_Info[i][1]=0;   //Rank of BP (default=0, reserved for MPI)
            BP_xz_Data[i][0]=0;
            BP_xz_Data[i][1]=0;
        }
        nstart=nbp1;
        //calculates the new normals of the moving boundary particles when Monaghan Boundary conditions are used.
        updateNormals_2D(1,nb);
        two_alpha_h = 2.0*1.3; rNum_h = 3.0*h;
        one_over_rNum_h = 1.0/rNum_h; one_over_r_ij_dot_n = 1.3/h;
    }
    else nstart=1;
    nstart_minus1=nstart-1;
    twopi=2.0*PI;

    //Different Kernel Smoothing Lengths
    two_h=2.0*h; four_h = 4.0*h; six_h=6.0*h; h_over_100=0.01*h; h_over_1000=0.001*h; h2=h*h; fourh2=4.0*h2;

    if(i_kernel==1)//2D Gaussian
    {
        gauss=1.0; gauss2=gauss*gauss; a1 = 1.0/PI; a2 = a1/(h*h); aa = a1/(h*h*h);
        adh=a2; f1 = -2.*aa/gauss2; deltap2=0.5; Wdeltap=a2*exp(-deltap2/gauss2); od_Wdeltap=1.0/Wdeltap;
    }
    else if(i_kernel==2) //2D Quadratic
    {
        a1 = 10.0/(PI*7.0); a2 = a1/(h*h); aa = a1/(h*h*h);
        b1 = 21.0*a2/20.0; e1 = 21.0*aa/20.0; adh=b1;
    }
    else if(i_kernel==3) //2D Cubic Spline
    {
        a1 =10.0/(PI*7.0); a2 =a1/(h*h); aa =a1/(h*h*h); adh=a2; a24=0.25*a2;
        c1 =-3.0*aa; d1=9.0*aa/4.0; c2 =-3.0*aa/4.0; deltap=1.0/1.5;
        Wdeltap=a2*(1.0-1.5*deltap*deltap+0.75*deltap*deltap*deltap); od_Wdeltap=1.0/Wdeltap;
    }
    else if(i_kernel==5)  //2D Wendland Quintic
    {
        Awen=0.557/(h*h); Bwen=-2.7852/(h*h*h); adh=Awen; deltap=0.5;
        Wdeltap=Awen*pow((1.0-deltap*0.5),4.0)*(2.0*deltap+1.0);
        od_Wdeltap=1.0/Wdeltap;
    }
    else {printf("i_kernel is not a valid options \n");exit(0);}

    //viscosity term eta2
    eta=0.1*h; eta2=eta*eta;
    //determine xmin,xmax,ymin,ymax,zmin,zmax
    xmin=xp[1]; xmax=xp[1];zmin=zp[1]; zmax=zp[1];
    for(i=2;i<=nb;i++)
    {
        if(xp[i]>xmax)xmax=xp[i]; if(xp[i]<xmin)xmin=xp[i];
        if(zp[i]>zmax)zmax=zp[i]; if(zp[i]<zmin)zmin=zp[i];
    }
    //Limits of domain -
    xmax_container=xmax; xmin_container=xmin;
    zmax_container=zmax; zmin_container=zmin;

    for(i=nbp1;i<=np;i++)
    {
        if(xp[i]>xmax)xmax=xp[i]; if(xp[i]<xmin)xmin=xp[i];
        if(zp[i]>zmax)zmax=zp[i]; if(zp[i]<zmin)zmin=zp[i];
    }
    xmax=xmax+0.1*h; xmin=xmin-0.1*h;
    zmax=zmax+0.1*h; zmin=zmin-0.1*h;

    //To place particles in a convenient place that leave the domain
    xtrash=xmax_container-0.005; ztrash=0.005;

    //determine number of cells
    one_over_h=1.0/h; one_over_2h=1.0/(2.0*h);
    ncx = (int)((xmax_container-xmin_container)*one_over_2h)+1;
    ncz = (int)((zmax_container-zmin_container)*one_over_2h)+1;

    nct=ncx*ncz; printf("ncx= %d, ncz= %d, nct= %d \n",ncx,ncz,nct);
    if(nct>nct_max){printf("ERROR: nct= %d >nct_max! \n",nct);exit(0);} //No. of cells exceeds nct_max

    pout("PART_00000.dat");
}

char num_char(int a)
{
	switch(a)
	{
		case 0: return '0'; break;
		case 1: return '1'; break;
		case 2: return '2'; break;
		case 3: return '3'; break;
		case 4: return '4'; break;
		case 5: return '5'; break;
		case 6: return '6'; break;
		case 7: return '7'; break;
		case 8: return '8'; break;
		case 9: return '9'; break;
		default: return '0';break;
	}
}

void get_save_name(void)
{
	int i,j;
	char a[5];
	i=ngrab;
	a[0]=num_char(i/10000);
	a[1]=num_char((i%10000)/1000);
	a[2]=num_char((i%1000)/100);
	a[3]=num_char((i%100)/10);
	a[4]=num_char(i%10);
	//shape
	for(j=0;j<5;j++)
	{ supp[j]=a[j]; detsupp[j]=a[j];}
}
/*
int main(void)
{
    getdata();
    pout("part_0000.dat");
    exit(0);
}
*/

