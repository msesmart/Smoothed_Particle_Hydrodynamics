/****************************/
//need to include "ac_KC_2D.c"
/****************************/
// detect  the  position  of  particles  outside  the  computational  domain  and  relocate  them
void check_limits(void)
{
    double dxp_dble,differ; int i,ncases=0;
    zmax=0;
    for(i=nbp1;i<=np;i++)   //Check if any fluid particle is over the box
    {
        if(zp[i]>zmax)zmax=zp[i];
        if(iflag[i]==1 && (xp[i]>xmax_container || xp[i]<xmin_container || zp[i]<zmin_container))
        {
            iflag[i]=0; ncases=ncases+1;
            xp[i]=xtrash; zp[i]=ztrash; up[i]=0.0; wp[i]=0.0;
            uo[i]=0.0; wo[i]=0.0;
        }
    }
    if(ncases>=100){printf("ERROR : Too many particles outside the domain \n");exit(0);}
}

double variable_time_step(void)
{
    //Calculate time step Conditions given in Monaghan JCP 1989, Monaghan & Koss JWPCOE 1999
    int i_fa=-1000,i,nstart_VTS; double Fa,Fa_max=-1000.0,cs_max=-1000.0,vel_magnitude,Fa_sqrt;
    double dt_1,dt_2;

    if(iBC==1)nstart_VTS=nstart; else nstart_VTS=nbp1;

    for(i=nstart_VTS;i<=np;i++)
    {
        Fa=sqrt(ax[i]*ax[i]+az[i]*az[i]);
        if(Fa>=Fa_max)
        {   i_fa=i; Fa_max=Fa; }
        if(cs[i]>=cs_max)cs_max=cs[i];
    }
    if(i_fa>0)
    {
        Fa_sqrt=pow((ax[i_fa]+grx)*(ax[i_fa]+grx)+(az[i_fa]+grz)*(az[i_fa]+grz),0.25);
        dt_1=sqrt(h)/Fa_sqrt;
    }
    else
    {
        Fa_sqrt=pow(grx*grx+grz*grz,0.25);
        dt_1=sqrt(h)/Fa_sqrt;
    }

    //dt_2=h/(cs_max+h*visc_dt)
    if(cs_max>0.0)dt_2=h/(cs_max+h*visc_dt);
    else dt_2=dt_1;

    return CFL_number*(dt_1<dt_2?dt_1:dt_2);
}

void step_Predictor(void)
{
    int i,iDBC;
    double dt2;
    dt2=0.5*dt;
    //divide particles into cells
	ini_divide(1); divide(nbp1,np,1);
    ac();
    //compute corrections to:   (a) rate of change of velocity due to: (i)  body forces, and (ii) boundary forces
    //                          (b) rate of change of position (free-surface, XSPH)
    correct();

    if(ivar_dt==1)ddt_p=variable_time_step();
    // predictor
    for(i=nbp1;i<=np;i++)
    {
        xo[i]=xp[i]; zo[i]=zp[i]; uo[i]=up[i]; wo[i]=wp[i]; po[i]=p[i];
        rhoo[i]=rhop[i]; rhop[i]=rhoo[i]+dt2*rdot[i]; pVol[i]=pm[i]/rhop[i];
        equation_of_state(rhop[i],TEp[i],&p[i],&cs[i],i_EoS);
    }

    for(i=nbp1;i<=np;i++)
    {
        xp[i] = xo[i] + dt2*xdot[i];
        zp[i] = zo[i] + dt2*zdot[i];
        up[i] = uo[i] + dt2*udot[i];
        wp[i] = wo[i] + dt2*wdot[i];
    }

    check_limits();
    ini_divide(1); divide(nbp1,np,1);
    ac();
    correct();

    if(ivar_dt==1)
    {
        ddt_c=variable_time_step();
        dt=ddt_p<ddt_c?ddt_p:ddt_c; //Time step to be used in next loop
    }

    for(i=nbp1;i<=np;i++)
    {
        rhop[i]=rhoo[i]+dt2*rdot[i];
        xp[i] = xo[i] + dt2*xdot[i];
        zp[i] = zo[i] + dt2*zdot[i];
        up[i] = uo[i] + dt2*udot[i];
        wp[i] = wo[i] + dt2*wdot[i];
    }

// Perform final integration correction --
    for(i=nbp1;i<=np;i++)
    {
        rhop[i]=2.0*rhop[i]-rhoo[i];
        pVol[i]=pm[i]/rhop[i];
        //TEp[i]=2.0*TEp[i]-TEo[i];
        equation_of_state(rhop[i],TEp[i],&p[i],&cs[i],i_EoS);

        xp[i]=2.0*xp[i]-xo[i]; zp[i]=2.0*zp[i]-zo[i];
        up[i]=2.0*up[i]-uo[i]; wp[i]=2.0*wp[i]-wo[i];
    }
}
