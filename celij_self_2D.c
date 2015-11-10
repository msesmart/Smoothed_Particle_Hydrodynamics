
void celij(int j1,int j2,int kind_p1,int kind_p2,int lx2)
{
    int i,ii,j,jj; double rr2,temp,drx,drz,dux,duz,robar,one_over_rhobar,cbar,dot,dot2;
    double p_v,fab,Ra,Rb,term2i,term2j,term1i,term1j,pmi_Wab_over_rhobar,pmj_Wab_over_rhobar,fxbMON,fzbMON;

    if(nc[j2][kind_p2]!=0)
    {
        for(ii=1;ii<=nc[j1][kind_p1];ii++)
        {
            i=ibox[j1][kind_p1][ii];
            for(jj=1;jj<=nc[j2][kind_p2];jj++)
            {
                j=ibox[j2][kind_p2][jj];
                drx=xp[i]-xp[j]; drz=zp[i]-zp[j];
                rr2=drx*drx+drz*drz;
                if(rr2<fourh2&&rr2>1.E-18)  //fourh2=4.0*h*h
                {
                    dux=up[i]-up[j]; duz=wp[i]-wp[j];
                    kernel(drx,drz,i,j,j1,j2,rr2); kernel_correction(i,j); //Calculating kernel & Normalized Kernel Gradient
                    // average density
                    robar=0.5*( rhop[i]+rhop[j]); one_over_rhobar=2.0/(rhop[i]+rhop[j]); cbar =0.5*(cs[i]+cs[j]);
                    dot=drx*dux+drz*duz;// inner product rv
                    //Used to calculate the time step duu to viscosity
                    visc_dt=(dot/(rr2+eta2))>visc_dt?(dot/(rr2+eta2)):visc_dt;
                    // pressure and viscous force (Monaghan 1992; Ann. Rev.Astron. Astrop. 30. Formula 3.3) pm(j) is mass of particle j
                    p_v=pr[i]+pr[j];
                    if (index_tensile==1)
                    {
                        //Tensile correction
                        fab=Wab*od_Wdeltap; //NOTE: We'll use a non-normalized kernel to calculate tensile correction
                        fab=fab*fab*fab*fab; //It's the (Wab/Wdeltap)**4  of Monaghan's paper
                        if(p[i]>0)Ra=0.01*pr[i]; else Ra=0.2*fabs(pr[i]);
                        if(p[j]>0)Rb=0.01*pr[j]; else Rb=0.2*fabs(pr[j]);
                        p_v=p_v+(Ra+Rb)*fab;//R=Ra+Rb;
                    }

                    if(i>nb&&j>nb)   //both fluid particles
                    {
                        ax[i] = ax[i]-pm[j]*p_v*frxi; az[i] = az[i]-pm[j]*p_v*frzi;
                        ax[j] = ax[j]+pm[i]*p_v*frxj; az[j] = az[j]+pm[i]*p_v*frzj;

                        //gradients_calc(i,j,dux,duz);
                        viscosity(dot,drx,drz,dux,duz,rr2,cbar,robar,one_over_rhobar,i,j,j1,j2,&term2i,&term2j);
                        //density acceleration, using the derivative of the kernel, not the kernel itself
                        dot2=dux*frxi+duz*frzi; ar[i]=ar[i]+pm[j]*dot2;
                        dot2=dux*frxj+duz*frzj; ar[j]=ar[j]+pm[i]*dot2;

                        //Thermal Energy not used
                        //term1i=0.5*p_v*(frxi*dux+frzi*duz);term1j=0.5*p_v*(frxj*dux+frzj*duz);
                        //aTE[i]=aTE[i]+pm[j]*(term1i+term2i);aTE[j]=aTE[j]+pm[i]*(term1j+term2j);

                        //XSPH correction (Monaghan 1994;  J. Comp. Phys. 110. Formula 2.6)
                        pmj_Wab_over_rhobar = pm[j]*Wab*one_over_rhobar;
                        ux[i] = ux[i]-dux*pmj_Wab_over_rhobar;  //pm[j] * dux * Wab / robar ! (2.6)
                        wx[i] = wx[i]-duz*pmj_Wab_over_rhobar;  //pm[j] * duz * Wab / robar

                        pmi_Wab_over_rhobar = pm[i]*Wab*one_over_rhobar;
                        ux[j] = ux[j]+dux*pmi_Wab_over_rhobar;   //pm[i] * dux * Wab / robar
                        wx[j] = wx[j]+duz*pmi_Wab_over_rhobar;   //pm[i] * duz * Wab / robar
                    }
                    else if(i>nb&&j<=nb)    //j is boundary particle
                    {
                        //accounts  for  Monaghan’s repulsive force between fluid and boundary particles
                        monaghanBC(i,j,drx,drz,dot,dux,duz,&fxbMON,&fzbMON);
                        ax[i]=ax[i]+fxbMON; az[i]=az[i]+fzbMON; //*uno_sum_wab_i
                        temp = pVol[i]*visc_wall*((drx*frxi+drz*frzi)/(rr2+eta2));
                        ax[i] = ax[i]+temp*dux; az[i]=az[i]+temp*duz;
                    }
                    else if(j>nb&&i<=nb)    //i is boundary particle
                    {
                        monaghanBC(j,i,-drx,-drz,dot,-dux,-duz,&fxbMON,&fzbMON);
                        ax[j]=ax[j]+fxbMON; az[j]=az[j]+fzbMON;
                        temp = pVol[j]*visc_wall*((drx*frxj+drz*frzj)/(rr2+eta2));
                        ax[j] = ax[j]-temp*dux; az[j]=az[j]-temp*duz;
                    }
                    else ;
                } //if q<2
            } // Box jj
        } // Box ii
    } // Box jj is not empty
}

void self(int j1,int kind_p1,int kind_p2)
{
    int jj_start,ii,i,jj,j; double rr2,temp,dux,duz,drx,drz,p_v,dot,robar,one_over_rhobar,cbar;
    double fab,Ra,Rb,term2i,term2j,term1i,term1j,dot2i,dot2j,pmi_Wab_over_rhobar,pmj_Wab_over_rhobar,fxbMON,fzbMON;

    jj_start=1;
    for(ii=1;ii<=nc[j1][kind_p1];ii++)
    {
        i=ibox[j1][kind_p1][ii];
        if(kind_p1==kind_p2)jj_start=jj_start+1;

        for(jj=jj_start;jj<=nc[j1][kind_p2];jj++)
        {
            j=ibox[j1][kind_p2][jj];
            drx = xp[i]-xp[j]; drz = zp[i]-zp[j];
            rr2 = drx*drx + drz*drz;
            if(rr2<fourh2&&rr2>1.0E-18)
            {
                dux=up[i]-up[j]; duz=wp[i]-wp[j];

                //Calculating kernel & Normalized Kernel Gradient
                kernel(drx,drz,i,j,j1,j1,rr2); kernel_correction(i,j);

                //...  average density
                robar=0.5*(rhop[i]+rhop[j]); one_over_rhobar = 2.0/(rhop[i]+rhop[j]);
                cbar = 0.5*(cs[i]+cs[j]);
                //...  inner product rv
                dot = drx*dux + drz*duz;

                //Used to calculate the time step due to viscosity
                visc_dt=dot/(rr2+eta2)>visc_dt?dot/(rr2+eta2):visc_dt;

                //...  pressure and viscous force (1992; 3.3)
                p_v = pr[i]+pr[j]; //+  pi_visc

                //for tensile correction (Monaghan JCP 2000) Only to be activated with cubic spline kernel
                if(index_tensile==1) // ... Tensile correction
                {
                    fab=Wab*od_Wdeltap;       //NOTE: We'll use a non-normalized kernel to calculate tensile correction
                    fab=fab*fab*fab*fab; //It's the (Wab/Wdeltap)**4  of Monaghan's paper
                    if(p[i]>0)Ra=0.01*pr[i];else Ra=0.2*fabs(pr[i]);
                    if(p[j]>0)Rb=0.01*pr[j];else Rb= 0.2*fabs(pr[j]);
                    p_v = p_v+(Ra+Rb)*fab;
                }

                if(i>nb&&j>nb)    //both fluid particles
                {
                    ax[i]=ax[i]-pm[j]*p_v*frxi; az[i]=az[i]-pm[j]*p_v*frzi;
                    ax[j]=ax[j]+pm[i]*p_v*frxj; az[j]=az[j]+pm[i]*p_v*frzj;

                    //gradients_calc(i,j,dux,duz);
                    viscosity(dot,drx,drz,dux,duz,rr2,cbar,robar,one_over_rhobar,i,j,j1,j1,&term2i,&term2j);

                    //density acceleration (1992; 3.9)
                    dot2i = dux*frxi + duz*frzi; dot2j = dux*frxj + duz*frzj;
                    ar[i] = ar[i] + pm[j]*dot2i; ar[j] = ar[j] + pm[i]*dot2j;

                    // ...  Thermal Energy Monaghan, JCP 110 (1994) 399- 406)
                    //term1i=0.5*p_v *(frxi*dux+frzi*duz); term1j=0.5*p_v *(frxj*dux+frzj*duz);
                    //aTE[i]=aTE[i]+pm[j]*(term1i+term2i); aTE[j]=aTE[j]+pm[i]*(term1j+term2j);

                    // ...  XSPH correction
                    pmj_Wab_over_rhobar=pm[j]*Wab*one_over_rhobar;
                    ux[i] = ux[i]-dux*pmj_Wab_over_rhobar;  //pm[j] * dux * Wab / robar
                    wx[i] = wx[i]-duz*pmj_Wab_over_rhobar;  //pm[j] * duz * Wab / robar
                    pmi_Wab_over_rhobar = pm[i]*Wab*one_over_rhobar;
                    ux[j] = ux[j] + dux*pmi_Wab_over_rhobar;   //pm[i] * dux * Wab / robar
                    wx[j] = wx[j] + duz*pmi_Wab_over_rhobar;   //pm[i] * duz * Wab / robar
                }
                else if(i>nb&&j<=nb)    //j is boundary particle
                {
                    //accounts  for  Monaghan’s repulsive force between fluid and boundary particles
                    monaghanBC(i,j,drx,drz,dot,dux,duz,&fxbMON,&fzbMON);
                    ax[i]=ax[i]+fxbMON; az[i]=az[i]+fzbMON; //*uno_sum_wab_i
                    temp = pVol[i]*visc_wall*((drx*frxi+drz*frzi)/(rr2+eta2));
                    ax[i] = ax[i]+temp*dux; az[i]=az[i]+temp*duz;
                }
                else if(j>nb&&i<=nb)    //i is boundary particle
                {
                    monaghanBC(j,i,-drx,-drz,dot,-dux,-duz,&fxbMON,&fzbMON);
                    ax[j]=ax[j]+fxbMON; az[j]=az[j]+fzbMON;
                    temp = pVol[j]*visc_wall*((drx*frxj+drz*frzj)/(rr2+eta2));
                    ax[j] = ax[j]-temp*dux; az[j]=az[j]-temp*duz;
                }
                else {printf("ERROR: self do not satisfy i>nb&&j>nb, codes stop \n");exit(0);}
            }
        }
    }
}

