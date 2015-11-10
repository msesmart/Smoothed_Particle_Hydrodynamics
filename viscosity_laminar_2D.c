
void viscosity(double dot,double drx,double drz,double dux,double duz,double rr2,double cbar,double robar,double one_over_rhobar,int i,int j,int j1,int j2,double *term2i,double *term2j)
{
    double tempi,tempj;
    //--- Viscous Diffusion terms (Lo & Shao 2002) ---
    tempi =2.0*viscos_val*one_over_rhobar*((drx*frxi+drz*frzi)/(rr2 + eta2));
    tempj =2.0*viscos_val*one_over_rhobar*((drx*frxj+drz*frzj)/(rr2 + eta2));

    ax[i] = ax[i] + pm[j]*tempi*dux; az[i] = az[i] + pm[j]*tempi*duz;
    ax[j] = ax[j] - pm[i]*tempj*dux; az[j] = az[j] - pm[i]*tempj*duz;

    *term2i=-0.5*tempi*(dux*dux+duz*duz); *term2j=-0.5*tempj*(dux*dux+duz*duz);
}


