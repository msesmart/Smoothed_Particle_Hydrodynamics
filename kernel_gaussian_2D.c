
void kernel(double drx,double drz,int i,int j,int j1,int j2,double rr2)
{
     double fac,qq,rad; //gauss2=1.0; a2=1.0/(PI*h*h); f1 = -2.0/(PI*h*h*h);
    //Gaussian 2D Kernel
    rad=sqrt(rr2); qq=rad/h; Wab=a2*exp(-qq*qq/gauss2);
    fac=f1*qq*exp(-qq*qq/gauss2)/rad;  //=partial derivative of Wab /rad;
    frx=fac*drx; frz=fac*drz;  //?
}

// kernel_correction_KC_2D.c
void kernel_correction(int i,int j)
{
    double V_j,V_i;
    //Sum_wab to calculate the normalized kernel
	if(j>nb)
    {
        V_j=pVol[j];  //=pm(j)/rhop(j)
        sum_wab[i]=sum_wab[i]+Wab*V_j;
        rhop_sum[i]=rhop_sum[i]+pm[j]*Wab;
    }
	if(i>nb)
    {
        V_i=pVol[i];   //=pm(i)/rhop(i)
        sum_wab[j]=sum_wab[j]+Wab*V_i;
        rhop_sum[j]=rhop_sum[j]+pm[i]*Wab;
    }
    frxi=frx; frzi=frz; frxj=frx; frzj=frz;
}
