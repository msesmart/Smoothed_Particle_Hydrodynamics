void pre_self_Shepard(int j1,int kind_p1,int ini_kind_p2)
{
    int kind_p2,jj_start,ii,i,jj,j,nstart_Shepard=nbp1;double drx,drz,rr2,V_i,V_j;
    for(kind_p2=ini_kind_p2;kind_p2<2;kind_p2++)
    {
        jj_start=1;
        if(kind_p1>kind_p2)jj_start=nplink_max+1;
        for(ii=1;ii<=nc[j1][kind_p1];ii++)
        {
            i = ibox[j1][kind_p1][ii];
            if(kind_p1==kind_p2)jj_start=jj_start+1;
            for(jj=jj_start;jj<=nc[j1][kind_p2];jj++)
            {
                j=ibox[j1][kind_p2][jj];
                if(i>=nstart_Shepard&&j>=nstart_Shepard)
                {
                    drx=xp[i]-xp[j]; drz=zp[i]-zp[j]; rr2=drx*drx+drz*drz;
                    if(rr2<fourh2&&rr2>1.0E-18)kernel(drx,drz,i,j,j1,j1,rr2);
                    //Sum_wab to calculate the normalized kernel
                    if(j>=nstart_Shepard)
                    {
                        V_j=pVol[j]; //=pm[j]/rhop[j]
                        sum_wab[i]=sum_wab[i]+ Wab*V_j;
                        rhop_sum[i]=rhop_sum[i]+pm[j]*Wab;
                    }
                    if(i>=nstart_Shepard)
                    {
                        V_i=pVol[i]; //=pm[i]/rhop[i]
                        sum_wab[j]=sum_wab[j]+Wab*V_i;
                        rhop_sum[j]=rhop_sum[j]+pm[i]*Wab;
                    }
                }
            }
        }
    }
}

void pre_celij_Shepard(int j1,int j2,int kind_p1,int ini_kind_p2,int lx2)
{
    int kind_p2,i,ii,j,jj,nstart_Shepard=nbp1;
    double rr2,drx,drz,V_j,V_i;
    for(kind_p2=ini_kind_p2;kind_p2<2;kind_p2++)
    {
        if(nc[j2][kind_p2]!=0)
        {
            for(ii=1;ii<=nc[j1][kind_p1];ii++)
            {
                i = ibox[j1][kind_p1][ii];
                for(jj=1;jj<=nc[j2][kind_p2];jj++)
                {
                    j = ibox[j2][kind_p2][jj];
                    if(i>=nstart_Shepard && j>=nstart_Shepard)
                    {
                        drx=xp[i]-xp[j]; drz=zp[i]-zp[j];
                        //periodicityCorrection(i,j,drx,drz,lx2);
                        rr2 = drx*drx + drz*drz;
                        if(rr2<fourh2&&rr2>1.0E-18)
                        {
                            kernel(drx,drz,i,j,j1,j2,rr2);
                            //Sum_wab to calculate the normalized kernel
                            if(j>=nstart_Shepard)
                            {
                                V_j=pVol[j];  //=pm(j)/rhop(j)
                                sum_wab[i]=sum_wab[i]+Wab*V_j;
                                rhop_sum[i]=rhop_sum[i]+pm[j]*Wab;
                            }

                            if(i>=nstart_Shepard)
                            {   V_i=pVol[i];   //=pm(i)/rhop(i)
                                sum_wab[j]=sum_wab[j]+Wab*V_i;
                                rhop_sum[j]=rhop_sum[j]+pm[i]*Wab;
                            }
                        }
                    }
                } //Box jj
            } //Box ii
        }  //Box jj is not empty
    }  //Kind of particle
}

void ac_Shepard(void)
{
    int i,kind_p1,ini_kind_p2,lz,lz2,lx,lx2,j1;

//***store useful arrays
    for(i=nstart;i<=np;i++)
    {   sum_wab[i]=0; rhop_sum[i]=0; }

    //ly2=0;   //Default Y-Periodic value

    for(kind_p1=0;kind_p1<2;kind_p1++)
    {
        ini_kind_p2=1-kind_p1;

        for(lz=1;lz<=ncz;lz++)
        for(lx=1;lx<=ncx;lx++)
        {
            j1=lx+(lz-1)*ncx;

            //if the cell is not empty, then loop over it and over neighboring cells
            if(nc[j1][kind_p1]>0)
            //Cells in the same XY sheet -
            {
                lx2=lx+1;
                if(lx2<=ncx)pre_celij_Shepard(j1,j1+1,kind_p1,ini_kind_p2,lx2);    //East

                lz2=lz+1;
                if(lz2<=ncz)
                {
                    pre_celij_Shepard(j1,j1+ncx,kind_p1,ini_kind_p2,lx2);   //North
                    lx2=lx+1;
                    if(lx2>2)pre_celij_Shepard(j1,j1+ncx-1,kind_p1,ini_kind_p2,lx2); //North-West
                    lx2=lx+1;
                    if(lx2<=ncx)pre_celij_Shepard(j1,j1+ncx+1,kind_p1,ini_kind_p2,lx2);   //North-East
                }
            }
        }

        for(j1=1;j1<=nct;j1++)
        {   if(nc[j1][kind_p1]>0)pre_self_Shepard(j1,kind_p1,ini_kind_p2); }


//-- Periodic Boundary Calls in X-Direction --
/*      if(i_periodicOBs[0]==1)
        {   //Special Treatment for Right Column Cells (lx=ncx) -lx = ncx
            lx2 = lx+1;
            for(lz=1;lz<=ncz-1;lz++)
            {
                j1=lx+(lz-1)*ncx;
                pre_celij_Shepard(j1,j1+ncall11,kind_p1,ini_kind_p2,lx2); //East X-Periodic
                pre_celij_Shepard(j1,j1+ncall13,kind_p1,ini_kind_p2,lx2); //North-East X-Periodic
            }

         //Special Treatment for Corner Cell (lx=ncx, lz=ncz) -
         //Note this is for X-Periodicity Only!               -
         //Norths & North-Easts of lx = 1
            lz = ncz;
            j1 = lx+lz*ncx;
            pre_celij_Shepard(j1,j1+ncall11,kind_p1,ini_kind_p2,lx2); //East X-Periodic
         //Special Treatment for Left Column Cells (lx=1) -
            lx = 0;
            lx2 = lx+1;
            for(lz=1;lz<=ncz-1;lz++)
            {
                j1=lx+(lz-1)*ncx;
                pre_celij_Shepard(j1,j1+ncall12,kind_p1,ini_kind_p2,lx2); //North-West X-Periodic
            }
        }*/
    }
}

void densityFilter(void)
{
    int i,nstart_Shepard;
//    itime_check=1E8;   //!0    !10e8   !0    !
//    i_ParticleCheck=23197; //  !41805   !231   !nb_local + 1   !3299   !
//    if(itime>=itime_check)rhop_old = rhop[i_ParticleCheck];
    nstart_Shepard=nbp1;   //nstart
    ac_Shepard();
/*
    if(itime>=itime_check)
    {
        i=i_ParticleCheck;
        sum_wab_check= sum_wab[i];
        rhop_sum_check= rhop_sum[i];
    }
*/
    for(i=nstart_Shepard;i<=np;i++)
    {
        sum_wab[i] = sum_wab[i]+ adh*pm[i]/rhop[i];  //self contribution
        rhop_sum[i] = rhop_sum[i]+ adh*pm[i];  //self contribution
        rhop[i] = rhop_sum[i]/sum_wab[i];
        pVol[i] = pm[i]/rhop[i];
        equation_of_state(rhop[i],TEp[i],&p[i],&cs[i],i_EoS);
    }
}
