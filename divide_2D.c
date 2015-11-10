
void ini_divide(int kind_p)
{
    int i,j;
    for(i=1;i<=nct;i++)
    {
        nc[i][kind_p]=0;
        for(j=1;j<=nplink_max;j++)
        ibox[i][kind_p][j]=0;
    }
}

void divide(int n_start,int n_end,int kind_p)
{
    int k,ii,j,icell,kcell;
    for(k=n_start;k<=n_end;k++)
    {
        if (iflag[k]!=0)
        {
            dx = xp[k]-0.0; dz = zp[k]-0.0; //dx = xp[k]-xmin; dz = zp[k]-zmin;
            icell=(int)(dx*one_over_2h )+1; kcell=(int)(dz*one_over_2h )+1;
//            if(i_periodicOBs[0]==1&&icell>ncx)icell=ncx;   //Periodic Open Boundary correction X-Direction
//            if(i_periodicOBs[2]==1&&kcell>ncz)kcell=ncz;   //Periodic Open Boundary correction Z-Direction

            ii=icell+(kcell-1)*ncx; //ii is the linear cell position in the matrix of 2h cells
            if(ii<1){printf("ERROR in divide_2D.f: ii<1 \n");exit(0);}
            if(ii>nct){printf("ERROR in divide_2D.f: ii>nct \n");exit(0);}

            nc[ii][kind_p]=nc[ii][kind_p]+1; //nc is the number of particles in cell ii
            if(nc[ii][kind_p]>nplink_max){printf("ERROR in divide_2D.f, nc(ii,kind_p)>= nplink_max \n");exit(0);}
            //Tells us that particle with array location k (i.e. xp(k)) is in box ii which, so far, contains "nc(ii,mkind)" particles
            j=nc[ii][kind_p]; ibox[ii][kind_p][j]=k;
        }
    }
}

void keep_list(void)
{
    int i;
    for(i=1;i<=nct;i++)keep_nc[i]=nc[i][0];
}






