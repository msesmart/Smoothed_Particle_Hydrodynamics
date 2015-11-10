void ac(void)   //ac_KC_2D.f
{
    int i,kind_p1,kind_p2,ini_kind_p2,lz,lx,lx2,lz2,j1; double one_over_sum_wab;

    //ac_main();
    for(i=nstart;i<=np;i++)  //nstart=nbp1
    {
        pr[i]=p[i]/(rhop[i]*rhop[i]);  // pressure force
        // -- Zeroing Variables --
        ax[i] = 0.0; az[i] = 0.0; ar[i] = 0.0; ux[i] = 0.0; wx[i] = 0.0; aTE[i]=0.0;
        sum_wab[i] = 0.0; rhop_sum[i] = 0.0;
    }

    kind_p1=0;kind_p2=1;  // B. Cells affect P. Cells
    for(lz=1;lz<=ncz;lz++)
    for(lx=1;lx<=ncx;lx++)
    {
        j1 = lx+(lz-1)*ncx;
        if(nc[j1][kind_p1]>0) // if B. cell is not empty, then loop over it and over neighboring cells
        {
            // -- Cells in the same XY sheet -
            lx2=lx+1; if(lx2<=ncx)celij(j1,j1+1,kind_p1,kind_p2,lx2); //East
            lz2=lz+1;
            if(lz2<=ncz)
            {
                celij(j1,j1+ncx,kind_p1,kind_p2,lx2);  //North or up
                if(lx2>2)celij(j1,j1+ncx-1,kind_p1,kind_p2,lx2);   //North-West
                if(lx2<=ncx)celij(j1,j1+ncx+1,kind_p1,kind_p2,lx2);   //North-East
            }
            self(j1,kind_p1,kind_p2);
        }
    }

    kind_p1=1;kind_p2=0;
    for(lz=1;lz<=ncz;lz++)
    for(lx=1;lx<=ncx;lx++)
    {
        j1 = lx+(lz-1)*ncx;
        if(nc[j1][kind_p1]>0) // if the cell is not empty, then loop over it and over neighboring cells
        {
            // -- Cells in the same XY sheet -
            lx2=lx+1; lz2=lz+1;
            if(lx2<=ncx)celij(j1,j1+1,kind_p1,kind_p1,lx2); //East
            if(lz2<=ncz)
            {
                celij(j1,j1+ncx,kind_p1,kind_p1,lx2);   //North
                if(lx2>2)celij(j1,j1+ncx-1,kind_p1,kind_p1,lx2);   //North-West
                if(lx2<=ncx)celij(j1,j1+ncx+1,kind_p1,kind_p1,lx2);   //North-East
            }
            self(j1,kind_p1,kind_p1);
        }
    }

    for(i=nstart;i<=np;i++)
    {
        udot[i]=ax[i]; wdot[i]=az[i]; rdot[i]=ar[i];
        xcor[i]=eps*ux[i]; zcor[i]=eps*wx[i];
    }
    // End of ac_main();
    for(i=nbp1;i<=np;i++)
    {
        sum_wab[i]=sum_wab[i]+adh*pm[i]/rhop[i];  //self contribution
        if(sum_wab[i]<0.1)one_over_sum_wab=1.0; else one_over_sum_wab=1.0/sum_wab[i];
        xcor[i]=eps*ux[i]*one_over_sum_wab; zcor[i]=eps*wx[i]*one_over_sum_wab;
    }
}
