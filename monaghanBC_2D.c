
//accounts  for  Monaghan’s repulsive force between fluid and boundary particles

void monaghanBC(int iip,int jjp,double drx,double drz,double dot,double duxp,double duzp,double *fxbp,double *fzbp)
{
    double xpp,q,pdot3,drx_proj,drz_proj,deltappt,u_normal,App,eps_min,epsilon_Z,dyn_head,epsilon_dyn,Appc,facp;
    //iip is a WaterParticle(WP), jjp is a BoundaryParticle(BP), drxp=xp(WP)-xp(BP), drzp=zp(WP)-zp(BP)
    //Projection of WP position onto normal
    pdot3= drx*xnb[jjp]+drz*znb[jjp];  //SPHysics Guide, psi in Eq. 2.13
    drx_proj=drx-(pdot3)*xnb[jjp]; drz_proj=drz-(pdot3)*znb[jjp];
    //Projection of plane-projected vector r(BP->WP projection) onto tangent, SPHysics Guide, ksi in Eq. 2.12 & 2.15
    xpp=drx_proj*xtb[jjp]+drz_proj*ztb[jjp];
    if(xpp>0.0)deltappt=deltaptb[jjp][1]; else deltappt=deltaptb[jjp][0]; //SPHysics Guide, delta_b in Eq. 2.15

    if(pdot3>0.0) q=pdot3*one_over_2h;       //particle inside, normalized distance from wall
    else if(pdot3<0.0) q=0.2*fabs(pdot3)*one_over_2h;   //Particle ouside, normalized distance from wall
    else q=h*(1.0E-3)*one_over_2h;      //Particle on the boundary, correcting for zero distance from wall

    if(q<1.0&&pdot3>0.0&&fabs(xpp)<deltappt)
    {
        u_normal=duxp*xnb[jjp]+duzp*znb[jjp];  //流体粒子速度在边界法向上的分量
        //u_normal < 0 if particles are approaching,   u_normal > 0 if particles are moving apart
        if(u_normal>0.0)
            App=(0.01*cs[iip]*cs[iip]*one_over_beta_coef2-cs[iip]*u_normal*one_over_beta_coef)*one_over_h;
        else
            App=0.01*(cs[iip]*cs[iip]*one_over_beta_coef2)*one_over_h;     //SPHysics Guide, Eq. 2.14
        //--- Depth Correction ---      !SPHysics Guide, Eq. 2.17
        eps_min=0.02; epsilon_Z=fabs(1.0-(zp[jjp]/h_SWL));
        dyn_head=0.0; //=(up[jjp]*up[jjp]+wp[jjp]*wp[jjp])/(2.0*fabs(grz)*h_SWL);

        if(epsilon_Z>1.0-eps_min)epsilon_Z=1.0; //Close to the bottom (No correction )
        else if(zp[iip]>h_SWL)epsilon_Z=eps_min+dyn_head; //over still water level
        else  epsilon_Z = eps_min+epsilon_Z+dyn_head;       //under still water level & far from bottom

        //--- Velocity Correction ---     !SPHysics Guide, Eq. 2.18
        if(u_normal>0.0)epsilon_dyn=0.0;
        else epsilon_dyn=20.0*fabs(u_normal)/cs0;     //factor 20 needed to keep scheme stable with Monaghan BC wavemaker
        if(epsilon_dyn>1.0)epsilon_dyn=1.0;

        Appc = App*(epsilon_Z+epsilon_dyn); //Force correction     !SPHysics Guide, Eq. 2.16
        facp=0.5*(1.0+cos(PI*fabs(xpp)/deltappt))*(Appc*(1.0-q)/sqrt(q)); //SPHysics Guide, Eq. 2.12, 2.13, 2.15

        *fxbp = xnb[jjp]*facp;  //SPHysics Guide, Eq. 2.12
        *fzbp = znb[jjp]*facp;  //SPHysics Guide, Eq. 2.12
        //if(zp[jjp]<0.02)
        //fprintf(fp,"%lf %lf %lf %lf %lf %lf \n",xp[jjp],zp[jjp],xp[iip],zp[iip],*fxbp,*fzbp);//b_x,b_z,p_x,p_z,f_x,f_z
    }
    else
    {   *fxbp=0.0; *fzbp=0.0; }
}
//c0 = coeff.Vmax (Vmax = sqrt(g.h_SWL); coef = beta_coef.10, beta_coef2 = beta_coef**2
             //beta_coef2 = Factor that B is changed by
