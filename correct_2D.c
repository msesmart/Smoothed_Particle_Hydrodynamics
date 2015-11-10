
//This subroutine is called by STEP every time step.
//It basically  accounts  for  the  body  forces  and  XSPH  correction.
//(and  SPS  terms  are calculated if i_visos = 3)

void correct(void)
{
    int i;
    for(i=nbp1;i<=np;i++)
    {
        //account for body forces
        udot[i]=udot[i]+grx*iflag[i]; wdot[i]=wdot[i]+grz*iflag[i];
        //account for XSPH
        xdot[i]=up[i]+xcor[i]; zdot[i]=wp[i]+zcor[i];
    }
}





