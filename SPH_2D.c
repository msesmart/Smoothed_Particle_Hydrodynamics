
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
//#include "common_2D.h"
#include "parameters.c"
#include "getdata_2D.c"
#include "divide_2D.c"  //have included ini_divide(kind_P), keep_list()
#include "kernel_gaussian_2D.c"  //have included "kernel_correction_KC_2D.c"
#include "viscosity_laminar_2D.c"
#include "monaghanBC_2D.c"
#include "celij_self_2D.c"
#include "ac_2D.c"
#include "correct_2D.c"
#include "step_predictor_corrector_2D.c"
#include "densityFilter_Shepard_2D.c"

//Axis
//			^
//			|      / Y
//			| Z   /
//			|    /
//			|   /
//			|  /
//			| /
//			|/
//			--------------> X

void main()
{
    char name[40];
    /*** Initialization ***********************/
    getdata(); ini_divide(0); divide(1,nbf,0); keep_list(); printf("initia finished \n");//exit(0);
    ngrab=0; ngrabdet=0; tdetail=0.0; grx=0.0; grz=-g;
    time(&time_begin);

    //fp=fopen("BoundaryForce.txt","w"); fprintf(fp,"b_x b_z p_x p_z f_x f_z \n");

    while(Time<tmax)
    {
        visc_dt=0.0;
        if(Time-trec_ini-out*ngrab>out)ipoute=1; else ipoute=0;

        step_Predictor();

        itime=itime+1; Time=Time+dt; //printf("Time = %lf \n",Time);
        if(ipoute==1){ ngrab++; get_save_name(); strcpy(name,"PART_"); strcat(name,supp); pout(name); }
        if(i_densityFilter>0&&(itime%ndt_FilterPerform==0))densityFilter();
    }
    time(&time_end); printf("Time of Operation is : %d seconds\n",time_end-time_begin);
    //fclose(fp);
}

