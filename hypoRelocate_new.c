/***********************************************************************************
 * This is a earthquake relocation code named hypoRelocate
 * developed by Li Sun, sunli@seis.ac.cn or sunli831022@gmail.com
 * For reference, see:
 * Sun, L., Zhang, M., & Wen, L. (2016). A new method for high‐resolution event 
 *     relocation and application to the aftershocks of Lushan Earthquake, China. 
 *     Journal of Geophysical Research: Solid Earth.
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
***********************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "myvar.h"
#include <malloc.h>
#include <omp.h>
typedef struct velocitys 
{
    float vel;
    float dep;
} MVAL;

typedef struct dvals 
{
    double abtrms1;
    double abtrms2;
    double trms;
    double trmscc;
    double trms_n;
    double tmean1;
    double tmean2;
    double coda_res;
} DVAL;

typedef struct initial_vals 
{
    double stla;
    double stlo;
    double DT1;
    double DT2;
    double TT;
    char stname[10];
    char phase[10];
    int ev1;
    int ev2;
    double coda;
    int rel;
    int mod_id;
    int gid;
} IVAL;

typedef struct sta_res 
{
    double DD_DT;
    double DD_DT2;
    double s_rms;
    double s_rmsc;
    double ev1_r;
    double ev2_r;
} SRES;

typedef struct initial_evs 
{
    double evla;
    double evlo;
    double evdp;
    int evid;
    int evflag;
    int date;
    float sec;
} EIVAL;
typedef struct final_evs 
{
    double evla[2];
    double evlo[2];
    double evdp[2];
    double deot[2];
} FEIVAL;

typedef struct mod_dt
{
    double dtdd[101][1001];
    double dtdh[101][1001];
    double tt[101][1001];
} MDT;

double dist(double evla, double evlo, double stla, double stlo);
int Calculate_Res (int ggid, int Nst, int Est, int evf, double coef1, double coef2, double coef3, double *dis_la, float *vel, float *dep, MDT *DT,IVAL *IV, EIVAL *EIV, FEIVAL *FEIV, DVAL *DD, SRES *SR);
int Calculate_Res1 (int ggid, int Nst, int Est, int evf, double coef1, double coef2, double coef3, double *dis_la, IVAL *IV, EIVAL *EIV, FEIVAL *FEIV, DVAL *DD, SRES *SR);
int Calculate_Res2 (int ggid, int Nst, int Est, double coef1, double coef2, double coef3, double *dis_la, MDT *DT, IVAL *IV, EIVAL *EIV, FEIVAL *FEIV, DVAL *DD, SRES *SR);
int Calculate_Res21 (int ggid, int Nst, int Est, double coef1, double coef2, double coef3, double *dis_la, float *vel, float *dep, IVAL *IV, EIVAL *EIV, FEIVAL *FEIV, DVAL *DD, SRES *SR);
extern void ttime_(float *delta, float *depth, int *nl, float *v, float *top, float *t, float *ain);
int ReadFileM(char *name, MVAL *MIV);
int show_usage(){
  fprintf(stderr,"Usage:\n\thypoRelocate eventfile= timefile= velfile= svelfile= velfile2= outfile= maxla= maxlo= maxed= itnum= nmax= sacoef= coef1= coef2= coef3= \n");
  fprintf(stderr,"\tstandard output: outfile\n");
	return (-1);
}

int main(int ac,char **av){

    int nlays, nlays2, Nst1, Nst2, mm, vv;
    float tt,tain;
    float *vel, *dep,*svel, *sdep, *vel2, *dep2;
    double pi=3.1415926, lad; 
    int ist, i, n, imax, Nst, Est, it, iter, nnum;
    double evla, evlo, evdp, evla1, evla2, evlo1, evlo2, evdp1, evdp2; 
    double maxlo, maxla, maxed, maxdo, coef1, coef2, coef3;
    double sacoef, coda_tt, total_rms, total_rms0, e_tmean;
    int nmax, ggid, k, nnd, groupnum, groupnum0;
    int depi, gcai, count, evf, grn, iitt, tnst, itnum;
    float depp, gcaa;
    double t,t0=1.0,cdla,cdlo,cddp,cdlt, pp, tp;
    long int s;
    double trms_min, vpa;
    double *dis_la;
    FILE *fop;

    DVAL *DD;
    IVAL *IV, *TIV;
    MVAL *MIV;
    EIVAL *EIV;
    FEIVAL *FEIV;
    SRES *SR;
    MDT *DT;
    char eventfile[50], timefile[50], outfile[50], velfile[100], svelfile[100], velfile2[100];
    int ReadFile(char *name, IVAL *IV);
    int ReadEF(char *name, EIVAL *EIV);
    int WriteFile(char *name, IVAL *IV, DVAL *DD, int imax, int Nst);
    int WriteFile2(char *name, EIVAL *EIV, FEIVAL *FEIV, int Est);
    double uniform(double a,double b,long int *seed);
    double cauchyrnd(double alpha, double beta, long int *s);
    if(ac==1) {
      show_usage();
      exit(1);
    }

    setpar(ac,av);
    mstpar("eventfile", "s",   eventfile);
    mstpar("timefile", "s",   timefile);
    mstpar("velfile", "s",   velfile);
    mstpar("svelfile", "s",   svelfile);
    mstpar("velfile2", "s",   velfile2);
    mstpar("outfile",  "s",   outfile);
    getpar("Maxla",    "F",   &maxla);
    getpar("Maxlo",    "F",   &maxlo);
    getpar("Maxed",    "F",   &maxed);
    getpar("itnum", "d",   &itnum);
    getpar("nmax", "d",   &nmax);
    getpar("sacoef",    "F",   &sacoef);
    getpar("coef1",    "F",   &coef1);
    getpar("coef2",    "F",   &coef2);
    getpar("coef3",    "F",   &coef3);
    endpar();
/****Read Velocity Model 1*****/
    nlays = 100;

    if((MIV = malloc(sizeof(MVAL)*nlays))==NULL){
        fprintf(stderr,"malloc memory error for IV\n");
        return(-1);
    }
    nlays = ReadFileM(velfile, MIV);

    if ((vel = malloc(nlays * sizeof(*vel))) == NULL) {
              fprintf(stderr, "allocation failed for vel.\n");
    }
    if ((dep = malloc(nlays * sizeof(*dep))) == NULL) {
              fprintf(stderr, "allocation failed for dep.\n");
    }
    for (i = 0; i < nlays; i++){
	vel[i] = MIV[i].vel;
	dep[i] = MIV[i].dep;

    }
free(MIV);
/**********************************/
/****Read S Velocity Model*****/
    nlays = 100;
    if((MIV = malloc(sizeof(MVAL)*nlays))==NULL){
        fprintf(stderr,"malloc memory error for IV\n");
        return(-1);
    }
    nlays = ReadFileM(svelfile, MIV);
    if ((svel = malloc(nlays * sizeof(*svel))) == NULL) {
              fprintf(stderr, "allocation failed for vel.\n");
    }
    if ((sdep = malloc(nlays * sizeof(*sdep))) == NULL) {
              fprintf(stderr, "allocation failed for dep.\n");
    }
    for (i = 0; i < nlays; i++){
	svel[i] = MIV[i].vel;
	sdep[i] = MIV[i].dep;
    }
free(MIV);
/**********************************/
/**********************************/
/****Read Velocity Model 2*****/
    nlays2 = 100;
    if((MIV = malloc(sizeof(MVAL)*nlays))==NULL){
        fprintf(stderr,"malloc memory error for IV\n");
        return(-1);
    }
    nlays2 = ReadFileM(velfile2, MIV);
   if ((vel2 = malloc(nlays2 * sizeof(*vel2))) == NULL) {
              fprintf(stderr, "allocation failed for vel2.\n");
    }
    if ((dep2 = malloc(nlays2 * sizeof(*dep2))) == NULL) {
              fprintf(stderr, "allocation failed for dep2.\n");
    }

    for (i = 0; i < nlays; i++){
	vel2[i] = MIV[i].vel;
	dep2[i] = MIV[i].dep;	

    }
free(MIV);
/**********************************/

  if((DT = malloc(sizeof(MDT)*3))==NULL){
        fprintf(stderr,"malloc memory error for DT1\n");
        return(-1);
    }
 
/*********calculate dtdd dtdh table 1******/
    for (depi = 0; depi < 101; depi++){
	for(gcai = 0; gcai < 1001; gcai++){
		depp = depi;
		gcaa = gcai;
		if (depi == 0) 
			depp = depi+0.001;

		for (vv = 0; vv < nlays; vv++){
			if ((depp < dep[vv])&&(depp >= dep[vv-1])&& (vv > 0))
			    vpa = vel[vv-1];
		}
    		ttime_(&gcaa, &depp, &nlays, vel, dep, &tt, &tain);
		DT[0].tt[depi][gcai] = tt;
		DT[0].dtdd[depi][gcai] = sin((tain*pi)/180)/vpa*111.2;
		DT[0].dtdh[depi][gcai] = cos((tain*pi)/180)/vpa;

	}
    }


/********************************************/
/*********calculate S wave dtdd dtdh table******/
    for (depi = 0; depi < 51; depi++){
	for(gcai = 0; gcai < 1001; gcai++){
		depp = depi;
		gcaa = gcai;
		if (depi == 0) 
			depp = depi+0.001;

		for (vv = 0; vv < nlays; vv++){
			if ((depp < sdep[vv])&&(depp >= sdep[vv-1])&& (vv > 0)){
			    vpa = svel[vv-1];
			    break;
			}
		}
    		ttime_(&gcaa, &depp, &nlays, svel, sdep, &tt, &tain);
		DT[2].tt[depi][gcai] = tt;
		DT[2].dtdd[depi][gcai] = sin((tain*pi)/180)/vpa*111.2;
		DT[2].dtdh[depi][gcai] = cos((tain*pi)/180)/vpa;
	}
    }
/********************************************/
/*********calculate dtdd dtdh table 2******/
    for (depi = 0; depi < 101; depi++){
	for(gcai = 0; gcai < 1001; gcai++){
		depp = depi;
		gcaa = gcai;
		if (depi == 0) 
			depp = depi+0.001;

		for (vv = 0; vv < nlays2; vv++){
			if ((depp < dep2[vv])&&(depp >= dep2[vv-1])&& (vv > 0))
			    vpa = vel2[vv-1];
		}
    		ttime_(&gcaa, &depp, &nlays, vel2, dep2, &tt, &tain);
		DT[1].tt[depi][gcai] = tt;
		DT[1].dtdd[depi][gcai] = sin((tain*pi)/180)/vpa*111.2;
		DT[1].dtdh[depi][gcai] = cos((tain*pi)/180)/vpa;

	}
    }
 
/********************************************/
/**********cal change of dist per deg with lat ******************/
   if ((dis_la = malloc(901 * sizeof(*dis_la))) == NULL) {
              fprintf(stderr, "allocation failed for dis_la.\n");
    }
	for(i = 0; i < 901; i++){
		lad = 0.1 * i;
		dis_la[i]=dist(lad,103,lad,104)*100;

	}

    Est = 2000;

    if((EIV = malloc(sizeof(EIVAL)*Est))==NULL){
        fprintf(stderr,"malloc memory error for EIV\n");
        return(-1);
    }
    Est = ReadEF(eventfile, EIV);

    if((FEIV = malloc(sizeof(FEIVAL)*Est))==NULL){
        fprintf(stderr,"malloc memory error for FEIV\n");
        return(-1);
    }

    Nst = 2000000;
    if((IV = malloc(sizeof(IVAL)*Nst))==NULL){
        fprintf(stderr,"malloc memory error for IV\n");
        return(-1);
    }
    Nst = ReadFile(timefile, IV);

    if((TIV = malloc(sizeof(IVAL)*(Nst*10/Est)))==NULL){
        fprintf(stderr,"malloc memory error for TIV\n");
        return(-1);
    }
    if((SR = malloc(sizeof(SRES)*Nst))==NULL){
        fprintf(stderr,"malloc memory error for SR\n");
        return(-1);
    }
    groupnum = IV[Nst-1].gid;
    if((DD = malloc(sizeof(DVAL)*(groupnum+1)))==NULL){
        fprintf(stderr,"malloc memory error for DD\n");
        return(-1);
    }

	for(nnd=0;nnd<Est;nnd++){
	    FEIV[nnd].evla[1] = EIV[nnd].evla;
	    FEIV[nnd].evlo[1] = EIV[nnd].evlo;
	    FEIV[nnd].evdp[1] = EIV[nnd].evdp;
	    FEIV[nnd].deot[1] = 0;
	    FEIV[nnd].evla[0] = EIV[nnd].evla;
	    FEIV[nnd].evlo[0] = EIV[nnd].evlo;
	    FEIV[nnd].evdp[0] = EIV[nnd].evdp;
	    FEIV[nnd].deot[0] = 0;
	}
/***********Iteration Begin*********************************************************/
      for (iitt= 0;iitt < itnum; iitt++){
		printf("Iteration %4d\n",iitt);
	for (nnd=0;nnd<Est;nnd++){
	        trms_min = 1.e25;
		evf = EIV[nnd].evid;
		tnst = 0;
		for(i=0;i<Nst;i++){
			if (IV[i].ev1 == evf || IV[i].ev2 == evf) {
				TIV[tnst] = IV[i];
				tnst += 1;
			}
		}		
		fflush(stdout);
		  srand(time(NULL));	
	      for(iter=0;iter<nmax;iter++){
			count =0;
		    t = t0*pow(sacoef,iter);
round:	  	    s = rand();
/************for fix depth********************************/
//deott:	    	    cdlt = cauchyrnd(0,1,&s);
//	    	    FEIV[nnd].deot[0] = FEIV[nnd].deot[1] + t*maxdo*cdlt;
//		    if (FEIV[nnd].deot[0]<(-maxdo) || FEIV[nnd].deot[0]> maxdo)
//			goto deott;	
lati:	    	    cdla = cauchyrnd(0,1,&s);
	    	    FEIV[nnd].evla[0] = FEIV[nnd].evla[1] + t*maxla*cdla;	
		    if (FEIV[nnd].evla[0]<(EIV[nnd].evla-maxla) || FEIV[nnd].evla[0]> (EIV[nnd].evla+maxla))
			goto lati;
longi:	    	    cdlo = cauchyrnd(0,1,&s);
	    	    FEIV[nnd].evlo[0] = FEIV[nnd].evlo[1] + t*maxlo*cdlo;
		    if (FEIV[nnd].evlo[0]<(EIV[nnd].evlo-maxlo) || FEIV[nnd].evlo[0]>(EIV[nnd].evlo+maxlo))
			goto longi;
		  if (EIV[nnd].evflag == 0){
evddp:	    	    cddp = cauchyrnd(0,1,&s);
	    	    FEIV[nnd].evdp[0] = FEIV[nnd].evdp[1] + t*maxed*cddp;
		    if (FEIV[nnd].evdp[0]< 0 || FEIV[nnd].evdp[0]>(EIV[nnd].evdp+maxed))
			goto evddp;
		  }
		 else if(EIV[nnd].evflag == 1){
	    	    FEIV[nnd].evdp[0] = FEIV[nnd].evdp[1];
		 }
/*
		 else if(EIV[nnd].evflag == 1){
evddp2:	    	    cddp = cauchyrnd(0,1,&s);
	    	    FEIV[nnd].evdp[0] = FEIV[nnd].evdp[1] + t*maxed*cddp;
		    if (FEIV[nnd].evdp[0]<(EIV[nnd].evdp-maxed) || FEIV[nnd].evdp[0]>(EIV[nnd].evdp+maxed))
			goto evddp2;
		 }

/***************************************************************/
/*******For fix lati longi depth*********************/
/*
		  if (EIV[nnd].evflag == 0){
lati:	    	    cdla = cauchyrnd(0,1,&s);
	    	    FEIV[nnd].evla[0] = FEIV[nnd].evla[1] + t*maxla*cdla;	
		    if (FEIV[nnd].evla[0]<(EIV[nnd].evla-maxla) || FEIV[nnd].evla[0]> (EIV[nnd].evla+maxla))
			goto lati;
longi:	    	    cdlo = cauchyrnd(0,1,&s);
	    	    FEIV[nnd].evlo[0] = FEIV[nnd].evlo[1] + 0*t*maxlo*cdlo;
		    if (FEIV[nnd].evlo[0]<(EIV[nnd].evlo-maxlo) || FEIV[nnd].evlo[0]>(EIV[nnd].evlo+maxlo))
			goto longi;
evddp:	    	    cddp = cauchyrnd(0,1,&s);
	    	    FEIV[nnd].evdp[0] = FEIV[nnd].evdp[1] + t*maxed*cddp;
		    if (FEIV[nnd].evdp[0]<(EIV[nnd].evdp-maxed) || FEIV[nnd].evdp[0]>(EIV[nnd].evdp+maxed))
			goto evddp;
		  }
		 else if(EIV[nnd].evflag == 1){
		    FEIV[nnd].evla[0] = FEIV[nnd].evla[1];	
	    	    FEIV[nnd].evlo[0] = FEIV[nnd].evlo[1];
	    	    FEIV[nnd].evdp[0] = FEIV[nnd].evdp[1];		
		 }

/******************************************************************/
	      
	      total_rms = 0;
		grn = 0;
		
	      for(ggid = 1;ggid <= groupnum;ggid++){
		
			Calculate_Res (ggid, tnst, Est, evf, coef1, coef2, coef3, dis_la, svel, sdep, DT, TIV, EIV, FEIV, DD, SR);
//			Calculate_Res1 (ggid, tnst, Est, evf, coef1, coef2, coef3, dis_la, vel, dep, DT, TIV, EIV, FEIV, DD, SR);
	     }

	      for(ggid = 1;ggid <= groupnum;ggid++){
		if (DD[ggid].trms_n > 0.000001 ){
		  if(DD[ggid].trms_n <= 0.2){
		      total_rms += DD[ggid].trms_n*10;
		      grn += 10;
		  }
		  else if((DD[ggid].trms_n > 0.2)&&(DD[ggid].trms_n <= 0.4)){
		      total_rms += DD[ggid].trms_n*5;
		      grn += 5;
		  }
		  else if((DD[ggid].trms_n > 0.4)&&(DD[ggid].trms_n <= 0.6)){
		      total_rms += DD[ggid].trms_n*3;
		      grn += 3;
		  }
		  else if((DD[ggid].trms_n > 0.6)&&(DD[ggid].trms_n <= 0.8)){
		      total_rms += DD[ggid].trms_n*2;
		      grn += 2;
		  }
		  else if((DD[ggid].trms_n > 0.8)){
		      total_rms += DD[ggid].trms_n*1;
		      grn += 1;
		  }
		}
		
	      }
	      if (grn == 0) 
		total_rms = 0;
	      else 
	      	total_rms = total_rms/(grn);

	      if (trms_min < total_rms) {
		pp = exp((trms_min - total_rms)/(t*maxla));
	  	s = rand();
		tp = uniform(0.0, 1.0, &s);
	      }
	      else {
		tp = 1.e10;
	      }


               if((trms_min >= total_rms)||(pp > fabs(tp))){
		   trms_min = total_rms;
			FEIV[nnd].evla[1] = FEIV[nnd].evla[0];
			FEIV[nnd].evlo[1] = FEIV[nnd].evlo[0];
			FEIV[nnd].evdp[1] = FEIV[nnd].evdp[0];
			FEIV[nnd].deot[1] = FEIV[nnd].deot[0];
			if (iitt == (itnum -1)){
				printf("Res %-5d %lf  %lf  %lf  %lf  %lf\n",EIV[nnd].evid, FEIV[nnd].deot[0],FEIV[nnd].evlo[0], FEIV[nnd].evla[0], FEIV[nnd].evdp[0], total_rms);
			}
               }
	       else {
		count ++;
		if (count > 1000)
			continue;
		else
			goto round;
	       }
	}
		printf("%-2d Ev %4d rms  %lf\n",iitt,evf,total_rms);
    }
		for(nnd= 0;nnd<Est;nnd++) {

			printf("event %-5d %lf  %lf  %lf  %lf  %lf  %lf  %lf\n",EIV[nnd].evid,EIV[nnd].evlo, EIV[nnd].evla, EIV[nnd].evdp, FEIV[nnd].deot[0],FEIV[nnd].evlo[0], FEIV[nnd].evla[0], FEIV[nnd].evdp[0]);
		}

/***********************residual print******************************************/
 		printf("\n");
 		printf("group_num  rela_rms  rela_rmscc  ev1_rms  ev2_rms gourp_total_rms\n");

    	       for(ggid = 1;ggid <= groupnum;ggid++){

			Calculate_Res2 (ggid, Nst, Est, coef1, coef2, coef3, dis_la, DT, IV, EIV, FEIV, DD, SR);
			Calculate_Res21 (ggid, Nst, Est, coef1, coef2, coef3, dis_la, svel, sdep, IV, EIV, FEIV, DD, SR);
	       }

		grn = 0;
		total_rms = 0;
	      for(ggid = 1;ggid <= groupnum;ggid++){
		if (DD[ggid].trms_n > 0.000001 ){

		  if(DD[ggid].trms_n <= 0.2){
		      total_rms += DD[ggid].trms_n*10;
		      grn += 10;
		  }
		  else if((DD[ggid].trms_n > 0.2)&&(DD[ggid].trms_n <= 0.4)){
		      total_rms += DD[ggid].trms_n*5;
		      grn += 5;
		  }
		  else if((DD[ggid].trms_n > 0.4)&&(DD[ggid].trms_n <= 0.6)){
		      total_rms += DD[ggid].trms_n*3;
		      grn += 3;
		  }
		  else if((DD[ggid].trms_n > 0.6)&&(DD[ggid].trms_n <= 0.8)){
		      total_rms += DD[ggid].trms_n*2;
		      grn += 2;
		  }
		  else if((DD[ggid].trms_n > 0.8)){
		      total_rms += DD[ggid].trms_n*1;
		      grn += 1;
		  }
		}


	      }
	      if (grn == 0) 
		total_rms = 0;
	      else 
	      	total_rms = total_rms/(grn);
	    printf("total rms %lf\n",total_rms);
		for(mm = 1;mm <= groupnum;mm++){
			printf("group %-5d    %lf  %lf  %lf  %lf  %lf\n",mm, DD[mm].trms, DD[mm].trmscc, DD[mm].abtrms1, DD[mm].abtrms2, DD[mm].coda_res);

		    for(nnd = 0; nnd < Nst; nnd++) {
                        if(IV[nnd].rel >= 0 && IV[nnd].gid == mm){
				printf("sta %s prms %10.6f  ev1rms %10.6f ev2rms %10.6f\n",IV[nnd].stname, SR[nnd].s_rms, SR[nnd].ev1_r, SR[nnd].ev2_r);
			}
		    }

		}
			printf("\n");

   }




/**********************************************************************************/
	WriteFile2(outfile, EIV, FEIV, Est);
//    free(EIV);
//    free(FEIV);
//    free(IV);
//    free(DD);
}





int Calculate_Res (int ggid, int Nst, int Est, int evf, double coef1, double coef2, double coef3, double *dis_la, float *vel, float *dep, MDT *DT, IVAL *IV, EIVAL *EIV, FEIVAL *FEIV, DVAL *DD, SRES *SR) {
	int Nst0, ist, i, j, k, iid, ev1, ev2, depi, gcai;
	int Nst1, Nst2, Nst0c,laa, cottn =0, nnl, vf;
	float nvel[100],ndep[100], ain, scoda_tt, hdis, GCarc;
	double lat0, lon0, lat1, lon1, lat2, lon2, lat3, lon3;
	double gca, a1, a2, dep1, depth1, depth2, dtdd, dtdh, tdep;
	double tt1[Nst], tt2[Nst], GCarc1, GCarc2, coda_tt;
            	DD[ggid].trms  = 0.0;
            	DD[ggid].trmscc  = 0.0;
           	DD[ggid].trms_n  = 0.0;
            	DD[ggid].tmean1 = 0.0;
            	DD[ggid].tmean2 = 0.0;
            	DD[ggid].abtrms1 = 0.0;
            	DD[ggid].abtrms2 = 0.0;
            	Nst0 = 0;
            	Nst0c = 0;
		Nst1 = 0;
		Nst2 = 0;


            	for(ist = 0; ist < Nst; ist++){ 
		   if(IV[ist].gid == ggid){
			if (IV[ist].ev1 != evf && IV[ist].ev2 != evf){
				return 0;
				break;
			}
			coda_tt = IV[ist].coda;
                	if(IV[ist].rel >= 1 ){
				ev1 = IV[ist].ev1;
				for(iid=0;iid<Est;iid++){
				   if (EIV[iid].evid == ev1){
                    			lat1 =  FEIV[iid].evla[0];
                    			lon1 =  FEIV[iid].evlo[0];
					depth1 = FEIV[iid].evdp[0];

					  dep1 = EIV[iid].evdp;
   					  lat0 = EIV[iid].evla;
					  lon0 = EIV[iid].evlo;
					break;
				   }
				}
				ev2 = IV[ist].ev2;
				for(i=0;i<Est;i++){
				   if (EIV[i].evid == ev2){
                    			lat3 =  FEIV[i].evla[0];
                    			lon3 =  FEIV[i].evlo[0];
					depth2 = FEIV[i].evdp[0];

					break;
				   }
				}
	     		        hdis = fabs(depth1-depth2);
				if (depth1 <= depth2)
					tdep = depth1;
				else
					tdep = depth2;
				nnl = sizeof(vel);
				k = 0;
				for(j=0;j<nnl;j++){
					if (dep[j] <= tdep)
						continue;
					else{
						if (k == 0){
							ndep[k] = 0;
							nvel[k] = vel[j-1];
							k++;
							ndep[k] = dep[j]-tdep;
							nvel[k] = vel[j];	
							
						}else{
							ndep[k] = dep[j]-tdep;
							nvel[k] = vel[j];
						}

						k++;
					}
				}
			       break;
                	}
		    }
		}
		
            	for(ist = 0; ist < Nst; ist++){ 
		   if(IV[ist].gid == ggid){
			if (IV[ist].ev1 != evf && IV[ist].ev2 != evf){
				return 0;
				break;
			}
			coda_tt = IV[ist].coda;
                	if(IV[ist].rel >= 1 ){
			    if (IV[ist].gid == ggid) {
 				lat2 = IV[ist].stla;
 				lon2 = IV[ist].stlo;
				ev1 = IV[ist].ev1;
				laa = floor((lat1+lat2)/2 *10);
				GCarc1 = sqrt(((lat1 -lat2)*111)*((lat1 -lat2)*111) + ((lon1 -lon2)*dis_la[laa])*((lon1 -lon2)*dis_la[laa]))/111;
				laa = floor((lat0+lat2)/2 *10);
				GCarc2 = sqrt(((lat0 -lat2)*111)*((lat0 -lat2)*111) + ((lon0 -lon2)*dis_la[laa])*((lon0 -lon2)*dis_la[laa]))/111;
				laa = floor((lat3+lat2)/2 *10);
				GCarc = sqrt(((lat3 -lat2)*111)*((lat3 -lat2)*111) + ((lon3 -lon2)*dis_la[laa])*((lon3 -lon2)*dis_la[laa]))/111;
				gca = GCarc1*111;
				depi = ((int) depth1) + 1;
				gcai = ((int) gca) + 1;
			
			if(IV[ist].phase[0] == 'P'){
				a1 = DT[0].dtdd[depi-1][gcai-1]+(DT[0].dtdd[depi-1][gcai]-DT[0].dtdd[depi-1][gcai-1])*(gca-gcai+1);
				a2 = DT[0].dtdd[depi][gcai-1]+(DT[0].dtdd[depi][gcai]-DT[0].dtdd[depi][gcai-1])*(gca-gcai+1);
				dtdd = a1 + (a2 - a1)*(depth1 - depi +1);
				a1 = DT[0].dtdh[depi-1][gcai-1]+(DT[0].dtdh[depi-1][gcai]-DT[0].dtdh[depi-1][gcai-1])*(gca-gcai+1);
				a2 = DT[0].dtdh[depi][gcai-1]+(DT[0].dtdh[depi][gcai]-DT[0].dtdh[depi][gcai-1])*(gca-gcai+1);
				dtdh = a1 + (a2 - a1)*(depth1 - depi +1);
				a1 = DT[IV[ist].mod_id].tt[depi-1][gcai-1]+(DT[IV[ist].mod_id].tt[depi-1][gcai]-DT[IV[ist].mod_id].tt[depi-1][gcai-1])*(gca-gcai+1);
				a2 = DT[IV[ist].mod_id].tt[depi][gcai-1]+(DT[IV[ist].mod_id].tt[depi][gcai]-DT[IV[ist].mod_id].tt[depi][gcai-1])*(gca-gcai+1);
				tt1[ist] = a1 + (a2 - a1)*(depth1 - depi +1);
				
				SR[ist].DD_DT   =  ((GCarc - GCarc1) * dtdd - (FEIV[i].evdp[0]-FEIV[iid].evdp[0])*dtdh);

				gca = GCarc2*111;
				depi = ((int) dep1) + 1;
				gcai = ((int) gca) + 1;

				a1 = DT[0].dtdd[depi-1][gcai-1]+(DT[0].dtdd[depi-1][gcai]-DT[0].dtdd[depi-1][gcai-1])*(gca-gcai+1);
				a2 = DT[0].dtdd[depi][gcai-1]+(DT[0].dtdd[depi][gcai]-DT[0].dtdd[depi][gcai-1])*(gca-gcai+1);
				dtdd = a1 + (a2 - a1)*(depth1 - depi +1);
				a1 = DT[0].dtdh[depi-1][gcai-1]+(DT[0].dtdh[depi-1][gcai]-DT[0].dtdh[depi-1][gcai-1])*(gca-gcai+1);
				a2 = DT[0].dtdh[depi][gcai-1]+(DT[0].dtdh[depi][gcai]-DT[0].dtdh[depi][gcai-1])*(gca-gcai+1);
				dtdh = a1 + (a2 - a1)*(depth1 - depi +1);
				gca = GCarc*111;
				depi = ((int) depth2) + 1;
				gcai = ((int) gca) + 1;
				a1 = DT[IV[ist].mod_id].tt[depi-1][gcai-1]+(DT[IV[ist].mod_id].tt[depi-1][gcai]-DT[IV[ist].mod_id].tt[depi-1][gcai-1])*(gca-gcai+1);
				a2 = DT[IV[ist].mod_id].tt[depi][gcai-1]+(DT[IV[ist].mod_id].tt[depi][gcai]-DT[IV[ist].mod_id].tt[depi][gcai-1])*(gca-gcai+1);
				tt2[ist] = a1 + (a2 - a1)*(depth2 - depi +1);
				SR[ist].DD_DT2   =  ((GCarc1 - GCarc2) * dtdd - (FEIV[iid].evdp[0]-EIV[iid].evdp)*dtdh);
                    		SR[ist].DD_DT  *=  -1.0;
                    		SR[ist].DD_DT2  *=  -1.0;
			}
			else if (IV[ist].phase[0] == 'S'){
				a1 = DT[2].dtdd[depi-1][gcai-1]+(DT[2].dtdd[depi-1][gcai]-DT[2].dtdd[depi-1][gcai-1])*(gca-gcai+1);
				a2 = DT[2].dtdd[depi][gcai-1]+(DT[2].dtdd[depi][gcai]-DT[2].dtdd[depi][gcai-1])*(gca-gcai+1);
				dtdd = a1 + (a2 - a1)*(depth1 - depi +1);
				a1 = DT[2].dtdh[depi-1][gcai-1]+(DT[2].dtdh[depi-1][gcai]-DT[2].dtdh[depi-1][gcai-1])*(gca-gcai+1);
				a2 = DT[2].dtdh[depi][gcai-1]+(DT[2].dtdh[depi][gcai]-DT[2].dtdh[depi][gcai-1])*(gca-gcai+1);
				dtdh = a1 + (a2 - a1)*(depth1 - depi +1);
				a1 = DT[2].tt[depi-1][gcai-1]+(DT[2].tt[depi-1][gcai]-DT[2].tt[depi-1][gcai-1])*(gca-gcai+1);
				a2 = DT[2].tt[depi][gcai-1]+(DT[2].tt[depi][gcai]-DT[2].tt[depi][gcai-1])*(gca-gcai+1);
				tt1[ist] = a1 + (a2 - a1)*(depth1 - depi +1);
				
				SR[ist].DD_DT   =  ((GCarc - GCarc1) * dtdd - (FEIV[i].evdp[0]-FEIV[iid].evdp[0])*dtdh);

				gca = GCarc2*111;
				depi = ((int) dep1) + 1;
				gcai = ((int) gca) + 1;

				a1 = DT[2].dtdd[depi-1][gcai-1]+(DT[2].dtdd[depi-1][gcai]-DT[2].dtdd[depi-1][gcai-1])*(gca-gcai+1);
				a2 = DT[2].dtdd[depi][gcai-1]+(DT[2].dtdd[depi][gcai]-DT[2].dtdd[depi][gcai-1])*(gca-gcai+1);
				dtdd = a1 + (a2 - a1)*(depth1 - depi +1);
				a1 = DT[2].dtdh[depi-1][gcai-1]+(DT[2].dtdh[depi-1][gcai]-DT[2].dtdh[depi-1][gcai-1])*(gca-gcai+1);
				a2 = DT[2].dtdh[depi][gcai-1]+(DT[2].dtdh[depi][gcai]-DT[2].dtdh[depi][gcai-1])*(gca-gcai+1);
				dtdh = a1 + (a2 - a1)*(depth1 - depi +1);
				gca = GCarc*111;
				depi = ((int) depth2) + 1;
				gcai = ((int) gca) + 1;
				a1 = DT[2].tt[depi-1][gcai-1]+(DT[2].tt[depi-1][gcai]-DT[2].tt[depi-1][gcai-1])*(gca-gcai+1);
				a2 = DT[2].tt[depi][gcai-1]+(DT[2].tt[depi][gcai]-DT[2].tt[depi][gcai-1])*(gca-gcai+1);
				tt2[ist] = a1 + (a2 - a1)*(depth2 - depi +1);
				SR[ist].DD_DT2   =  ((GCarc1 - GCarc2) * dtdd - (FEIV[iid].evdp[0]-EIV[iid].evdp)*dtdh);
                    		SR[ist].DD_DT  *=  -1.0;
                    		SR[ist].DD_DT2  *=  -1.0;


			}

			    }
                	}
		    }
		}
		
               for(ist = 0; ist < Nst; ist++){ 
                  if(IV[ist].rel >= 0 && IV[ist].gid == ggid){

			if ((IV[ist].DT2 > 0)&&(IV[ist].phase[0] == 'P')) {
			   SR[ist].ev2_r = IV[ist].DT2 - tt2[ist];
			   DD[ggid].abtrms2    += IV[ist].rel *SR[ist].ev2_r;
                           Nst2 += IV[ist].rel;
			}
			if ((IV[ist].DT1 > 0)&&(IV[ist].phase[0] == 'P')) {
			   SR[ist].ev1_r = IV[ist].DT1 - tt1[ist];
			   DD[ggid].abtrms1    += IV[ist].rel *SR[ist].ev1_r;
                           Nst1 += IV[ist].rel;
			}
                  }
              }
		if (Nst2 > 2){
			DD[ggid].tmean2 = DD[ggid].abtrms2/Nst2;
			FEIV[i].deot[0] = DD[ggid].tmean2;
			Nst2 = 0;
			DD[ggid].abtrms2 = 0;
		}
		if (Nst1 > 2){
			DD[ggid].tmean1 = DD[ggid].abtrms1/Nst1;
			FEIV[iid].deot[0] = DD[ggid].tmean1;
			Nst1 = 0;
			DD[ggid].abtrms1 = 0;
		}
               for(ist = 0; ist < Nst; ist++){ 
                  if(IV[ist].rel >= 0 && IV[ist].gid == ggid){
		  	if (IV[ist].TT > -10) {
				SR[ist].s_rmsc = SR[ist].DD_DT + (IV[ist].TT + DD[ggid].tmean1-DD[ggid].tmean2);
			   if(fabs(SR[ist].s_rmsc)<=0.4){
                   	        DD[ggid].trmscc += IV[ist].rel*10 * (SR[ist].s_rmsc * SR[ist].s_rmsc); 
                    		   Nst0c += IV[ist].rel*10;
			   }
			   else if((fabs(SR[ist].s_rmsc)>0.4)&&(fabs(SR[ist].s_rmsc)<=0.8)){
                   	        DD[ggid].trmscc += IV[ist].rel*3* (SR[ist].s_rmsc * SR[ist].s_rmsc); 
                    		   Nst0c += IV[ist].rel*3;
			   }
			   else if(fabs(SR[ist].s_rmsc)>0.8){
                   	        DD[ggid].trmscc += IV[ist].rel * (SR[ist].s_rmsc * SR[ist].s_rmsc); 
                    		   Nst0c += IV[ist].rel;
			   }
		  	}
			if ((IV[ist].DT2 > 0)&&(IV[ist].DT1 > 0)){
				SR[ist].s_rms = SR[ist].DD_DT + ((IV[ist].DT2 - IV[ist].DT1) + DD[ggid].tmean1-DD[ggid].tmean2);
			   if(fabs(SR[ist].s_rms)<=0.4){
                   	        DD[ggid].trms += IV[ist].rel*10 * (SR[ist].s_rms * SR[ist].s_rms); 
                    		   Nst0 += IV[ist].rel*10;
			   }
			   else if((fabs(SR[ist].s_rms)>0.4)&&(fabs(SR[ist].s_rms)<=0.8)){
                   	        DD[ggid].trms += IV[ist].rel*3 * (SR[ist].s_rms * SR[ist].s_rms); 
                    		   Nst0 += IV[ist].rel*3;
			   }
			   else if(fabs(SR[ist].s_rms)>0.8){
                   	        DD[ggid].trms += IV[ist].rel * (SR[ist].s_rms * SR[ist].s_rms); 
                    		   Nst0 += IV[ist].rel;
			   }
			}
			if (IV[ist].DT2 > 0) {
			   SR[ist].ev2_r = IV[ist].DT2 - tt2[ist] - DD[ggid].tmean2;
			   if(fabs(SR[ist].ev2_r)<=0.4){
			  	 DD[ggid].abtrms2 += IV[ist].rel*10 *SR[ist].ev2_r*SR[ist].ev2_r;
                         	  Nst2 += IV[ist].rel*10;
			   }
			   else if((fabs(SR[ist].ev2_r)>0.4)&&(fabs(SR[ist].ev2_r)<=0.8)){
			  	 DD[ggid].abtrms2 += IV[ist].rel*3 *SR[ist].ev2_r*SR[ist].ev2_r;
                         	  Nst2 += IV[ist].rel*3;
			   }
			   else if(fabs(SR[ist].ev2_r)>0.8){
			  	 DD[ggid].abtrms2 += IV[ist].rel *SR[ist].ev2_r*SR[ist].ev2_r;
                         	  Nst2 += IV[ist].rel;
			   }
				
			}
			if (IV[ist].DT1 > 0) {
			   SR[ist].ev1_r = IV[ist].DT1 - tt1[ist] - DD[ggid].tmean1;
			   if(fabs(SR[ist].ev1_r)<=0.4){
			   	DD[ggid].abtrms1    += IV[ist].rel*10 *SR[ist].ev1_r*SR[ist].ev1_r;
	                           Nst1 += IV[ist].rel*10;
			   }
			   if((fabs(SR[ist].ev1_r)>0.4)&&(fabs(SR[ist].ev1_r)<=0.8)){
			   	DD[ggid].abtrms1    += IV[ist].rel*3 *SR[ist].ev1_r*SR[ist].ev1_r;
	                           Nst1 += IV[ist].rel*3;
			   }
			   if(fabs(SR[ist].ev1_r)>0.8){
			   	DD[ggid].abtrms1    += IV[ist].rel *SR[ist].ev1_r*SR[ist].ev1_r;
	                           Nst1 += IV[ist].rel;
			   }
			}
                  }
              }
		if (Nst0 < 1)
			Nst0 = 1;
		if (Nst0c < 1)
			Nst0c = 1;
		if (Nst1 < 1)
			Nst1 = 1;
		if (Nst2 < 1)
			Nst2 = 1;

		if (Nst0 >= 3)
              		DD[ggid].trms = sqrt(DD[ggid].trms/Nst0);
		else
              		DD[ggid].trms = 0;
		if (Nst0c >= 3)
	              DD[ggid].trmscc = sqrt(DD[ggid].trmscc/Nst0c);
		else
			DD[ggid].trmscc = 0;
		if (Nst1 >= 3)
	              DD[ggid].abtrms1 = sqrt(DD[ggid].abtrms1/Nst1);
		else
			DD[ggid].abtrms1 = 0;
		if (Nst2 >= 3)
	              DD[ggid].abtrms2 = sqrt(DD[ggid].abtrms2/Nst2);
		else
			DD[ggid].abtrms2 = 0;
				
	      GCarc = sqrt(((lat3 -lat1)*111)*((lat3 -lat1)*111) + ((lon3 -lon1)*dis_la[laa])*((lon3 -lon1)*dis_la[laa]))/111;
		if (fabs(coda_tt) < 0.00001) {
       	            DD[ggid].trms_n = coef1*(DD[ggid].trms*0.3 + DD[ggid].trmscc*0.7) + coef2 *(DD[ggid].abtrms1+DD[ggid].abtrms2);
		}
		else {
    		    ttime_(&GCarc, &hdis, &k, nvel, ndep, &scoda_tt, &ain);
		    DD[ggid].coda_res = fabs(fabs(coda_tt+DD[ggid].tmean1-DD[ggid].tmean2) - scoda_tt);
       	            DD[ggid].trms_n = coef1*(DD[ggid].trms*0.3 + DD[ggid].trmscc*0.7) + coef2*(DD[ggid].abtrms1+DD[ggid].abtrms2) + coef3*DD[ggid].coda_res;
		}
}

int Calculate_Res1 (int ggid, int Nst, int Est, int evf, double coef1, double coef2, double coef3, double *dis_la, IVAL *IV, EIVAL *EIV, FEIVAL *FEIV, DVAL *DD, SRES *SR) {
	int Nst0, ist, i, iid, ev1, ev2, depi, gcai;
	int Nst1, Nst2, Nst0c,laa, cottn =0 ;
	double rms1, rmsc1, ab1, ab2;
	double lat0, lon0, lat1, lon1, lat2, lon2, lat3, lon3;
	double gca, a1, a2, dep1, depth1, depth2, dtdd, dtdh;
	double tt1[Nst], tt2[Nst], GCarc, GCarc1, GCarc2, coda_tt;

		rms1 = DD[ggid].trms;
		rmsc1 = DD[ggid].trmscc;
		ab1 = DD[ggid].abtrms1;
		ab2 = DD[ggid].abtrms2;
            	DD[ggid].trms  = 0.0;
            	DD[ggid].trmscc  = 0.0;
           	DD[ggid].trms_n  = 0.0;
            	DD[ggid].abtrms1 = 0.0;
            	DD[ggid].abtrms2 = 0.0;
            	Nst0 = 0;
            	Nst0c = 0;
		Nst1 = 0;
		Nst2 = 0;

               for(ist = 0; ist < Nst; ist++){ 
                  if(IV[ist].rel >= 0 && IV[ist].gid == ggid){
			coda_tt = IV[ist].coda;
		  	if ((IV[ist].TT > -10) && (fabs(SR[ist].s_rmsc) < 2.0*rmsc1)){
                   	        DD[ggid].trmscc += IV[ist].rel * (SR[ist].s_rmsc * SR[ist].s_rmsc); 
                    		   Nst0c += IV[ist].rel;
		  	}

			if ((IV[ist].DT2 > 0)&&(IV[ist].DT1 > 0)&& (fabs(SR[ist].s_rms) < 2.0*rms1)){
                   	        DD[ggid].trms += IV[ist].rel * (SR[ist].s_rms * SR[ist].s_rms); 
                   		   Nst0 += IV[ist].rel;
			}

			if ((IV[ist].DT2 > 0) && (fabs(SR[ist].ev2_r) < 2.0*ab2)){
			   DD[ggid].abtrms2    += IV[ist].rel *SR[ist].ev2_r*SR[ist].ev2_r;
                           Nst2 += IV[ist].rel;
			}
			if ((IV[ist].DT1 > 0) && (fabs(SR[ist].ev1_r) < 2.0*ab1)){
			   DD[ggid].abtrms1    += IV[ist].rel *SR[ist].ev1_r*SR[ist].ev1_r;
                           Nst1 += IV[ist].rel;
			}
                  }
              }
		if (Nst0 < 1)
			Nst0 = 1;
		if (Nst0c < 1)
			Nst0c = 1;
		if (Nst1 < 1)
			Nst1 = 1;
		if (Nst2 < 1)
			Nst2 = 1;
		if (Nst0 >= 3)
              		DD[ggid].trms = sqrt(DD[ggid].trms/Nst0);
		else
              		DD[ggid].trms = 0;
		if (Nst0c >= 3)
	              DD[ggid].trmscc = sqrt(DD[ggid].trmscc/Nst0c);
		else
			DD[ggid].trmscc = 0;
		if (Nst1 >= 3)
	              DD[ggid].abtrms1 = sqrt(DD[ggid].abtrms1/Nst1);
		else
			DD[ggid].abtrms1 = 0;
		if (Nst2 >= 3)
	              DD[ggid].abtrms2 = sqrt(DD[ggid].abtrms2/Nst2);
		else
			DD[ggid].abtrms2 = 0;
				
		laa = floor((lat3+lat1)/2 *10);

	      GCarc = sqrt(((lat3 -lat1)*111)*((lat3 -lat1)*111) + ((lon3 -lon1)*dis_la[laa])*((lon3 -lon1)*dis_la[laa]))/111;
		if (fabs(coda_tt) < 0.00001) {
       	            DD[ggid].trms_n = coef1*(DD[ggid].trms*0.3 + DD[ggid].trmscc*0.7) + coef2 *(DD[ggid].abtrms1+DD[ggid].abtrms2);
		}
		else {
       	            DD[ggid].trms_n = coef1*(DD[ggid].trms*0.3 + DD[ggid].trmscc*0.7) + coef2*(DD[ggid].abtrms1+DD[ggid].abtrms2) + coef3*DD[ggid].coda_res;
		}
}

int Calculate_Res2 (int ggid, int Nst, int Est, double coef1, double coef2, double coef3, double *dis_la, MDT *DT, IVAL *IV, EIVAL *EIV, FEIVAL *FEIV, DVAL *DD, SRES *SR) {
	int Nst0, ist, i, iid, ev1, ev2, depi, gcai;
	int Nst1, Nst2, Nst0c, laa;
	double rms1, rmsc1, ab1, ab2;
	double lat0, lon0, lat1, lon1, lat2, lon2, lat3, lon3;
	double gca, a1, a2, dep1, depth1, depth2, dtdd, dtdh;
	double tt1[Nst], tt2[Nst], GCarc, GCarc1, GCarc2, coda_tt;
            	DD[ggid].trms  = 0.0;
            	DD[ggid].trmscc  = 0.0;
           	DD[ggid].trms_n  = 0.0;
            	DD[ggid].tmean1 = 0.0;
            	DD[ggid].tmean2 = 0.0;
            	DD[ggid].abtrms1 = 0.0;
            	DD[ggid].abtrms2 = 0.0;
            	Nst0 = 0;
            	Nst0c = 0;
		Nst1 = 0;
		Nst2 = 0;
            	for(ist = 0; ist < Nst; ist++){ 
		   if(IV[ist].gid == ggid){
			coda_tt = IV[ist].coda;
                	if((IV[ist].rel >= 1 ) && (IV[ist].TT > -10)){
			    if (IV[ist].gid == ggid) {
 				lat2 = IV[ist].stla;
 				lon2 = IV[ist].stlo;
				ev1 = IV[ist].ev1;
				for(iid=0;iid<Est;iid++){
				   if (EIV[iid].evid == ev1){
                    			lat1 =  FEIV[iid].evla[0];
                    			lon1 =  FEIV[iid].evlo[0];
					depth1 = FEIV[iid].evdp[0];
					  dep1 = EIV[iid].evdp;
   					  lat0 = EIV[iid].evla;
					  lon0 = EIV[iid].evlo;
					break;
				   }
				}
//				GCarc1 = dist(lat1,lon1,lat2,lon2);
				laa = floor((lat1+lat2)/2 *10);
				GCarc1 = sqrt(((lat1 -lat2)*111)*((lat1 -lat2)*111) + ((lon1 -lon2)*dis_la[laa])*((lon1 -lon2)*dis_la[laa]))/111;
				ev2 = IV[ist].ev2;
				for(i=0;i<Est;i++){
				   if (EIV[i].evid == ev2){
                    			lat3 =  FEIV[i].evla[0];
                    			lon3 =  FEIV[i].evlo[0];
					depth2 = FEIV[i].evdp[0];

					break;
				   }
				}
				laa = floor((lat0+lat2)/2 *10);
				GCarc2 = sqrt(((lat0 -lat2)*111)*((lat0 -lat2)*111) + ((lon0 -lon2)*dis_la[laa])*((lon0 -lon2)*dis_la[laa]))/111;
				laa = floor((lat3+lat2)/2 *10);
				GCarc = sqrt(((lat3 -lat2)*111)*((lat3 -lat2)*111) + ((lon3 -lon2)*dis_la[laa])*((lon3 -lon2)*dis_la[laa]))/111;

				   gca = GCarc1*111;
				   depi = ((int) depth1) + 1;
				   gcai = ((int) gca) + 1;
				if(IV[ist].phase[0] == 'P'){
				   a1 = DT[0].dtdd[depi-1][gcai-1]+(DT[0].dtdd[depi-1][gcai]-DT[0].dtdd[depi-1][gcai-1])*(gca-gcai+1);
				   a2 = DT[0].dtdd[depi][gcai-1]+(DT[0].dtdd[depi][gcai]-DT[0].dtdd[depi][gcai-1])*(gca-gcai+1);
				   dtdd = a1 + (a2 - a1)*(depth1 - depi +1);
				   a1 = DT[0].dtdh[depi-1][gcai-1]+(DT[0].dtdh[depi-1][gcai]-DT[0].dtdh[depi-1][gcai-1])*(gca-gcai+1);
				   a2 = DT[0].dtdh[depi][gcai-1]+(DT[0].dtdh[depi][gcai]-DT[0].dtdh[depi][gcai-1])*(gca-gcai+1);
				   dtdh = a1 + (a2 - a1)*(depth1 - depi +1);
				   a1 = DT[IV[ist].mod_id].tt[depi-1][gcai-1]+(DT[IV[ist].mod_id].tt[depi-1][gcai]-DT[IV[ist].mod_id].tt[depi-1][gcai-1])*(gca-gcai+1);
				   a2 = DT[IV[ist].mod_id].tt[depi][gcai-1]+(DT[IV[ist].mod_id].tt[depi][gcai]-DT[IV[ist].mod_id].tt[depi][gcai-1])*(gca-gcai+1);
   				   tt1[ist] = a1 + (a2 - a1)*(depth1 - depi +1);
				   SR[ist].DD_DT   =  ((GCarc - GCarc1) * dtdd - (FEIV[i].evdp[0]-FEIV[iid].evdp[0])*dtdh);

				   gca = GCarc2*111;
				   depi = ((int) dep1) + 1;
				   gcai = ((int) gca) + 1;
				   a1 = DT[0].dtdd[depi-1][gcai-1]+(DT[0].dtdd[depi-1][gcai]-DT[0].dtdd[depi-1][gcai-1])*(gca-gcai+1);
				   a2 = DT[0].dtdd[depi][gcai-1]+(DT[0].dtdd[depi][gcai]-DT[0].dtdd[depi][gcai-1])*(gca-gcai+1);
				   dtdd = a1 + (a2 - a1)*(depth1 - depi +1);
				   a1 = DT[0].dtdh[depi-1][gcai-1]+(DT[0].dtdh[depi-1][gcai]-DT[0].dtdh[depi-1][gcai-1])*(gca-gcai+1);
				   a2 = DT[0].dtdh[depi][gcai-1]+(DT[0].dtdh[depi][gcai]-DT[0].dtdh[depi][gcai-1])*(gca-gcai+1);
				   dtdh = a1 + (a2 - a1)*(depth1 - depi +1);
				   gca = GCarc*111;
				   depi = ((int) depth2) + 1;
				   gcai = ((int) gca) + 1;
				   a1 = DT[IV[ist].mod_id].tt[depi-1][gcai-1]+(DT[IV[ist].mod_id].tt[depi-1][gcai]-DT[IV[ist].mod_id].tt[depi-1][gcai-1])*(gca-gcai+1);
				   a2 = DT[IV[ist].mod_id].tt[depi][gcai-1]+(DT[IV[ist].mod_id].tt[depi][gcai]-DT[IV[ist].mod_id].tt[depi][gcai-1])*(gca-gcai+1);
				   tt2[ist] = a1 + (a2 - a1)*(depth2 - depi +1);
				   SR[ist].DD_DT2   =  ((GCarc1 - GCarc2) * dtdd - (FEIV[iid].evdp[0]-EIV[iid].evdp)*dtdh);
                    		   SR[ist].DD_DT  *=  -1.0;
                    		   SR[ist].DD_DT2  *=  -1.0;
				   DD[ggid].tmean1 = FEIV[iid].deot[0];
				   DD[ggid].tmean2 = FEIV[i].deot[0];
				}
				else if(IV[ist].phase[0] == 'S'){
				   a1 = DT[2].dtdd[depi-1][gcai-1]+(DT[2].dtdd[depi-1][gcai]-DT[2].dtdd[depi-1][gcai-1])*(gca-gcai+1);
				   a2 = DT[2].dtdd[depi][gcai-1]+(DT[2].dtdd[depi][gcai]-DT[2].dtdd[depi][gcai-1])*(gca-gcai+1);
				   dtdd = a1 + (a2 - a1)*(depth1 - depi +1);
				   a1 = DT[2].dtdh[depi-1][gcai-1]+(DT[2].dtdh[depi-1][gcai]-DT[2].dtdh[depi-1][gcai-1])*(gca-gcai+1);
				   a2 = DT[2].dtdh[depi][gcai-1]+(DT[2].dtdh[depi][gcai]-DT[2].dtdh[depi][gcai-1])*(gca-gcai+1);
				   dtdh = a1 + (a2 - a1)*(depth1 - depi +1);
				   a1 = DT[2].tt[depi-1][gcai-1]+(DT[2].tt[depi-1][gcai]-DT[2].tt[depi-1][gcai-1])*(gca-gcai+1);
				   a2 = DT[2].tt[depi][gcai-1]+(DT[2].tt[depi][gcai]-DT[2].tt[depi][gcai-1])*(gca-gcai+1);
   				   tt1[ist] = a1 + (a2 - a1)*(depth1 - depi +1);
				   SR[ist].DD_DT   =  ((GCarc - GCarc1) * dtdd - (FEIV[i].evdp[0]-FEIV[iid].evdp[0])*dtdh);

				   gca = GCarc2*111;
				   depi = ((int) dep1) + 1;
				   gcai = ((int) gca) + 1;
				   a1 = DT[2].dtdd[depi-1][gcai-1]+(DT[2].dtdd[depi-1][gcai]-DT[2].dtdd[depi-1][gcai-1])*(gca-gcai+1);
				   a2 = DT[2].dtdd[depi][gcai-1]+(DT[2].dtdd[depi][gcai]-DT[2].dtdd[depi][gcai-1])*(gca-gcai+1);
				   dtdd = a1 + (a2 - a1)*(depth1 - depi +1);
				   a1 = DT[2].dtdh[depi-1][gcai-1]+(DT[2].dtdh[depi-1][gcai]-DT[2].dtdh[depi-1][gcai-1])*(gca-gcai+1);
				   a2 = DT[2].dtdh[depi][gcai-1]+(DT[2].dtdh[depi][gcai]-DT[2].dtdh[depi][gcai-1])*(gca-gcai+1);
				   dtdh = a1 + (a2 - a1)*(depth1 - depi +1);
				   gca = GCarc*111;
				   depi = ((int) depth2) + 1;
				   gcai = ((int) gca) + 1;
				   a1 = DT[2].tt[depi-1][gcai-1]+(DT[2].tt[depi-1][gcai]-DT[2].tt[depi-1][gcai-1])*(gca-gcai+1);
				   a2 = DT[2].tt[depi][gcai-1]+(DT[2].tt[depi][gcai]-DT[2].tt[depi][gcai-1])*(gca-gcai+1);
				   tt2[ist] = a1 + (a2 - a1)*(depth2 - depi +1);
				   SR[ist].DD_DT2   =  ((GCarc1 - GCarc2) * dtdd - (FEIV[iid].evdp[0]-EIV[iid].evdp)*dtdh);
                    		   SR[ist].DD_DT  *=  -1.0;
                    		   SR[ist].DD_DT2  *=  -1.0;
				   DD[ggid].tmean1 = FEIV[iid].deot[0];
				   DD[ggid].tmean2 = FEIV[i].deot[0];
				}
			    }
                	}
		    }
		}



               for(ist = 0; ist < Nst; ist++){ 
                  if(IV[ist].rel >= 0 && IV[ist].gid == ggid){
		  	if (IV[ist].TT > -10) {

				SR[ist].s_rmsc = SR[ist].DD_DT + (IV[ist].TT + DD[ggid].tmean1-DD[ggid].tmean2);
			   if(fabs(SR[ist].s_rmsc)<=0.4){
                   	        DD[ggid].trmscc += IV[ist].rel*10 * (SR[ist].s_rmsc * SR[ist].s_rmsc); 
                    		   Nst0c += IV[ist].rel*10;
			   }
			   else if((fabs(SR[ist].s_rmsc)>0.4)&&(fabs(SR[ist].s_rmsc)<=0.8)){
                   	        DD[ggid].trmscc += IV[ist].rel*3* (SR[ist].s_rmsc * SR[ist].s_rmsc); 
                    		   Nst0c += IV[ist].rel*3;
			   }
			   else if(fabs(SR[ist].s_rmsc)>0.8){
                   	        DD[ggid].trmscc += IV[ist].rel * (SR[ist].s_rmsc * SR[ist].s_rmsc); 
                    		   Nst0c += IV[ist].rel;
			   }


		  	}
			if ((IV[ist].DT2 > 0)&&(IV[ist].DT1 > 0)){
				SR[ist].s_rms = SR[ist].DD_DT + ((IV[ist].DT2 - IV[ist].DT1) + DD[ggid].tmean1-DD[ggid].tmean2);
			   if(fabs(SR[ist].s_rms)<=0.4){
                   	        DD[ggid].trms += IV[ist].rel*10 * (SR[ist].s_rms * SR[ist].s_rms); 
                    		   Nst0 += IV[ist].rel*10;
			   }
			   else if((fabs(SR[ist].s_rms)>0.4)&&(fabs(SR[ist].s_rms)<=0.8)){
                   	        DD[ggid].trms += IV[ist].rel*3 * (SR[ist].s_rms * SR[ist].s_rms); 
                    		   Nst0 += IV[ist].rel*3;
			   }
			   else if(fabs(SR[ist].s_rms)>0.8){
                   	        DD[ggid].trms += IV[ist].rel * (SR[ist].s_rms * SR[ist].s_rms); 
                    		   Nst0 += IV[ist].rel;
			   }
			}
			if (IV[ist].DT2 > 0) {
			   SR[ist].ev2_r = IV[ist].DT2 - tt2[ist] - DD[ggid].tmean2;
			   if(fabs(SR[ist].ev2_r)<=0.4){
			  	 DD[ggid].abtrms2 += IV[ist].rel*10 *SR[ist].ev2_r*SR[ist].ev2_r;
                         	  Nst2 += IV[ist].rel*10;
			   }
			   else if((fabs(SR[ist].ev2_r)>0.4)&&(fabs(SR[ist].ev2_r)<=0.8)){
			  	 DD[ggid].abtrms2 += IV[ist].rel*3 *SR[ist].ev2_r*SR[ist].ev2_r;
                         	  Nst2 += IV[ist].rel*3;
			   }
			   else if(fabs(SR[ist].ev2_r)>0.8){
			  	 DD[ggid].abtrms2 += IV[ist].rel *SR[ist].ev2_r*SR[ist].ev2_r;
                         	  Nst2 += IV[ist].rel;
			   }
			}
			if (IV[ist].DT1 > 0) {
			   SR[ist].ev1_r = IV[ist].DT1 - tt1[ist] - DD[ggid].tmean1;
			   if(fabs(SR[ist].ev1_r)<=0.4){
			   	DD[ggid].abtrms1    += IV[ist].rel*10 *SR[ist].ev1_r*SR[ist].ev1_r;
	                           Nst1 += IV[ist].rel*10;
			   }
			   if((fabs(SR[ist].ev1_r)>0.4)&&(fabs(SR[ist].ev1_r)<=0.8)){
			   	DD[ggid].abtrms1    += IV[ist].rel*3 *SR[ist].ev1_r*SR[ist].ev1_r;
	                           Nst1 += IV[ist].rel*3;
			   }
			   if(fabs(SR[ist].ev1_r)>0.8){
			   	DD[ggid].abtrms1    += IV[ist].rel *SR[ist].ev1_r*SR[ist].ev1_r;
	                           Nst1 += IV[ist].rel;
			   }
			}
                  }
              }
		if (Nst0 < 1)
			Nst0 = 1;
		if (Nst0c < 1)
			Nst0c = 1;
		if (Nst1 < 1)
			Nst1 = 1;
		if (Nst2 < 1)
			Nst2 = 1;

		if (Nst0 >= 3)
              		DD[ggid].trms = sqrt(DD[ggid].trms/Nst0);
		else
              		DD[ggid].trms = 0;
		if (Nst0c >= 3)
	              DD[ggid].trmscc = sqrt(DD[ggid].trmscc/Nst0c);
		else
			DD[ggid].trmscc = 0;
		if (Nst1 >= 3)
	              DD[ggid].abtrms1 = sqrt(DD[ggid].abtrms1/Nst1);
		else
			DD[ggid].abtrms1 = 0;
		if (Nst2 >= 3)
	              DD[ggid].abtrms2 = sqrt(DD[ggid].abtrms2/Nst2);
		else
			DD[ggid].abtrms2 = 0;
/*
	      GCarc = sqrt(((lat3 -lat1)*111)*((lat3 -lat1)*111) + ((lon3 -lon1)*dis_la[laa])*((lon3 -lon1)*dis_la[laa]))/111;
		if (fabs(coda_tt) < 0.00001) {
       	            DD[ggid].trms_n = coef1*(DD[ggid].trms*0.3 + DD[ggid].trmscc*0.7) + coef2 *(DD[ggid].abtrms1+DD[ggid].abtrms2);
		}
		else {
		    DD[ggid].coda_res = fabs(fabs(coda_tt+DD[ggid].tmean1-DD[ggid].tmean2)-scoda_tt);
       	            DD[ggid].trms_n = coef1*(DD[ggid].trms*0.3 + DD[ggid].trmscc*0.7) + coef2*(DD[ggid].abtrms1+DD[ggid].abtrms2) + coef3*DD[ggid].coda_res;
		}
*/
}

int Calculate_Res21 (int ggid, int Nst, int Est, double coef1, double coef2, double coef3, double *dis_la, float *vel, float *dep, IVAL *IV, EIVAL *EIV, FEIVAL *FEIV, DVAL *DD, SRES *SR) {
	int Nst0, ist, i, j, k, iid, ev1, ev2, depi, gcai;
	int Nst1, Nst2, Nst0c,laa, cottn =0, nnl, vf;
	double rms1, rmsc1, ab1, ab2;
	float nvel[100],ndep[100], ain, scoda_tt, hdis, GCarc;
	double lat0, lon0, lat1, lon1, lat2, lon2, lat3, lon3;
	double gca, a1, a2, dep1, depth1, depth2, dtdd, dtdh, tdep;
	double tt1[Nst], tt2[Nst], GCarc1, GCarc2, coda_tt;

		rms1 = DD[ggid].trms;
		rmsc1 = DD[ggid].trmscc;
		ab1 = DD[ggid].abtrms1;
		ab2 = DD[ggid].abtrms2;
            	DD[ggid].trms  = 0.0;
            	DD[ggid].trmscc  = 0.0;
           	DD[ggid].trms_n  = 0.0;
            	DD[ggid].abtrms1 = 0.0;
            	DD[ggid].abtrms2 = 0.0;
            	Nst0 = 0;
            	Nst0c = 0;
		Nst1 = 0;
		Nst2 = 0;
               for(ist = 0; ist < Nst; ist++){ 
                  if(IV[ist].rel >= 0 && IV[ist].gid == ggid){
				ev1 = IV[ist].ev1;
				ev2 = IV[ist].ev2;
			coda_tt = IV[ist].coda;
				for(i=0;i<Est;i++){
				   if (EIV[i].evid == ev1){
                    			lat1 =  FEIV[i].evla[0];
                    			lon1 =  FEIV[i].evlo[0];
					depth1 = FEIV[i].evdp[0];

					break;
				   }
				}
				for(i=0;i<Est;i++){
				   if (EIV[i].evid == ev2){
                    			lat3 =  FEIV[i].evla[0];
                    			lon3 =  FEIV[i].evlo[0];
					depth2 = FEIV[i].evdp[0];

					break;
				   }
				}
			break;
		  }
		}

				if (depth1 <= depth2)
					tdep = depth1;
				else
					tdep = depth2;
				nnl = sizeof(vel);
				k = 0;
				for(j=0;j<nnl;j++){
					if (dep[j] <= tdep)
						continue;
					else{
						if (k == 0){
							ndep[k] = 0;
							nvel[k] = vel[j-1];
							k++;
							ndep[k] = dep[j]-tdep;
							nvel[k] = vel[j];							
						}else{
							ndep[k] = dep[j]-tdep;
							nvel[k] = vel[j];
						}
					}
				}

               for(ist = 0; ist < Nst; ist++){ 
                  if(IV[ist].rel >= 0 && IV[ist].gid == ggid){

		  	if ((IV[ist].TT > -10) && (fabs(SR[ist].s_rmsc) < 2.0*rmsc1)){
                   	        DD[ggid].trmscc += IV[ist].rel * (SR[ist].s_rmsc * SR[ist].s_rmsc); 
                    		   Nst0c += IV[ist].rel;
		  	}

			if ((IV[ist].DT2 > 0)&&(IV[ist].DT1 > 0)&& (fabs(SR[ist].s_rms) < 2.0*rms1)){
                   	        DD[ggid].trms += IV[ist].rel * (SR[ist].s_rms * SR[ist].s_rms); 
                    		   Nst0 += IV[ist].rel;
			}

			if ((IV[ist].DT2 > 0) && (fabs(SR[ist].ev2_r) < 2.0*ab2)){
			   DD[ggid].abtrms2    += IV[ist].rel *SR[ist].ev2_r*SR[ist].ev2_r;
                           Nst2 += IV[ist].rel;
			}
			if ((IV[ist].DT1 > 0) && (fabs(SR[ist].ev1_r) < 2.0*ab1)){
			   DD[ggid].abtrms1    += IV[ist].rel *SR[ist].ev1_r*SR[ist].ev1_r;
                           Nst1 += IV[ist].rel;
			}
                  }
              }
		if (Nst0 < 1)
			Nst0 = 1;
		if (Nst0c < 1)
			Nst0c = 1;
		if (Nst1 < 1)
			Nst1 = 1;
		if (Nst2 < 1)
			Nst2 = 1;
		if (Nst0 >= 3)
              		DD[ggid].trms = sqrt(DD[ggid].trms/Nst0);
		else
              		DD[ggid].trms = 0;
		if (Nst0c >= 3)
	              DD[ggid].trmscc = sqrt(DD[ggid].trmscc/Nst0c);
		else
			DD[ggid].trmscc = 0;
		if (Nst1 >= 3)
	              DD[ggid].abtrms1 = sqrt(DD[ggid].abtrms1/Nst1);
		else
			DD[ggid].abtrms1 = 0;
		if (Nst2 >= 3)
	              DD[ggid].abtrms2 = sqrt(DD[ggid].abtrms2/Nst2);
		else
			DD[ggid].abtrms2 = 0;



		laa = floor((lat1+lat3)/2 *10);
	      GCarc = sqrt(((lat3 -lat1)*111)*((lat3 -lat1)*111) + ((lon3 -lon1)*dis_la[laa])*((lon3 -lon1)*dis_la[laa]))/111;
		if (fabs(coda_tt) < 0.00001) {
       	            DD[ggid].trms_n = coef1*(DD[ggid].trms*0.3 + DD[ggid].trmscc*0.7) + coef2*(DD[ggid].abtrms1+DD[ggid].abtrms2);
		}
		else {
		    hdis = fabs(depth1-depth2);
    		    ttime_(&GCarc, &hdis, &k, nvel, ndep, &scoda_tt, &ain);
		    DD[ggid].coda_res = fabs(fabs(coda_tt+DD[ggid].tmean1-DD[ggid].tmean2)-scoda_tt);
       	            DD[ggid].trms_n = coef1*(DD[ggid].trms*0.3 + DD[ggid].trmscc*0.7) + coef2*(DD[ggid].abtrms1+DD[ggid].abtrms2) + coef3*DD[ggid].coda_res;
		}
}


double dist(double evla, double evlo, double stla, double stlo)
{
    double delta;
    double s1, s2, c1, c2, cd, lon_dist;

    evla *= M_PI / 180.0;
    stla *= M_PI / 180.0;
    lon_dist = (evlo - stlo) * M_PI / 180.0;

    /* Calculate geocentric latitude.
     * Use GRS80 reference ellipsoid.
     * a = 6378137.0 meters (equatorial radius)
     * b = 6356752.0 meters (polar radius)
     * f = (b/a)^2 = 0.99330552180201
     */
    s1 = 0.99330552180201 * sin(evla);
    s2 = 0.99330552180201 * sin(stla);
    evla = atan2(cos(evla), s1);
    stla = atan2(cos(stla), s2);

    s1 = sin(evla);
    s2 = sin(stla);
    c1 = cos(evla);
    c2 = cos(stla);
    cd = cos(lon_dist);

    delta = acos(c1 * c2 + s1 * s2 * cd) / M_PI * 180.0;

    return delta;
}
int ReadFile(char *name, IVAL *IV)
{
        int i, lmax0, test;
        FILE *infile;

        lmax0=2000000;
        test=0;

        while((infile=fopen(name,"r"))==NULL){
            fprintf(stdout,"Can not open file in ReadFile %s\n", name);
            exit(-1);
        }

        for (i = 0; i <= lmax0; i++){
            if(fscanf(infile, "%lf %lf %d %d %lf %lf %lf %lf %d %s %s %d %d\n",&IV[i].stlo,&IV[i].stla,&IV[i].ev1,&IV[i].ev2,&IV[i].DT1,&IV[i].DT2,&IV[i].TT,&IV[i].coda,&IV[i].rel,IV[i].stname,IV[i].phase,&IV[i].mod_id,&IV[i].gid)==EOF)test=1;
            if(test==1)break;
        }
        fclose(infile);
        return i;
}
int ReadEF(char *name, EIVAL *EIV)
{
        int i, lmax0, test;
        FILE *infile;

        lmax0=10000;
        test=0;

        while((infile=fopen(name,"r"))==NULL){
            fprintf(stdout,"Can not open file in ReadFile %s\n", name);
            exit(-1);
        }

        for (i = 0; i <= lmax0; i++){
            if(fscanf(infile, "%d %f %lf %lf %lf %d %d\n",&EIV[i].date,&EIV[i].sec,&EIV[i].evla,&EIV[i].evlo,&EIV[i].evdp,&EIV[i].evid,&EIV[i].evflag)==EOF)test=1;
            if(test==1)break;
        }
        fclose(infile);
        return i;
}
int ReadFileM(char *name, MVAL *MIV)
{
        int i, lmax0, test;
        FILE *infile;

        lmax0=10000;
        test=0;

        while((infile=fopen(name,"r"))==NULL){
            fprintf(stdout,"Can not open file in ReadFile %s\n", name);
            exit(-1);
        }

        for (i = 0; i <= lmax0; i++){
            if(fscanf(infile, "%f %f\n",&MIV[i].dep,&MIV[i].vel)==EOF)test=1;
            if(test==1)break;
        }
        fclose(infile);
        return i;
}

int WriteFile2(char *name, EIVAL *EIV, FEIVAL *FEIV, int Est)
{
        int i, test;
        FILE *outfile;
        test=0;

        while((outfile=fopen(name,"w"))==NULL){
            fprintf(stdout,"Can not open file in WriteFile %s\n", name);
            exit(-1);
        }

        fprintf(outfile, "%s\n","oelo oela oedp felo fela evdp evid evflag");
        for (i = 0; i < Est; i++){
            fprintf(outfile, "%d %f %lf %lf %lf %lf %lf %lf %lf %d %d\n",EIV[i].date,EIV[i].sec,EIV[i].evlo,EIV[i].evla,EIV[i].evdp,FEIV[i].deot[1],FEIV[i].evlo[1],FEIV[i].evla[1],FEIV[i].evdp[1],EIV[i].evid,EIV[i].evflag);
        }
        fclose(outfile);
        return 1;
}


double uniform(double a,double b,long int *seed)
{

  double t;
  *seed=2045*(*seed)+1;
  *seed=*seed-(*seed/1048576)*1048576;
  t=(*seed)/1048576.0;
  t=a+(b-a)*t;
  return(t);
}

double cauchyrnd(double alpha, double beta, long int *s)
{
	double u,y;
	u = uniform(0.0,1.0,s);
	y = alpha - beta/tan(PI*u);
	
	return(y);
}

