
/* sets initial and boundary conditions */

#include <cstdio>
#include "common.h"

void read_inlet_2D(const std::string& finlet,int ig)
{
  int i;
  int inb;
  int inu;
  int Nint;
  int Nint_temp;
  int Nint_turb1;
  int Nint_turb2;
  int ndirn;
  int tdirn;
  double alpha;
  double Vin;
  double turbin;
  double tempin;
  double **inlet_data_vel;
  double **inlet_data_temp;
  double **inlet_data_turb1;
  double **inlet_data_turb2;
  char dummystr[300];

  /* processor 0 reads inlet data */
  FILE *fid = fopen(finlet.c_str(),"r");
  if (fid==NULL)
    nrerror("Inlet file needed for FIXV-R inlet option, not found");

  printf("\n*** Reading inlet-data file ... ");

  fgets(dummystr,300,fid);
  fscanf(fid,"%d",&Nint);
  fgets(dummystr,300,fid);

  inlet_data_vel=dmatrix(0,1,1,Nint);

  for (i=1; i<=Nint; i++) {
    fscanf(fid,"%lf %le",&inlet_data_vel[0][i],&inlet_data_vel[1][i]);
    fgets(dummystr,300,fid);
    if (i>1)
      if (inlet_data_vel[0][i]<=inlet_data_vel[0][i-1])
        nrerror("Incorrect ordering in 'inlet.dat' file: coordinate must be increasing");
  }

  if (turmod) {
    fscanf(fid,"%d",&Nint_turb1);
    fgets(dummystr,300,fid);
    inlet_data_turb1=dmatrix(0,1,1,Nint_turb1);
    for (i=1; i<=Nint_turb1; i++) {
      fscanf(fid,"%le %le",&inlet_data_turb1[0][i],&inlet_data_turb1[1][i]);
      fgets(dummystr,300,fid);
      if (i>1)
        if (inlet_data_turb1[0][i]<=inlet_data_turb1[0][i-1])
          nrerror("Incorrect ordering in 'inlet.dat' file: coordinate must be increasing");
    }

    fscanf(fid,"%d",&Nint_turb2);
    fgets(dummystr,300,fid);
    inlet_data_turb2=dmatrix(0,1,1,Nint_turb2);
    for (i=1; i<=Nint_turb2; i++) {
      fscanf(fid,"%le %le",&inlet_data_turb2[0][i],&inlet_data_turb2[1][i]);
      fgets(dummystr,300,fid);
      if (i>1)
        if (inlet_data_turb2[0][i]<=inlet_data_turb2[0][i-1])
          nrerror("Incorrect ordering in 'inlet.dat' file: coordinate must be increasing");
    }
  }

  if (temperature) {
    fscanf(fid,"%d",&Nint_temp);
    fgets(dummystr,300,fid);
    inlet_data_temp=dmatrix(0,1,1,Nint_temp);
    for (i=1; i<=Nint_temp; i++) {
      fscanf(fid,"%le %le",&inlet_data_temp[0][i],&inlet_data_temp[1][i]);
      fgets(dummystr,300,fid);
      if (i>1)
        if (inlet_data_temp[0][i]<=inlet_data_temp[0][i-1])
          nrerror("Incorrect ordering in 'inlet.dat' file: coordinate must be increasing");
    }
  }

  fclose(fid);
  printf("done\n");

/* Define normal and tangential directions */
  if (std::abs(BCgroup[ig].n[0])>=std::abs(BCgroup[ig].n[1]))
    ndirn=0;
  else
    ndirn=1;
  tdirn = (ndirn+1)%Ndim;


/* loop over nodes in inlet plane and interpolate values */
  for (inb=1; inb<=BCgroup[ig].nnode; inb++) {
    inu = Nobg[ig][inb].node;
    const double x = M.vv[tdirn][inu]-BCgroup[ig].invals[1];

    for (i=2; i<=Nint; i++)
      if (inlet_data_vel[0][i]>=x) break;

    i--;

    alpha = (x-inlet_data_vel[0][i])/(inlet_data_vel[0][i+1]-inlet_data_vel[0][i]);

    Vin = alpha*inlet_data_vel[1][i+1] + (1.-alpha)*inlet_data_vel[1][i];
    No_W[ndirn+1][inu] = std::abs(BCgroup[ig].invals[0])*Vin;
    No_W[tdirn+1][inu] = 0.;

    if (turmod) {
      for (i=2; i<=Nint_turb1; i++)
        if (inlet_data_turb1[0][i]>=x)
          break;
      i--;
      alpha = (x-inlet_data_turb1[0][i])/(inlet_data_turb1[0][i+1]-inlet_data_turb1[0][i]);
      turbin = alpha*inlet_data_turb1[1][i+1] + (1.-alpha)*inlet_data_turb1[1][i];
      No_W[iv_turb1][inu] = BCgroup[ig].invals[0]*BCgroup[ig].invals[0]*turbin;
      No_W[iv_turb1][inu] = std::max< double >(1.e-20,No_W[iv_turb1][inu]);

      for (i=2; i<=Nint_turb2; i++)
        if (inlet_data_turb2[0][i]>=x) break;
      i--;
      alpha = (x-inlet_data_turb2[0][i])/(inlet_data_turb2[0][i+1]-inlet_data_turb2[0][i]);
      turbin = alpha*inlet_data_turb2[1][i+1] + (1.-alpha)*inlet_data_turb2[1][i];
      No_W[iv_turb2][inu] = BCgroup[ig].invals[0]*BCgroup[ig].invals[0]*turbin;
      No_W[iv_turb2][inu] = std::max< double >(1.e-20,No_W[iv_turb2][inu]);

      /*
       * turbin = alpha*inlet_data_turb[2][i+1]  + (1.-alpha)*inlet_data_turb[2][i];
       * if (turmod/10==ITMGKE)
       *   No_W[iv_turb2][inu] = Cmu*pow(No_W[iv_turb1][inu],1.5)/turbin;
       * else if (turmod/10==ITMGKW)
       *   No_W[iv_turb2][inu] = sqrt(No_W[iv_turb1][inu])/turbin;
       */
    }

    if (temperature) {
      for (i=2; i<=Nint_temp; i++)
        if (inlet_data_temp[0][i]>=x) break;
      i--;
      alpha = (x-inlet_data_temp[0][i])/(inlet_data_temp[0][i+1]-inlet_data_temp[0][i]);
      tempin = alpha*inlet_data_temp[1][i+1] + (1.-alpha)*inlet_data_temp[1][i];
      No_W[iv_temp][inu] = tempin;
    }

  }

  free_dmatrix(inlet_data_vel,0,1,1,Nint);
  if (temperature)
    free_dmatrix(inlet_data_temp,0,1,1,Nint_temp);
  if (turmod) {
    free_dmatrix(inlet_data_turb1,0,1,1,Nint_turb1);
    free_dmatrix(inlet_data_turb2,0,1,1,Nint_turb2);
  }

  return;
}

