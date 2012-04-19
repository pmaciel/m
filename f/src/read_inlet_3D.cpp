
/* sets initial and boundary conditions */

#include <cstdio>
#include "common.h"

/*
 * read and interpolate velocity profile from inlet file
 *
 * all dependent variables must be present in file.
 *
 * data contained in xplot file "inlet.xpl" as non-uniform cartesian grid
 * first line: "unstprim grid data Nxid Nyid"
 *
 * input data from start file:
 * BCgroup[ig].invals[0]=scaling factor for velocity components
 * BCgroup[ig].invals[1]=grid direction (0/1/2) corresponding to x in inlet.xpl
 * BCgroup[ig].invals[2]=grid direction (0/1/2) corresponding to y in inlet.xpl
 * BCgroup[ig].invals[3]=offset in x direction of inlet data
 * BCgroup[ig].invals[4]=offset in y direction of inlet data
 */
void read_inlet_3D(const std::string& finlet,int ig)
{
  int i;
  int j;
  int ij;
  int ic;
  int id;
  int inb;
  int inu;
  int Nnid;
  int Nxid;
  int Nyid;
  int Ncid;
  int idirn;
  int jdirn;
  int kdirn;
  double alpha;
  double beta;
  double Us;
  double Un;
  double Vin;
  double turbn;
  double turbs;
  double turbin;
  double temps;
  double tempn;
  double tempin;
  double rdum;
  double *inlet_data_x;
  double *inlet_data_y;
  double **inlet_data_ve;
  double *inlet_data_temp  = NULL;
  double **inlet_data_turb = NULL;
  char dummystr[100];
  int   ret = 0;
  char* rec = NULL;


  std::cout << "read_inlet_3D: reading file \"" << finlet << "\"..." << std::endl;
  FILE *fid = fopen(finlet.c_str(),"r");
  if (fid==NULL)
    nrerror("Inlet file needed for FIXV-R inlet option, not found");


  ret = fscanf(fid,"unstprim grid data %d %d\n",&Nxid,&Nyid);
  ret = fscanf(fid,"%d %d",&Ncid,&Nnid);
  rec = fgets(dummystr,100,fid);

  if (Nnid!=(Nxid*Nyid))
    nrerror("In file 'inlet.xpl' number of nodes not equal to Nx*Ny");

  for (ic=0; ic<Ncid+2; ic++)
    rec = fgets(dummystr,100,fid);

  inlet_data_x=dvector(1,Nxid);
  inlet_data_y=dvector(1,Nyid);
  inlet_data_ve=dmatrix(1,3,1,Nxid*Nyid);
  if (temperature)
    inlet_data_temp=dvector(1,Nxid*Nyid);
  if (turmod)
    inlet_data_turb=dmatrix(1,2,1,Nxid*Nyid);

  ij=0;
  for (i=1; i<=Nxid; i++)
  for (j=1; j<=Nyid; j++) {
    ij++;
    ret = fscanf(fid,"%lf %lf %le %le %le %le",&inlet_data_x[i],&inlet_data_y[j],
      &inlet_data_ve[1][ij],&inlet_data_ve[2][ij],&rdum,&inlet_data_ve[3][ij]);
    if (temperature)
      ret = fscanf(fid,"%le",&inlet_data_temp[ij]);
    if (turmod)
      ret = fscanf(fid,"%le %le",&inlet_data_turb[1][ij],&inlet_data_turb[2][ij]);
    rec = fgets(dummystr,100,fid);
    if (j>1)
      if (inlet_data_y[j]<=inlet_data_y[j-1])
        nrerror("Incorrect ordering in 'inlet.xpl' file: must be structured i j (column by column)");
  }

  fclose(fid);
  std::cout << "read_inlet_3D: reading file." << std::endl;


/* find coordinate direction normal to inlet plane (x,y or z) */
  kdirn=0;
  ij=0;
  for (id=1; id<Ndim; id++)
    if (std::abs(BCgroup[ig].n[id])>=std::abs(BCgroup[ig].n[kdirn]))
      kdirn=id;

/* assign tangential directions according to input settings */
  idirn = (int)(BCgroup[ig].invals[1]+0.0001);
  jdirn = (int)(BCgroup[ig].invals[2]+0.0001);

/* loop over nodes in inlet plane and interpolate values */
  for (inb=1; inb<=BCgroup[ig].nnode; inb++) {
    inu = Nobg[ig][inb].node;
    const double x = M.vv[idirn][inu]-BCgroup[ig].invals[3];
    const double y = M.vv[jdirn][inu]-BCgroup[ig].invals[4];

    for (i=2; i<=Nxid; i++)
      if (inlet_data_x[i]>=x) break;
    for (j=2; j<=Nyid; j++)
      if (inlet_data_y[j]>=y) break;

    i--; j--; ij = (i-1)*Nyid + j;

    alpha = (x-inlet_data_x[i])/(inlet_data_x[i+1]-inlet_data_x[i]);
    beta  = (y-inlet_data_y[j])/(inlet_data_y[j+1]-inlet_data_y[j]);

    Us = alpha*inlet_data_ve[1][ij+Nyid]   + (1.-alpha)*inlet_data_ve[1][ij];
    Un = alpha*inlet_data_ve[1][ij+Nyid+1] + (1.-alpha)*inlet_data_ve[1][ij+1];
    Vin = beta*Un + (1.-beta)*Us;
    No_W[idirn+1][inu] = std::abs(BCgroup[ig].invals[0])*Vin;

    Us = alpha*inlet_data_ve[2][ij+Nyid]   + (1.-alpha)*inlet_data_ve[2][ij];
    Un = alpha*inlet_data_ve[2][ij+Nyid+1] + (1.-alpha)*inlet_data_ve[2][ij+1];
    Vin = beta*Un + (1.-beta)*Us;
    No_W[jdirn+1][inu] = std::abs(BCgroup[ig].invals[0])*Vin;

    Us = alpha*inlet_data_ve[3][ij+Nyid]   + (1.-alpha)*inlet_data_ve[3][ij];
    Un = alpha*inlet_data_ve[3][ij+Nyid+1] + (1.-alpha)*inlet_data_ve[3][ij+1];
    Vin = beta*Un + (1.-beta)*Us;
    No_W[kdirn+1][inu] = BCgroup[ig].invals[0]*Vin;

    if (temperature) {
      temps=alpha*inlet_data_temp[ij+Nyid]  +(1.-alpha)*inlet_data_temp[ij];
      tempn=alpha*inlet_data_temp[ij+Nyid+1]+(1.-alpha)*inlet_data_temp[ij+1];
      tempin = beta*tempn + (1.-beta)*temps;
      No_W[iv_temp][inu] = tempin;
    }

    if (turmod) {
      turbs=alpha*inlet_data_turb[1][ij+Nyid]  +(1.-alpha)*inlet_data_turb[1][ij];
      turbn=alpha*inlet_data_turb[1][ij+Nyid+1]+(1.-alpha)*inlet_data_turb[1][ij+1];
      turbin = beta*turbn + (1.-beta)*turbs;
      No_W[iv_turb1][inu] = BCgroup[ig].invals[0]*BCgroup[ig].invals[0]*turbin;

      turbs=alpha*inlet_data_turb[2][ij+Nyid]  +(1.-alpha)*inlet_data_turb[2][ij];
      turbn=alpha*inlet_data_turb[2][ij+Nyid+1]+(1.-alpha)*inlet_data_turb[2][ij+1];
      turbin = beta*turbn + (1.-beta)*turbs;
      No_W[iv_turb2][inu] = turbin;
    }

    /*
     * modification for non-normal inflow:
     * No_W[jdirn+1][inu] = 0.05*No_W[kdirn+1][inu];
     */
  }

  free_dvector(inlet_data_x,1,Nxid);
  free_dvector(inlet_data_y,1,Nyid);
  free_dmatrix(inlet_data_ve,1,3,1,Nxid*Nyid);
  if (turmod)
    free_dmatrix(inlet_data_turb,1,2,1,Nxid*Nyid);

  if (ret) ret=0;
  if (rec) rec=NULL;
}

