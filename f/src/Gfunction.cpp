
/* generation function */

double Gfunction(const int Ndim, const double gradv[4][3])
{
  double G = 0.;
  for (int iv=1; iv<=Ndim; ++iv)
    G += 2.*gradv[iv][iv-1]*gradv[iv][iv-1];

  for (int i=1; i<Ndim; ++i)
    for (int j=i+1; j<=Ndim; ++j)
      G += (gradv[j][i-1]+gradv[i][j-1])*(gradv[j][i-1]+gradv[i][j-1]);

  return G;
}

