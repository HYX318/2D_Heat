#include <mpi.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>
#include "HeatUtils.h"
#include "Interfaces.h"
void ReadParam(int &Nx, int &Ny, int &Nt, double &StabP){
  std::ifstream myfile;
  myfile.open("Param.in");
  if(!myfile) { // file couldn't be opened
    std::cerr << "Param.in file could not be opened" << std::endl;
  }
  myfile >> Nx >> Ny ;
  myfile >> Nt;
  myfile >> StabP ; 
 }

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////   2D Contigous Array  ///////////////////////////////////////////
/////                                         Allocation                                          ///////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

// 2D contiguous array memory allocation
void Contiguous2D(double  **Tab, int rowLength, int colLength){
  Tab[0]=new double[colLength*rowLength];
  for (int i=1; i < colLength; i++)
    Tab[i]= Tab[0] + i * rowLength ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////     Write to File    ///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

void Export(double **Sol, int I[], int J[], double x[], double y[], int Nx, int Ny, std::string SolFile){
	std::ofstream outfile;
	outfile.open(SolFile);
   	for (int j = 1; j < Ny+1; j++)
   	{
       for (int i = 1; i < Nx+1; i++){        
        outfile << I[i] << " " << J[j] << " " << x[i]<< " " << y[j] << " " << std::setprecision(15) <<  Sol[j][i] << std::endl;
       //	outfile << Sol[j][i] << " " ;
       }
       //outfile << std::endl;
   }
   outfile.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////          Init        ///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

void Init(double **Sol, double x[], double y[], int Nx, int Ny){
	double t0=0.;
	ExactSol(Sol,x,y,t0,Nx,Ny);
}

double Analytic(double x, double y, double t){
  double pi = atan(1)*4;
  double k0 = -5*pi*pi;
  double Val = sin(pi*x)*sin(2*pi*y)*exp(t*k0);
  return Val;
}

void ExactSol(double **Sol, double x[], double y[], double t, int Nx, int Ny){
for (int j=0;j<Ny+2;j++) 
    for (int i=0;i<Nx+2;i++){
      Sol[j][i] = Analytic(x[i],y[j],t);
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////   Error  /////////////////////////////////////////////////////
/////                                        Norms                                                ///////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

double TwoNorm(double **Err, int Nx, int Ny, double h){
  double e_2;
  e_2=0.;
  for (int j=1;j<Ny+1;j++) 
      for (int i=1;i<Nx+1;i++)
        e_2 = e_2 + pow(Err[j][i],2);
  return sqrt(e_2)*h;
}

double InftyNorm(double **Err, int Nx, int Ny){
  double e_infty;
  e_infty = 0.;
  for (int j=1;j<Ny+1;j++) 
      for (int i=1;i<Nx+1;i++)
        e_infty = fmax(e_infty , fabs(Err[j][i]));
      return e_infty;
    }

void Error(double **Sol, double **Ref, double **Err, int Nx, int Ny){
  for (int j=0;j<Ny+2;j++) 
      for (int i=0;i<Nx+2;i++)
        Err[j][i] = Sol[j][i] - Ref[j][i];
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////  Implicit Scheme  ////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

void Copy(double **Sol, double **Sol0, int Nx, int Ny){
  for (int j=0;j<Ny+2;j++) 
      for (int i=0;i<Nx+2;i++)
        Sol0[j][i] = Sol[j][i];
}

void Jacobi(double **x, double **x0, double **b, double &Residu, double Tol, 
            int &iConv, int Nx, int Ny, double lambda){
  double Coeff;
  Coeff = double(1)/(1+4*lambda);
  double Res;
  Residu = 1.0;
  iConv = 0;
  while (Residu > Tol ){
    iConv = iConv + 1;
    Res= 0.0;
    for (int j=1; j<Ny+1; j++)
      for (int i=1; i<Nx+1; i++){
        x[j][i] = Coeff* ( lambda* (x0[j+1][i]+x0[j][i+1]
                                   +x0[j-1][i]+x0[j][i-1]) + b[j][i] );
        Res = fmax(abs(x[j][i] - x0[j][i]),Res);
    }
    Residu = Res;  
    Copy(x,x0,Nx,Ny);
  }
}
