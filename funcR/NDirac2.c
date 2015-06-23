#include<stdio.h>
#include<math.h>
#include "mex.h"
#define PI 3.1415926
#define X_BP0 prhs[0]
#define Y_BP0 prhs[1]
#define A0 prhs[2]
#define X prhs[3]
#define Y prhs[4]

void NumericalDirac2D(
double x_bp0[],double y_bp0[],
double x[],double y[],
double a0[],
double x_min,double x_max,double y_min,double y_max,
double dx,double dy,
int M,int m,int n,
double chem[])
{
   int i;
   int ix,jy,i1,i2,j1,j2;
   double tmpx,tmpy;
   double x0,y0;

  for(i=0;i<M;i++)
  {
     x0=x_bp0[i];
     y0=y_bp0[i];

     i1=(x0-x_min-2*dx)/dx+1;
     i2=(x0-x_min+2*dx)/dx;
     j1=(y0-y_min-2*dy)/dy+1;
     j2=(y0-y_min+2*dy)/dy;
     
     for(ix=i1;ix<i2+1;ix++)
     {
	     for(jy=j1;jy<j2+1;jy++)
	     {
		     tmpx=(x[ix+m*jy]-x0)*PI/(2*dx);
		     tmpy=(y[ix+m*jy]-y0)*PI/(2*dy);
		     chem[ix+m*jy]+=a0[i]*(1+cos(tmpx))*(1+cos(tmpy))/(16*dx*dy);
	     }
      }

    }

				     
}

/*int main()
{
   int M;
   int m=100,n=100;
   int i,j;
   double dx,dy,x_min,x_max,y_min,y_max;
    
   M=200;
   dx=40.0/m; 
   dy=40.0/n;

   double chem[m*n],x_bp0[M],y_bp0[M],x[m*n],y[m*n],a0[M];
   
   for(i=0;i<M;i++)
   {
	 x_bp0[i]=12*cos(2*PI*i/M);
         y_bp0[i]=12*sin(2*PI*i/M);
         a0[i]=2*cos(2*PI*i/M)*cos(2*PI*i/M);
   }

   for(i=0;i<m;i++)
   {
	   for(j=0;j<n;j++)
	   {
		   x[i+n*j]=-20.0+dx*i;
		   y[i+n*j]=-20.0+dy*j;
		   chem[i+n*j]=0;
	   }
   }

     x_min=-20.0;
     x_max=20.0-dx;
     y_min=-20.0;
     y_max=20.0-dy;


   NumericalDirac2D(x_bp0,y_bp0,x,y,a0,x_min,x_max,y_min,y_max,dx,dy,M,m,n,chem);
    
   for(i=0;i<m;i++)
   {
	   for(j=0;j<n;j++)
	   {
		   printf("%.1f ",chem[i+n*j]);
	   }
	   printf("\n");
   }



  
}     

*/    

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  double x_min,x_max,y_min,y_max;
  double dx,dy;
  double *chem;
  int m,n,M;
  int i,j;
  double *x,*y,*x_bp0,*y_bp0,*a0;

  x=mxGetPr(X);
  y=mxGetPr(Y);
  a0=mxGetPr(A0);
  x_bp0=mxGetPr(X_BP0);
  y_bp0=mxGetPr(Y_BP0);

  m=mxGetM(X);
  n=mxGetN(X);
  M=mxGetN(X_BP0);

  x_min=x[0];
  x_max=x[(n-1)*m];
  y_min=y[0];
  y_max=y[m-1];

  dx=x[n]-x[0];
  dy=y[1]-y[0];

  plhs[0]=mxCreateDoubleMatrix(m,n,mxREAL);
  chem=mxGetPr(plhs[0]); 

  for(i=0;i<m;i++)
  {
	  for(j=0;j<n;j++)
	  {
		  chem[i+m*j]=0.0;
	  }
  }

  NumericalDirac2D(x_bp0,y_bp0,x,y,a0,x_min,x_max,y_min,y_max,dx,dy,M,m,n,chem);
} 



