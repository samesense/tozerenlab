#include "mex.h"
#include "matrix.h"

void simplegap(double scoredMatchMat[],
const double gap, int n, int m, double output_F[], double output_P[])
{
    // Standard Needleman-Wunsch algorithm
    
    double up, left, diagonal, best, pos;
    
    
    int i,j;
    

    for(i=0;i<m;i++)
    {
        output_F[i]=i*gap;
        output_P[i]=2;
    }
    
    
    for(j=0; j<n; j++)  //put initial values in
    {
        output_F[j*m]=j*gap;
        output_P[j*m]=4;
    }
    output_P[0]=1;
    

   
    for(j=1;j<n;j++) //cycle through columns
    {
         best=output_F[(j)*m];

         for(i=1; i<m; i++) //cycle through the rows
         {
             up=best+gap;
             left=output_F[i+(j-1)*m]+gap;
             diagonal=output_F[i-1+(j-1)*m]+scoredMatchMat[i-1+(j-1)*m];
             
             if (up>left)
             {
                 best=up;
                 pos=2;
             }
             else
            {
                best=left;
                pos=4;
            }
            
            if (diagonal>=best)
            {
                best=diagonal;
                output_P[i+j*m]=1;
            }
            else
            {
                output_P[i+j*m]=pos;
            }
            output_F[i+j*m]=best;
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[])
{
    
    double *gap, *output_F, *output_P;
    double *scoredMatchMat;
    int n, m;
    mwSize i;
    
    //double *F_col, *ptr_col;
    m=mxGetM(prhs[0]);
    n=mxGetN(prhs[0]);
    gap=mxGetPr(prhs[1]);
    
    plhs[0]=mxCreateDoubleMatrix(m,n,mxREAL);
    plhs[1]=mxCreateDoubleMatrix(m,n,mxREAL);
    
    output_F=mxGetPr(plhs[0]);
    output_P=mxGetPr(plhs[1]);
    
    
    scoredMatchMat=mxGetPr(prhs[0]);
    //
    
    
    simplegap(scoredMatchMat,*gap,n,m,output_F,output_P);
    
}

