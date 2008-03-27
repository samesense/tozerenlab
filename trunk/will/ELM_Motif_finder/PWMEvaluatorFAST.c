#include "mex.h"
#include "matrix.h"


void PWMEval(mxArray *SeqValCell_ptr, mxArray *AA_SEQ_CELL_ptr, double PWM_ptr[], int NumSeqs, int PWM_cols, int PWM_rows)
{
    double *this_seq,*this_val, current_val;
    mwSize i,j,k,this_seq_length;
    mxArray *temp_array;
    
    for (i=0; i<NumSeqs; i++)
    {
        //get the data out of the AA_SEQ_CELL_ptr
        this_seq_length=mxGetN(mxGetCell(AA_SEQ_CELL_ptr,i));
        this_seq=mxGetPr(mxGetCell(AA_SEQ_CELL_ptr,i));

        //create a place to put the result in SeqValCell_ptr
        temp_array=mxCreateDoubleMatrix(1,this_seq_length,mxREAL);
        mxSetCell(SeqValCell_ptr,i,temp_array);

        //get the ptr to the new location
        this_val=mxGetPr(mxGetCell(SeqValCell_ptr,i));
        
        for (k=0; k<(this_seq_length-PWM_cols); k++)
        {
            current_val=0;
            for (j=0; j<PWM_cols; j++)
                current_val+=PWM_ptr[PWM_rows*j+(int)this_seq[k+j]-1];
            this_val[k]=current_val;
        }
        
    }
}


void mexFunction(   int nlhs, mxArray *plhs[],
                    int nrhs, const mxArray *prhs[])
{
    //PWMEvaluator(PWM,AA_SEQ_CELL)
    //
    //RETURNS = SeqValCell
    //
    
    //Initialize Inputs
    double *PWM_data;

    mxArray *AA_SEQ_CELL_ptr;
    mwSize NumSeqs,PWM_cols,PWM_rows;
    
    //Initialize Outputs
    mxArray *SeqValCell_ptr;
    
    //Initialize Internal Variables
    
    
    //Create Gateway for Inputs
    
    AA_SEQ_CELL_ptr=prhs[1];
    
    PWM_data=mxGetPr(prhs[0]);
    
    PWM_cols=mxGetN(prhs[0]);
    PWM_rows=mxGetM(prhs[0]);
    
    NumSeqs=mxGetM(AA_SEQ_CELL_ptr);

    
    ///////////////////////////
    
    //Create Output Matrices
    
    
    //Output Matrices are the same size as the input ones;

    SeqValCell_ptr=mxCreateCellMatrix(NumSeqs,1);
    plhs[0]=SeqValCell_ptr;
    
    
    
   // mexPrintf("finished gateway \n");
    PWMEval(SeqValCell_ptr,AA_SEQ_CELL_ptr,PWM_data,NumSeqs,PWM_cols,PWM_rows);
    

    //Create Gateway for Output Matrices
    return;
    
}
