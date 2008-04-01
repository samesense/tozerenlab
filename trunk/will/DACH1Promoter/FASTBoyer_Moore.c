#include "mex.h"
#include "matrix.h"

struct pattern_data
{
    int m;
    char *pat;
    int skip[127];
} ;




void Boyer_Moore(double STARTING_POS_ptr[], mxArray *PATTERN_ptr, mxArray *SEARCH_ptr, int  num_patterns, int num_searches)
{
    
    
    /////Matlab implementation variables
    int pattern_counter, search_counter;
    mxArray *this_search, *this_pattern;
    
    /////Boyer_Moore Variables
    //int MAXCHAR=127;
    
    
    
    char *pat, *text;
    int n, test_spot;
    int i, j, k, m;
    
    struct pattern_data *super_pat;
    
    bool *been_found;
    
    super_pat=mxCalloc(num_patterns,(2*sizeof(int)+2*sizeof(char)+127*sizeof(int)));
    
    been_found=mxCalloc(num_searches,sizeof(bool));
    
    
    for(pattern_counter=0; pattern_counter<num_patterns; pattern_counter++)
    {
        //        mexPrintf("%s%d\n","Checking pattern: ",pattern_counter);
        this_pattern=mxGetCell(PATTERN_ptr,pattern_counter);
        m=mxGetN(this_pattern);
        
        pat=mxCalloc(m,sizeof(char));

        //mexPrintf("%s%d\n","finished allocating pattern array: ",pattern_counter);
        
        mxGetString(this_pattern,pat,m+1);
        super_pat[pattern_counter].pat=pat;

        super_pat[pattern_counter].m=m;
        
        //mexPrintf("%s\n%s\n","finished putting char in array: ",super_pat[pattern_counter].pat);
        
        
        //Generate Skipping matrix
        for( k=0; k<127; k++ ) super_pat[pattern_counter].skip[k] = m;
        for( k=0; k<m-1; k++ ) super_pat[pattern_counter].skip[pat[k]] = m-k-1;

        //mexPrintf("%s\n%s\n","finished putting char in array: ",pat);
    }
    
    //mexPrintf("%s\n%i\n","finished creating skipping matrix: ",skip);
    for(search_counter=0;search_counter<num_searches;search_counter++)
    {
        //if STARTING_POS_ptr[search_counter]==0
        //continue
        
       // mexPrintf("%s%d\n","starting search: ",search_counter);
        this_search=mxGetCell(SEARCH_ptr,search_counter);
        
        n=mxGetN(this_search);
        text=mxCalloc(n,sizeof(char));
        //mexPrintf("%s%d\n","Finished allocating search array: ",search_counter);
        mxGetString(this_search,text,n+1);
        //            mexPrintf("%s\n%d\n","finished putting char search array: ",n);
        
        for(pattern_counter=0;pattern_counter<num_patterns; pattern_counter++)
        {
            for( k=super_pat[pattern_counter].m-1; k < n; k += super_pat[pattern_counter].skip[text[k] ] )//& (127-1)
            {
                for( j=super_pat[pattern_counter].m-1, i=k; j>=0 && text[i] == super_pat[pattern_counter].pat[j]; j-- )
                {
                    i--;
                    //mexPrintf("%s%d%s%d\n","Matches: ",j," poss: ",i);
                    
                }
                
                if( j == (-1) )
                {
                    
//                     mexPrintf("pattern:%d\t loc:%d\t search:%d\t putting loc:%d\n%s\n",pattern_counter+1,i+2,search_counter+1,search_counter+num_searches*pattern_counter ,pat);
//                     for(test_spot=0;test_spot<m;test_spot++)
//                         mexPrintf("%c",text[test_spot+i+1]);
//                     mexPrintf("\n");
                    STARTING_POS_ptr[search_counter+num_searches*pattern_counter]=(i+1)+1;
                    been_found[search_counter]=true;
                    break;
                }
            }
            
        }
         mxFree(text);
        //mexPrintf("%s%d\n","Finished de-allocating search array: ",search_counter);
            
    }
 
    mxFree(been_found);
    mxFree(super_pat);
}






void mexFunction(int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[])
{
    //FASTBoyer_Moore(PATTERN_STRINGS,SEARCH_STRINGS)
    //
    //RETURNS = [STARTING_POS]
    //
    //
    //
    //
    //
    
    // STARTING_POS is the UNIQUE starting position.  0 if the string is
    //      not found or is found multiple times.
    
    
    
    //Initialize Inputs
    mxArray *PATTERN_ptr;
    mxArray *SEARCH_ptr;
    
    int num_patterns, out_size;
    
    
    
    //Initialize Outputs
    double *STARTING_POS_ptr;
    
    
    //Initialize Internal Variables
    
    
    //Create Gateway for Inputs
    
    PATTERN_ptr=(prhs[0]);
    SEARCH_ptr=(prhs[1]);
    
    
    num_patterns=mxGetM(prhs[0]);
    out_size=mxGetM(prhs[1]);
    
    
    ///////////////////////////
    
    //Create Output Matrices
    
    
    //Output Matrices are the same size as the input ones;
    plhs[0]=mxCreateDoubleMatrix(out_size,num_patterns,mxREAL);
    
    STARTING_POS_ptr=mxGetPr(plhs[0]);
    //  mexPrintf("finished gateway");
    
    
    Boyer_Moore(STARTING_POS_ptr,PATTERN_ptr,SEARCH_ptr,num_patterns,out_size);
    
    //Create Gateway for Output Matrices
    
    ////////////////////////////////
}









