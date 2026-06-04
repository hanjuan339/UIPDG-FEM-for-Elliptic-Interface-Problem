/*=================================================================
 *
 * aij.c, aij.mex: 
 *
 * The calling syntax is:
 *
 *		b = pdeerr_ccc(A,vi,vj)
 *
 * where A should be a sparse matrix, vi and vj be integer vertors.
 * b is a row vector satisfying b_m=A(vi_m,vj_m). 
 * This is a MEX-file for MATLAB.  
 *
 *=================================================================*/
/* $Revision: 1.0 $ */
#include "mex.h"

/* Input Arguments */

#define	A_IN	prhs[0]
#define	vi_IN	prhs[1]
#define vj_IN   prhs[2]


/* Output Arguments */

#define	b_OUT	plhs[0]

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
{ 
    double  *pr, *pi, *br, *bi, *vi, *vj;
    mwIndex  *ir, *jc;
    mwSize  ni, nj, l, row, col, k;
   
    /* Check for proper number of arguments */
    
    if (nrhs != 3) { 
	mexErrMsgTxt("Three input arguments required."); 
    } else if (nlhs > 1) {
	mexErrMsgTxt("Too many output arguments."); 
    } 
    
    ni = mxGetN(vi_IN)*mxGetM(vi_IN);
    nj = mxGetN(vj_IN)*mxGetM(vj_IN);
    if (ni != nj)
        mexErrMsgTxt("The lengths of the last two inputs should are the same");
    
    pr = mxGetPr(A_IN);
    pi = mxGetPi(A_IN);
    ir = mxGetIr(A_IN);
    jc = mxGetJc(A_IN);
    vi = mxGetPr(vi_IN);
    vj = mxGetPr(vj_IN);
    
    if (!mxIsComplex(A_IN)){
    /* Create a matrix for the return argument */ 
        b_OUT = mxCreateDoubleMatrix(1, ni, mxREAL); 
    
    /* Assign pointers to the various parameters */ 
        br = mxGetPr(b_OUT);
    
        for (l=0; l<ni; l++){
            row = *(vi+l);row--;
            col = *(vj+l);
            for (k=*(jc+col-1); k<*(jc+col); k++){
    	        if (*(ir+k)==row)
    	            *(br+l) = *(pr+k);
    	    } 
        }
    }else{
    /* Create a matrix for the return argument */ 
        b_OUT = mxCreateDoubleMatrix(1, ni, mxCOMPLEX); 
    
    /* Assign pointers to the various parameters */ 
        br = mxGetPr(b_OUT);
        bi = mxGetPi(b_OUT);
    
        for (l=0; l<ni; l++){
            row = *(vi+l);row--;
            col = *(vj+l);
            for (k=*(jc+col-1); k<*(jc+col); k++){
    	        if (*(ir+k)==row){
    	            *(br+l) = *(pr+k);
    	            *(bi+l) = *(pi+k);
    	        }
    	    } 
        }
    }
    return;
}


