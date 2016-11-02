/*
*
* This file was auto-generated by MCWRAP
* https://github.com/magland/mcwrap
*
* You should not edit this file.
* You might not even want to read it.
* 
*/ 

#include "mex.h"

#include "../isocut5.h"
#include "../isocut5.h"
#include "../jisotonic.h"

//====================================================================
//====================================================================
        
int mcwrap_size(const mxArray *X,int j);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

      //mexPrintf("test A\n");
//   Check the number of inputs/outputs
      if (nlhs==0) nlhs=1;
      if (nrhs!=1)
         mexErrMsgTxt("Incorrect number of inputs"); 
      else if (nlhs>2)
         mexErrMsgTxt ("Too many outputs.");

      // mexPrintf("test A.2\n");
//    Setup the set inputs
        //N
        int input_N=(int)(mcwrap_size(prhs[1-1],2));

      //mexPrintf("test B\n");
//    Setup the inputs
        //samples
        //Check that we have the correct dimensions!
        {
            int numdims=mxGetNumberOfDimensions(prhs[1-1]);
            if (numdims!=2) {
              mexErrMsgTxt("Incorrect number of dimensions in input: samples");
            }
            const mwSize *dims2=mxGetDimensions(prhs[1-1]);
            int dims[]={ 1,input_N };
            for (long ii=0; ii<numdims; ii++) {
              if (dims[ii]!=dims2[ii]) {
                mexErrMsgTxt("Incorrect size of input: samples");
              }
            }
        }
        double *input_samples=mxGetPr(prhs[1-1]);
        
    
      //mexPrintf("test C\n");
//    Setup the outputs
        //dipscore
        double *output_dipscore;
        if (1<=nlhs) {
        if ((2<1)||(2>20)) {
          mexErrMsgTxt("Bad number of dimensions for my taste: 2"); 
        }
        {
            int dims2[]={ 1,1 };
            for (long ii=0; ii<2; ii++) {
                if ((dims2[ii]<1)||(dims2[ii]>10000000000.0)) {
                  mexErrMsgTxt ("Bad array size for my taste: 1,1"); 
                }
            }
        }
        
            mwSize dims[]={ 1,1 };
            plhs[1-1]=mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
            output_dipscore=mxGetPr(plhs[1-1]);
        }
        //cutpoint
        double *output_cutpoint;
        if (2<=nlhs) {
        if ((2<1)||(2>20)) {
          mexErrMsgTxt("Bad number of dimensions for my taste: 2"); 
        }
        {
            int dims2[]={ 1,1 };
            for (long ii=0; ii<2; ii++) {
                if ((dims2[ii]<1)||(dims2[ii]>10000000000.0)) {
                  mexErrMsgTxt ("Bad array size for my taste: 1,1"); 
                }
            }
        }
        
            mwSize dims[]={ 1,1 };
            plhs[2-1]=mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
            output_cutpoint=mxGetPr(plhs[2-1]);
        }

    
      //mexPrintf("test D\n");
//    Run the subroutine
        isocut5_mex(
        output_dipscore,
        output_cutpoint,
        input_N,
        input_samples

        );
   
      //mexPrintf("test E\n");
//    Free the inputs
        //samples

      //mexPrintf("test F\n");
//    Set the outputs
        //dipscore
        //cutpoint

      //mexPrintf("test G\n");

/**** We are done *******/
}

int mcwrap_size(const mxArray *X,int j) {
    mwSize numdims=mxGetNumberOfDimensions(X);
    if ((j<1)||(j>numdims)) return 1;
    const mwSize *dims2=mxGetDimensions(X);
    return dims2[j-1];
}

        //$pname$
//CC Scalar output not yet supported!
