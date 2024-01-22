
/**
 * OPO_mex_ode23.cpp
 * Testing implementation of functions within  
 * main simulation loop (cOPO_gpu_prism_single_mex_test_script.m)
 * on GPU?
 * Starting with nonlinear function.
 */

#include "mex.h"
#include "gpu/mxGPUArray.h"
#include <math_constants.h>
#include "OPO_ODE23_HEADER.hpp"


/**
 * MEX gateway
 */
void mexFunction(int /*nlhs*/, mxArray *plhs[],
                 int nrhs, mxArray const *prhs[])
{
	/* Declare all variables.*/
	mxGPUArray const *alpha_w, *G, *beta_op;
    mxGPUArray *A0, *A1, *NL, *K1, *K2, *K3, *K4, *E, *Ap, *stepmods;// *EFT;
	float const *d_beta_op; 
    float const  *d_alpha_w, *d_G;
    float2  *d_A0, *d_A1, *d_NL, *d_K1, *d_K2, *d_K3, *d_K4, *d_E, *d_Ap;// *d_EFT;
	float w0, h, max_pct_err, min_pct_err, *d_stepmods;
    double bdiffw0, dt;
    int Npoints, Nsteps, sel, Nbatches;
    char const * const errId = "parallel:gpu:mexGPUExample:InvalidInput";
    char const * const errMsg = "Invalid input to MEX file.";

    // Initialize the MathWorks GPU API.
    mxInitGPU();

	// We expect to receive as input 
	if (nrhs!=14) {
	        mexErrMsgIdAndTxt(errId, errMsg);
	    }

    E           = mxGPUCopyFromMxArray(prhs[0]);
    alpha_w     = mxGPUCreateFromMxArray(prhs[1]);
    G           = mxGPUCreateFromMxArray(prhs[2]);
    w0          = (float)mxGetScalar(prhs[3]);
    bdiffw0     = (double)mxGetScalar(prhs[4]);
    h           = (float)mxGetScalar(prhs[5]);
    Nsteps      = (int)mxGetScalar(prhs[6]);
    dt          = (double)mxGetScalar(prhs[7]);
    max_pct_err = (float)mxGetScalar(prhs[9]);
    min_pct_err = (float)mxGetScalar(prhs[10]);
    sel         = (int)mxGetScalar(prhs[11]);

    beta_op = mxGPUCreateFromMxArray(prhs[8]);
    Ap      = mxGPUCopyFromMxArray(prhs[12]);
    stepmods = mxGPUCopyFromMxArray(prhs[13]);

	 /*
     * Extract pointers to the input
     * data on the device.
     */

    d_E         = (float2 *)(mxGPUGetData(E));
    d_Ap        = (float2 *)(mxGPUGetData(Ap));
    d_stepmods  = (float *)(mxGPUGetData(stepmods));
    d_beta_op   = (float const *)(mxGPUGetDataReadOnly(beta_op));
    d_G         = (float const *)(mxGPUGetDataReadOnly(G));
    d_alpha_w   = (float const *)(mxGPUGetDataReadOnly(alpha_w));


    Npoints = (int)(mxGPUGetNumberOfElements(beta_op));
    Nbatches = (int)(mxGPUGetNumberOfElements(E)/Npoints);

    // Create intermediate arrays and get pointers
    /*
    A0 = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(E),
                            mxGPUGetDimensions(E),
                            mxGPUGetClassID(E),
                            mxGPUGetComplexity(E),
                            MX_GPU_DO_NOT_INITIALIZE);
    */
    A0           = mxGPUCopyFromMxArray(prhs[0]);
    d_A0 = (float2 *)(mxGPUGetData(A0));
    

    A1 = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(E),
                            mxGPUGetDimensions(E),
                            mxGPUGetClassID(E),
                            mxGPUGetComplexity(E),
                            MX_GPU_DO_NOT_INITIALIZE);
    d_A1 = (float2 *)(mxGPUGetData(A1));

    NL = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(E),
                            mxGPUGetDimensions(E),
                            mxGPUGetClassID(E),
                            mxGPUGetComplexity(E),
                            MX_GPU_DO_NOT_INITIALIZE);
    d_NL = (float2 *)(mxGPUGetData(NL));

    K1 = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(E),
                            mxGPUGetDimensions(E),
                            mxGPUGetClassID(E),
                            mxGPUGetComplexity(E),
                            MX_GPU_DO_NOT_INITIALIZE);
    d_K1 = (float2 *)(mxGPUGetData(K1));

    K2 = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(E),
                            mxGPUGetDimensions(E),
                            mxGPUGetClassID(E),
                            mxGPUGetComplexity(E),
                            MX_GPU_DO_NOT_INITIALIZE);
    d_K2 = (float2 *)(mxGPUGetData(K2));

    K3 = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(E),
                            mxGPUGetDimensions(E),
                            mxGPUGetClassID(E),
                            mxGPUGetComplexity(E),
                            MX_GPU_DO_NOT_INITIALIZE);
    d_K3 = (float2 *)(mxGPUGetData(K3));

    K4 = mxGPUCreateGPUArray(mxGPUGetNumberOfDimensions(E),
                            mxGPUGetDimensions(E),
                            mxGPUGetClassID(E),
                            mxGPUGetComplexity(E),
                            MX_GPU_DO_NOT_INITIALIZE);
    d_K4 = (float2 *)(mxGPUGetData(K4));

    // Call .cu function
 OPO_TEST_CUDA(d_E, d_A0, d_A1, d_NL, d_K1, d_K2, d_K3, d_K4, d_G, d_alpha_w,
                    w0,bdiffw0,h,Nsteps,Npoints,Nbatches,dt,d_beta_op,
                    max_pct_err, min_pct_err, sel, d_Ap, d_stepmods);

    /* Wrap the result up as a MATLAB gpuArray for return. */
    plhs[0] = mxGPUCreateMxArrayOnGPU(A0);
    plhs[1] = mxGPUCreateMxArrayOnGPU(Ap);
    plhs[2] = mxGPUCreateMxArrayOnGPU(stepmods);

    /*
     * The mxGPUArray pointers are host-side structures that refer to device
     * data. These must be destroyed before leaving the MEX function.
     */
    mxGPUDestroyGPUArray(A0);
    mxGPUDestroyGPUArray(A1);
    mxGPUDestroyGPUArray(NL);
    mxGPUDestroyGPUArray(K1);
    mxGPUDestroyGPUArray(K2);
    mxGPUDestroyGPUArray(K3);
    mxGPUDestroyGPUArray(K4);
    mxGPUDestroyGPUArray(E);
    mxGPUDestroyGPUArray(Ap);
    mxGPUDestroyGPUArray(G);
    mxGPUDestroyGPUArray(beta_op);
    mxGPUDestroyGPUArray(alpha_w);
    mxGPUDestroyGPUArray(stepmods);

    
}