/* 
 *OPO_CUDA_ODE23.cu
 *
 */

#include <cuda_runtime_api.h> 
#include <math_constants.h>
#include <cufft.h>
#include "OPO_ODE23_HEADER.hpp"
#include <cuda_runtime.h>
#include "device_launch_parameters.h"
#include <limits>

/*Define block size.*/
#define threadsPerBlock 256

/* Define imaginary constant */
#define I make_float2(0.0,1.0)
#define c0 make_float2(0.0,0.0)
//const float2 I = make_float2(0.0,1.0);
//const float2 c0 = make_float2(0.0,0.0);

__constant__ float w;
__constant__ float d_invdt2;  
__constant__ int Np = 1; 

// complex math functions
__device__ inline
float2 conjugate(float2 arg)
{
    return make_float2(arg.x, -arg.y);
}

__device__ inline
float2 complex_exp(float arg)
{
    return make_float2(cosf(arg), sinf(arg));
}

__device__ inline
float2 complex_pow(float2 arg, int n)
{
	// return make_float2(cosf(acosf(arg.x)*n), sinf(asinf(arg.y)*n));
	return make_float2(cosf(acosf(arg.x) * n), -sinf(asinf(-arg.y)*n));
}

__device__ inline
float2 complex_add(float2 const a, float2 b)
{
    return make_float2(a.x + b.x, a.y + b.y);
}

__device__ inline
float2 complex_minus(float2 a, float2 b)
{
    return make_float2(a.x - b.x, a.y - b.y);
}

__device__ inline
float2 complex_mult(float2 ab, float2 cd)
{
    return make_float2((ab.x * cd.x) - (ab.y * cd.y), (ab.x * cd.y) + (ab.y * cd.x));
}

__device__ inline
float2 complex_scale(float2 a, float s) 
{
	return make_float2(s * a.x, s * a.y);
}


__host__ __device__ inline
float fabsf2(const float2 &a) 
{ 
	return hypot( a.x , a.y );
	//return ; 
}


__device__ inline
float2 f2absf2(const float2 &a) 
{ 
	return make_float2(fabsf2(a),0.0); 
}


__global__
//__device__ inline 
 void complex_AddVec(	float2 * const A,
             			float2 const * const B,
                        float2 const * const C,
                        float const scale,
             		//	int const gi)
 						int const N)
{
	
    /* Calculate the global linear index, assuming a 1-d grid. */
    /*
    */
    int const gi = blockDim.x * blockIdx.x + threadIdx.x;
    if (gi < N) 
    {
	
        A[gi] = complex_add(B[gi] , complex_scale(C[gi],scale));
    }
}

 __global__ 
 void complex_ScaleVec(	float2 * const A,
             			float const b,
             			int const N,
             			//float2 const c0,
             			int const gi1,
             			int const gi2)
{
    /* Calculate the global linear index, assuming a 1-d grid. */
   	// 2D grid?
    int const gx = blockDim.x * blockIdx.x + threadIdx.x;
	int const gi = (blockDim.x * blockIdx.x) + threadIdx.x + (N * blockIdx.y);

    float pos, g, bc;
    //float a = 24.0;
	float a = 24.0;

    if (gx < N) 
    {
    	
    	if (gx > gi1 && gx < gi2) // Remove temporal aliasing     		
    	{
    		pos = (a*2.0/N) * (fabs( (float) gx - (N/2.0) ));
        	g = coshf(pos);
    		bc = (1 - (1/(g*g))) * b;
    		A[gi] = complex_scale(A[gi] , bc);
        	//A[gi] = complex_scale(A[gi] , 0.0*b);
    		//float bc = (cbrt(cbrt(cbrt(cbrt(cbrt( (float) pos)))))) * b;
    	}
    	else
    	{
    	
			A[gi] = complex_scale(A[gi] , b);
    	}
    	
    }
}
// Simple all in one NLFN kernel to start with.
/* NLfn = ((A0.^2).*exp(i(wt-(bdiffw0)*z)))
 * +(2*expi(-wt+(bdiffw0)*z).*abs(A0).^2)
 */
__global__
//__device__ inline 
void NLFN_kernel(float2 * const NL, 
				 float2 const * const Ax,
				 float2 const * const Ay,
				 double const dt,
				 float const w0, 
				 double bdwz,
				 //int gi) 
				 int Npoints)
{
	/* Calculate the global linear index, assuming a 1-d grid. */
	// 2D grid?
    int const gx = blockDim.x * blockIdx.x + threadIdx.x;
	int const gi = (blockDim.x * blockIdx.x) + threadIdx.x + (Npoints * blockIdx.y);
    if (gx < Npoints)
	/*
	for (int gi =  blockDim.x * blockIdx.x + threadIdx.x;
		gi < Npoints;
		gi += blockDim.x * gridDim.x)
	*/	
    {
    
		/* NL[gi] = ((A0[gi] * A0[gi]) * (expfi(wt[gi] - (bdiffw0 * z)))) +
    				(2 * expfi(-wt[gi] + (bdiffw0 * z)) * (abs(A0[gi]) * abs(A0[gi]) )); */

    	//float arg = ((w0 * t[gx]) - bdwz);
		float arg = ((w0 * (dt*gx)) - bdwz);

    	float2 AxAy		= complex_mult(Ax[gi],Ay[gi]);
    	float2 abAxAy 	= complex_mult(f2absf2(Ax[gi]),f2absf2(Ay[gi]));
    	float2 expon 	= complex_exp(arg);
    	float2 dbcnjex	= complex_add((conjugate(expon)),(conjugate(expon))); 
    	NL[gi] 			= complex_add(complex_mult(AxAy, expon), complex_mult(dbcnjex,abAxAy));
		// if (NL[gi].x == (FLT_MAX) || NL[gi].y == (FLT_MAX))
		// {
		// 	NL[gi] = make_float2(0.0,0.0);
		// }

    }

}

// Runge-Kutta KN Step Kernel
// K1=-h*1i*G*(NLfn-1/(2*dt)*(1i/w0)*(NLfn(ind1)-NLfn(ind2)))
__global__
//__device__ inline 
void KN_kernel (float2 const * const A0,
				float2 * const KN,
				float2 * const KN2,
				float2 * const NL,
				float const * const G,
				float const w0,
				float const h,
				int const k,
				float invdt2,
				//float2 const I,
				int Npoints,
				float2 * const AN,
				float b1,
                float b2)
{
	int in1, in2;

	/* Calculate the global linear index, assuming a 1-d grid. */
 	// 2D grid?
    int const gx = blockDim.x * blockIdx.x + threadIdx.x;
	int const gi = (blockDim.x * blockIdx.x) + threadIdx.x + (Npoints * blockIdx.y);
    
    int const d = 1;

	float2 iby2dtw0 = complex_mult(make_float2(invdt2/(d * w0), 0.0), I);
    float2 hiG = complex_mult(I, make_float2((-1.0 * h * G[k]),0.0));

    // if cell is not a boundary
    if (gx > d - 1 && gx < Npoints - d)
    {
		in1 = gi+d;
	    in2 = gi-d;
	}   
	else
	{
	    if (gx < d)
	    {
	    //in1 = gi+(2*d);
	    in1 = gi + d;	
	    in2 = Npoints - d + gi;
	    //in2 = gi;
	    }
	    if (gx > Npoints - (d + 1))
    	{
    	in1 = gi - (Npoints - d);
    	//in1 = gi;
    	//in2 = gi-(2*d);
    	in2 = gi - d;
	    }
	}

	float2 dNL = complex_minus(NL[in1], NL[in2]);
	KN[gi] = complex_mult(hiG, (complex_minus(NL[gi], complex_mult(iby2dtw0, dNL))));
	AN[gi] = complex_add(A0[gi],complex_add(complex_scale(KN[gi],b1),complex_scale(KN2[gi],b2)));
	// if (AN[gi].x == (FLT_MAX) || AN[gi].y == (FLT_MAX))
	// {
	// 	AN[gi] = make_float2(0.0,0.0);
	// }
	//AN[gi] = complex_add(A0[gi],complex_scale(complex_add(KN[gi],KN2[gi]),fracstep));
}


__global__ void Dispersion_kernel (	float2 * const E,
									float2 * const EFT,
									float const * const alpha_w,
									float const * const beta_op,
									int Npoints,
									int stepmod,
									int chunk,
									int Nchunks,
									float * const stepmods,
									int const k)									
{
	/* Calculate the global linear index, assuming a 1-d grid. */
    //int const gi = blockDim.x * blockIdx.x + threadIdx.x;
    //int sc;
    float2 bophmod = make_float2(1.0,0.0);
	//float2 bophmod;

	// 2D grid?
    int const gx = blockDim.x * blockIdx.x + threadIdx.x;
	int const gy = blockIdx.y;
	int const gi = (blockDim.x * blockIdx.x) + threadIdx.x + (Npoints * blockIdx.y);

	if (gx == k)
	{
		stepmods[k] = stepmod;
	}

    if (gx < Npoints)
    {
    	
    	// for (sc=0;sc<stepmod;sc++)
		// {
		// 	bophmod = complex_mult(bophmod,beta_op[gi]);
		// }
	
    	//bophmod = complex_pow(beta_op[gi],stepmod);

    	bophmod = complex_exp(-beta_op[gx] * stepmod);
		bophmod = complex_scale(bophmod, powf(alpha_w[gx],stepmod));
    	EFT[gi] = complex_mult(EFT[gi],bophmod);
/*
    	if (gi > 200 && gi < 3000)
    	{
        	EFT[gi] = complex_scale(EFT[gi] , 0.001);
    	}
*/
    	//E[chunk + (gi * Nchunks)] = EFT[gi];
        E[(gy * (Nchunks-1) * Npoints) + (chunk * Npoints) + gi] = EFT[gi];
    }

}
									

__global__ void Error_kernel (	float2 * const K1,
								float2 * const K2,
								float2 * const E,
								int k,
								int stepmod,
								int Npoints,
								float const maxpe,
								float const minpe,
								bool * max_min_flag)
{
	/* Calculate the global linear index, assuming a 1-d grid. */
    //int const gi = blockDim.x * blockIdx.x + threadIdx.x;
    float errfv, pcterrfv;

  	// 2D grid?
    int const gx = blockDim.x * blockIdx.x + threadIdx.x;
	int const gi = (blockDim.x * blockIdx.x) + threadIdx.x + (Npoints * blockIdx.y);
    if (gx < Npoints)
    {
    	if(fabsf2(E[gi]) > 0)
    	{
	    	errfv 	 = 0.5 * (fabsf2(complex_minus(K2[gi],K1[gi])));
			//errfv 	 =  0.01 * fabsf2(make_float2(1.0,0.0));
			//errfv 	 =  0.01;
	    	pcterrfv = 4 * 100.0 * (errfv / (fabsf2(E[gi])+1));
			//pcterrfv = 0.01;

	    	if (pcterrfv >= maxpe)
	    	{
	    		max_min_flag[0] = true;
	    	} 
	    	else
	    	{
		    	if (pcterrfv >= minpe)
		    	{
		    		max_min_flag[1] = true;
		    	} 
    		}
	 	}
    }

}


__global__ void copyKernel (float2 * Eft,
							float2 * ApFT,
							int Npoints,
							int chunk,
							int Nchunks)
{
	/* Calculate the global linear index, assuming a 1-d grid. */
    //int const gi = blockDim.x * blockIdx.x + threadIdx.x;

    //if (gi < Npoints)
	for (int gi =  blockDim.x * blockIdx.x + threadIdx.x;
		gi < Npoints;
		gi += blockDim.x * gridDim.x)	
    {
    	ApFT[chunk + (gi * Nchunks)] = Eft[gi];
    }

} 

// Host function called by MEX gateway.
 void OPO_TEST_CUDA(float2 * const d_E,
					float2 * const d_A0,
					float2 * const d_A1,
					float2 * const d_NL,
					float2 * d_K1,
					float2 * const d_K2,
                    float2 * const d_K3,
                    float2 * d_K4,    
				   	float const * const d_G,
					float const * const d_alpha_w,
				   	float const w0, 
				   	double const bdiffw0,
				   	float const h,
				   	int const Nsteps,
				   	int Npoints,
					int Nbatches,
				   	double const dt,
					float const * const d_beta_op,
					float const max_pct_err,
					float const min_pct_err,
					int const sel,
					float2 * const d_Ap,
					float * const d_stepmods)
{	
	// GPU variables
	float2* d_temp;
	float stepmod, stepend, stepprev, k, invdt2, invSize, hmod;
	double bdwz;
	int  chunk, Nchunks, gi1, gi2;
	bool* max_min_err_exc;

	// Precompute scalars:
	invdt2 	= 1.0 / (2.0 * dt);
	invSize = 1.0 / Npoints;
	k = 0.0;
	stepprev = 4.0;
	Nchunks = Nsteps / sel;
	// Constants defining cut off for temporal aliasing reduction
	gi1 = 6.0*(Npoints/16);
	gi2 = 10.0*(Npoints/16);

	cudaMallocHost(&max_min_err_exc, 2*sizeof(bool));

	// CUFFT plan simple API
  	cufftHandle plan;
  	cufftPlan1d(&plan, Npoints, CUFFT_C2C, Nbatches);
	//cufftPlanMany(&plan, 1, Npoints, CUFFT_C2C, 1);

	//dim3 const blocksPerGrid(((Npoints) + threadsPerBlock - 1) / threadsPerBlock, Nbatches);
    dim3 const blocksPerGrid(Npoints / threadsPerBlock, Nbatches);
	//err_reset_kernel<<<1, 1 >>>();

    // Compute scalars
	const float butcher_as[] = {1.0/2.0, 3.0/4.0, 4.0/9.0, 1.0/3.0, 2.0/9.0};
	//const float butcher_as[] = {1.0/2.0, 1.0/4.0, 0.0, 1.0/4.0, 1.0/4.0};
	const float butcher_bs[] = {1.0/8.0, 1.0/3.0, 1.0/4.0, 7.0/24.0};
	//const float butcher_bs[] = {1.0/4.0, 1.0/4.0, 1.0/4.0, 1.0/4.0};
    stepmod = stepprev;
	hmod = stepmod * h;
    bdwz = bdiffw0 * k * h;

    NLFN_kernel<<<blocksPerGrid, threadsPerBlock>>>(d_NL,d_A0,d_A0,dt,w0,bdwz,Npoints);
	//KN_kernel<<<blocksPerGrid, threadsPerBlock>>>(d_A0,d_K1,d_K1,d_NL,d_G,w0,hmod,k,invdt2,Npoints,d_A1,0.5,0.0);
	KN_kernel<<<blocksPerGrid, threadsPerBlock>>>(d_A0,d_K1,d_K1,d_NL,d_G,w0,hmod,k,invdt2,Npoints,d_A1,butcher_as[0],0.0);

	for (chunk = 0; chunk < Nchunks; chunk++)
	{
		stepmod = stepprev;

		hmod = stepmod * h;

		for (k = k; k < (chunk + 1) * sel; k = k+stepmod)
		{	
            max_min_err_exc[0] = max_min_err_exc[1] = false;

			bdwz = bdiffw0 * (k + (butcher_as[0] * stepmod)) * h;// half-step co-ordinate coeff

			NLFN_kernel<<<blocksPerGrid, threadsPerBlock>>>(d_NL,d_A1,d_A1,dt,w0,bdwz,Npoints);
			KN_kernel<<<blocksPerGrid, threadsPerBlock>>>(d_A0,d_K2,d_K2,d_NL,d_G,w0,hmod,k,invdt2,Npoints,d_A1,butcher_as[1],0.0);	
			
            bdwz = bdiffw0 * (k + (butcher_as[1] * stepmod)) * h;// three-quarter-step co-ordinate coeff

            NLFN_kernel<<<blocksPerGrid, threadsPerBlock>>>(d_NL,d_A1,d_A1,dt,w0,bdwz,Npoints);
            KN_kernel<<<blocksPerGrid, threadsPerBlock>>>(d_A0,d_K3,d_K2,d_NL,d_G,w0,hmod,k,invdt2,Npoints,d_E,butcher_as[2],butcher_as[3]);	
             
			complex_AddVec<<<blocksPerGrid, threadsPerBlock>>>(d_E,d_E,d_K1,butcher_as[4],Npoints);

            bdwz = bdiffw0 * (k + stepmod) * h;// full-step co-ordinate coeff

            NLFN_kernel<<<blocksPerGrid, threadsPerBlock>>>(d_NL,d_E,d_E,dt,w0,bdwz,Npoints);
			KN_kernel<<<blocksPerGrid, threadsPerBlock>>>(d_A0,d_K4,d_K3,d_NL,d_G,w0,hmod,k,invdt2,Npoints,d_A1,butcher_bs[0],butcher_bs[1]);
             
			complex_AddVec<<<blocksPerGrid, threadsPerBlock>>>(d_A1,d_A1,d_K2,butcher_bs[2],Npoints);
			complex_AddVec<<<blocksPerGrid, threadsPerBlock>>>(d_A1,d_A1,d_K1,butcher_bs[3],Npoints);

            cudaDeviceSynchronize();
			//Error_kernel<<<blocksPerGrid, threadsPerBlock>>>(d_K1, d_K2, d_E, k, stepmod, Npoints, 
			Error_kernel<<<blocksPerGrid, threadsPerBlock>>>(d_A1, d_E, d_E, k, stepmod, Npoints, 
																max_pct_err, min_pct_err, max_min_err_exc);															
			//max_min_err_exc[1] = true;															
			cudaDeviceSynchronize();													
			
			if (max_min_err_exc[0] && stepmod > 1 && (k - stepmod) > 0)
			{
				k = k - stepmod;
				stepmod = stepmod - 1;
				hmod 	= stepmod * h;
				stepprev = stepmod;
				
				bdwz = bdiffw0 * (k + stepmod) * h;

    			NLFN_kernel<<<blocksPerGrid, threadsPerBlock>>>(d_NL,d_A0,d_A0,dt,w0,bdwz,Npoints);
				KN_kernel<<<blocksPerGrid, threadsPerBlock>>>(d_A0,d_K1,d_K1,d_NL,d_G,w0,hmod,k,invdt2,Npoints,d_A1,butcher_as[0],0.0);
				//complex_AddVec<<<blocksPerGrid, threadsPerBlock>>>(d_A1,d_A0,d_K1,butcher_as[0],Npoints);
			}
			else 
			{
				complex_ScaleVec<<<blocksPerGrid, threadsPerBlock>>>(d_E, 1, Npoints, gi1, gi2);
	
				// Transform signal
				cufftExecC2C(plan, (cufftComplex *) d_E,
								(cufftComplex *) d_E,
									CUFFT_FORWARD);						

				//Dispersion_kernel<<<Npoints/threadsPerBlock, threadsPerBlock>>>(d_Ap+(i * Nchunks), d_E+i, d_beta_op, Npoints, stepmod, chunk, Nchunks);	
				Dispersion_kernel<<<blocksPerGrid, threadsPerBlock>>>(d_Ap, d_E, d_alpha_w, d_beta_op, Npoints, stepmod, chunk, Nchunks, d_stepmods, k);
				//cudaDeviceSynchronize();

				// Transform signal back
				cufftExecC2C(plan, (cufftComplex *) d_E,
								(cufftComplex *) d_A0,
									CUFFT_INVERSE);

				complex_ScaleVec<<<blocksPerGrid, threadsPerBlock>>>(d_A0, invSize, Npoints, 0, 0); 
/*
                */

                d_temp = d_K1;
                d_K1 = d_K4;
                d_K4 = d_temp;
				cudaDeviceSynchronize();
                complex_AddVec<<<blocksPerGrid, threadsPerBlock>>>(d_A1,d_A0,d_K1,butcher_as[0],Npoints);
 				
				
				if (k + stepmod >= (chunk + 1) * sel && k < ((chunk + 1) * sel) - 1)
				{
					
					stepend = (((chunk + 1) * sel) - k) - 1; 
					hmod 	= stepend * h;
					stepprev = stepmod;
					stepmod = stepend;
				}
				else
				{
					if (!max_min_err_exc[1] && !max_min_err_exc[0] && (k + stepmod + stepmod) < ((chunk + 1) * sel) - 1)
					{
						stepmod = stepmod * 2;
						//stepmod = stepmod + 1;
						hmod 	= stepmod * h;
					}
					
					stepprev = stepmod;

				}
			

			} // dispersion
			
		
		} // chunk
		cudaDeviceSynchronize();

	} // crystal		
	cudaDeviceSynchronize();

  	// Destroy CUFFT context
	cufftDestroy(plan);

	// Release resources
	cudaFreeHost(max_min_err_exc);
	//cudaFree(max_min_err_exc);

}