
/* OPO_ODE23_HEADER.h */
#include <cufft.h>
#include <matrix.h>
#ifndef OPO_ODE23_HEADER_HPP
#define OPO_ODE23_HEADER_HPP


	void OPO_TEST_CUDA(float2* const d_E,
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
					float * const d_stepmods);

 #endif 