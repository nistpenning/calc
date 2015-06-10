#ifdef USE_DOUBLE
#ifdef cl_khr_fp64
#pragma OPENCL EXTENSION cl_khr_fp64: enable
#define REAL double
#define REAL3 double3
#define REAL4 double4
#define EPS2 1.0e-30
#endif
#else
#define REAL float
#define REAL3 float3
#define REAL4 float4
#define EPS2 1.0e-18f
#endif

// Advances particles moving in a uniform magnetic field along z
__kernel void advance_ptcls_cooling_along_x(
// Positions of particles
    __global REAL* xGlob,
    __global REAL* yGlob,
    __global REAL* zGlob,
// Velocities of particles
    __global REAL* vxGlob,
    __global REAL* vyGlob,
    __global REAL* vzGlob,
// 3 * nPtcls normally distributed random numbers
    __global REAL* randNumbers,
// Charge of particles
    const __global REAL* charge,
// Mass of particles
    const __global REAL* mass,
// Cooling rate at the intensity peak
    REAL peakCoolingRate,
// Doppler heating diffusion constant at the peak
    REAL peakDiffusionConstant,
// Offset of cooling laser from y = 0
    REAL offset,
// 1/e width of cooling laser beam (Gaussian intensity profile is assumed
    REAL width,
// Time step
    REAL dt,
// Number of particles
    int nPtcls)
{
  __private int n = get_global_id(0);
  if (n < nPtcls) {
    __private REAL vx = vxGlob[n];
    __private REAL vy = vyGlob[n];
    __private REAL vz = vzGlob[n];
    __private REAL y = yGlob[n];

    REAL profileFactor = exp(-pown((y - offset) / width, 2));

    vx *= exp(-peakCoolingRate * profileFactor * dt);

    REAL diffusionFactor =
      peakDiffusionConstant * sqrt(profileFactor * dt)/ 3.0; 
    vx += diffusionFactor * randNumbers[n + 0 * nPtcls];
    vy += diffusionFactor * randNumbers[n + 1 * nPtcls];
    vz += diffusionFactor * randNumbers[n + 2 * nPtcls];

    vxGlob[n] = vx;
    vyGlob[n] = vy;
    vzGlob[n] = vz;
  }
}


