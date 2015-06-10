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
__kernel void advance_ptcls_axial_damping(
// Positions of particles
    __global REAL* xGlob,
    __global REAL* yGlob,
    __global REAL* zGlob,
// Velocities of particles
    __global REAL* vxGlob,
    __global REAL* vyGlob,
    __global REAL* vzGlob,
// Charge of particles
    const __global REAL* charge,
// Mass of particles
    const __global REAL* mass,
// Magnetic field strength along z in T
    REAL eToTheMinusKappaDt,
// Number of particles
    int nPtcls)
{
  __private int n = get_global_id(0);
  if (n < nPtcls) {
    vzGlob[n] *= eToTheMinusKappaDt;
  }
}

