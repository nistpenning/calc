#ifdef USE_DOUBLE
#ifdef cl_khr_fp64
#pragma OPENCL EXTENSION cl_khr_fp64: enable
#define REAL double
#define REAL2 double2
#define REAL3 double3
#define REAL4 double4
#endif
#else
#define REAL float
#define REAL2 float2
#define REAL3 float3
#define REAL4 float4
#endif

__kernel void advance_ptcls_angular_damping(
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
// Damping coefficient
    REAL eToTheMinusKappaDT,
// Desired angular rotation frequency
    REAL omega,
    REAL minRadius,
// Number of particles
    int nPtcls)
{
  __private int n = get_global_id(0);
  if (n < nPtcls) {
    __private REAL2 v = (REAL2)(vxGlob[n], vyGlob[n]);
    __private REAL2 r = (REAL2)(xGlob[n], yGlob[n]);
    __private REAL radius = length(r);

    if (radius > minRadius) {
      __private REAL2 vHat = normalize((REAL2)(-r.y, r.x));
      __private REAL2 vPar = vHat * dot(vHat, v);
      __private REAL2 targetVelocity = omega * radius * vHat;
      __private REAL2 vPerp = v - vPar;

// update the velocities
      v = vPerp + targetVelocity + (vPar - targetVelocity) *
        eToTheMinusKappaDT;
      vxGlob[n] = v.x;
      vyGlob[n] = v.y;
    }
  }
}


