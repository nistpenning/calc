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

__kernel void compute_trap_acceleration(
// positions of particles
    const __global REAL* xGlob,
    const __global REAL* yGlob,
    const __global REAL* zGlob,
// charge of particles
    const __global REAL* charge,
// mass of particles
    const __global REAL* mass,
// Trap parameters
    REAL kx,
    REAL ky,
    REAL kz,
// cos and sin of theta, the angle of the rotating wall potential
// relative to the x-y coordinate system
    REAL cosTheta, REAL sinTheta,
// Number of particles
    int nPtcls,
// accelerations (output)
    __global REAL* axGlob,
    __global REAL* ayGlob,
    __global REAL* azGlob)
{
  __private int n = get_global_id(0);
  if (n < nPtcls) {
    __private REAL x = xGlob[n];
    __private REAL y = yGlob[n];
    __private REAL z = zGlob[n];

    __private REAL chargeOverMass = charge[n] / mass[n];

    __private REAL ax =
      (-kx * pown(cosTheta, 2) - ky * pown(sinTheta, 2)) * x +
      cosTheta * sinTheta * (ky - kx) * y;
    __private REAL ay =
      cosTheta * sinTheta * (-kx + ky) * x +
      (-kx * pown(sinTheta, 2) - ky * pown(cosTheta, 2)) * y;
    __private REAL az = - kz * z;

    ax *= chargeOverMass;
    ay *= chargeOverMass;
    az *= chargeOverMass;

    axGlob[n] += ax;
    ayGlob[n] += ay;
    azGlob[n] += az;
  }
}

