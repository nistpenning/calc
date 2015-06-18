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


REAL3 two_body_interaction(REAL4 ptcl1, REAL4 ptcl2, REAL impactFact);

__kernel void compute_coulomb_acceleration(
    // positions of particles
    const __global REAL* xGlob,
    const __global REAL* yGlob,
    const __global REAL* zGlob,
    // velocities of particles
    const __global REAL* vxGlob,
    const __global REAL* vyGlob,
    const __global REAL* vzGlob,
    // charge of particles
    const __global REAL* charge,
    // mass of particles
    const __global REAL* mass,
    REAL k,
    REAL impactFact,
    int nPtcls,
    // accelerations (output)
    __global REAL* axGlob,
    __global REAL* ayGlob,
    __global REAL* azGlob
    ) {

  __local REAL4 xyzqCache[BLOCK_SIZE];

  __private int n = get_global_id(0) + get_global_id(1) * get_global_size(0);
  __private REAL4 ptcl1[PTCL_UNROLL_FACTOR];
  __private REAL m1[PTCL_UNROLL_FACTOR];
  __private REAL3 acceleration[PTCL_UNROLL_FACTOR];

  for (int i = 0; i < PTCL_UNROLL_FACTOR; ++i) {
    int ptclId = n + i * get_global_size(0) * get_global_size(1);
    ptcl1[i].x = (ptclId < nPtcls) ? xGlob[ptclId] : 0;
    ptcl1[i].y = (ptclId < nPtcls) ? yGlob[ptclId] : 0;
    ptcl1[i].z = (ptclId < nPtcls) ? zGlob[ptclId] : 0;
    ptcl1[i].w = (ptclId < nPtcls) ? charge[ptclId] : 0;
    m1[i] = (ptclId < nPtcls) ? mass[ptclId] : 0;
    acceleration[i] = (REAL3)0;
  }

  int numTiles = nPtcls / BLOCK_SIZE;
  for (int tile = get_group_id(0); tile < get_group_id(0) + numTiles; ++tile) {
    int globInd;
    globInd = get_local_id(0) + (tile % numTiles) * BLOCK_SIZE;
    xyzqCache[get_local_id(0)].x = xGlob[globInd];
    xyzqCache[get_local_id(0)].y = yGlob[globInd];
    xyzqCache[get_local_id(0)].z = zGlob[globInd];
    xyzqCache[get_local_id(0)].w = charge[globInd];
    barrier(CLK_LOCAL_MEM_FENCE);
#pragma unroll 8
    for (int jj = 0; jj < BLOCK_SIZE; ++jj) {
      int jjStaggered = (jj + get_local_id(0)) & (BLOCK_SIZE - 1);
      REAL4 ptcl2 = xyzqCache[jjStaggered];
      for (int i = 0; i < PTCL_UNROLL_FACTOR; ++i) {
        acceleration[i] += two_body_interaction(ptcl1[i], ptcl2, impactFact);
      }
    }
    barrier(CLK_LOCAL_MEM_FENCE);
  }

  for (int j = numTiles * BLOCK_SIZE; j < nPtcls; ++j) {
    REAL4 ptcl2 = (REAL4)(xGlob[j], yGlob[j], zGlob[j], charge[j]);
    for (int i = 0; i < PTCL_UNROLL_FACTOR; ++i) {
      acceleration[i] += two_body_interaction(ptcl1[i], ptcl2, impactFact);
    }
  }
  for (int i = 0; i < PTCL_UNROLL_FACTOR; ++i) {
    int ptclId = n + i * get_global_size(0) * get_global_size(1);
    if (ptclId < nPtcls) {
      REAL myFact = k * ptcl1[i].w / m1[i];
      axGlob[ptclId] = myFact * acceleration[i].x;
      ayGlob[ptclId] = myFact * acceleration[i].y;
      azGlob[ptclId] = myFact * acceleration[i].z;
    }
  }
}

REAL3 two_body_interaction(REAL4 ptcl1, REAL4 ptcl2, REAL impactFact) {
  REAL3 dr = ptcl2.xyz - ptcl1.xyz;
  REAL distSqr = fma(dr.x, dr.x, fma(dr.y, dr.y, fma(dr.z, dr.z, EPS2 + impactFact)));
  REAL invDist;
#ifdef USE_DOUBLE
  invDist = rsqrt(distSqr);
#else
  invDist = native_rsqrt(distSqr);
#endif
  REAL invDistCube = invDist * invDist * invDist;
  REAL s = ptcl2.w * invDistCube;
  REAL3 acc = -s * dr;
  return acc;
}

