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

REAL lineShape(REAL gamma, REAL delta, REAL S);
REAL peakScatteringRate(REAL3 k, REAL3 v,
    REAL gamma, REAL delta, REAL S);
REAL gaussianBeamProfile(REAL3 k, REAL3 r0, REAL sigma, REAL3 r);

__kernel void compute_mean_scattered_photons_homogeneous_beam(
// positions of particles
    __global REAL* xGlob,
    __global REAL* yGlob,
    __global REAL* zGlob,
    __global REAL* vxGlob,
    __global REAL* vyGlob,
    __global REAL* vzGlob,
// Wave vector of cooling laser
    REAL kx, REAL ky, REAL kz,
// Transition linewidth
    REAL gamma,
// Detuning at zero velocity
    REAL delta0,
// Saturation parameter
    REAL S,
// Time step size
    REAL dt,
// Number of particles
    int nPtcls,
// Number of scattered photons for each particle 
    __global REAL* nbar
    )
{
  int i = get_global_id(0);
  if (i < nPtcls) {
    REAL3 v = (REAL3)(vxGlob[i], vyGlob[i], vzGlob[i]);
    REAL3 k = (REAL3)(kx, ky, kz);

    nbar[i] = dt * peakScatteringRate(k, v, gamma, delta0, S);
  }
}

__kernel void compute_mean_scattered_photons_gaussian_beam(
// positions of particles
    __global REAL* xGlob,
    __global REAL* yGlob,
    __global REAL* zGlob,
    __global REAL* vxGlob,
    __global REAL* vyGlob,
    __global REAL* vzGlob,
// Wave vector of cooling laser
    REAL kx, REAL ky, REAL kz,
// A position on the fluorescence peak
    REAL x0, REAL y0, REAL z0,
// 1/e radius of Gaussian cooling laser
    REAL sigma,
    REAL gamma,
    REAL delta0,
    REAL S,
    REAL dt,
// Number of particles
    int nPtcls,
// Number of scattered photons for each particle 
    __global REAL* nbar
    )
{
  int i = get_global_id(0);
  if (i < nPtcls) {
    REAL3 v = (REAL3)(vxGlob[i], vyGlob[i], vzGlob[i]);
    REAL3 r = (REAL3)(xGlob[i], yGlob[i], zGlob[i]);
    REAL3 k = (REAL3)(kx, ky, kz);
    REAL3 r0 = (REAL3)(x0, y0, z0);

    REAL nPeak = peakScatteringRate(k, v, gamma, delta0, S);
    REAL geomFact = gaussianBeamProfile(k, r0, sigma, r);
    nbar[i] = dt * geomFact * nPeak;
  }
}

__kernel void countEmissions(
    __global REAL* meanN,
    __global REAL* auxRandNumbers,
    int nMax,
    __global int* numEmissions,
    int numPtcls
    )
{
  int i = get_global_id(0);
  if (i < numPtcls) {
    REAL z = meanN[i] / nMax;
    auxRandNumbers += i;
    int actualN = 0;
    for (int j = 0; j < nMax; ++j) {
      if (*auxRandNumbers < z) {
        actualN += 1;
      }
      auxRandNumbers += numPtcls;
    }
    numEmissions[i] = actualN;
  }
}

__kernel void computeKicks(
    __global REAL* mass,
    __global int* actualNs,
    int nMax,
    __global REAL* dirs,
    REAL kx0, REAL ky0, REAL kz0,
    REAL recoilMomentum,
    REAL dt,
    __global REAL* ax, __global REAL* ay, __global REAL* az,
    int nPtcls
    )
{
  int i = get_global_id(0);
  if (i < nPtcls) {
    REAL m = mass[i];
    REAL recoilVelocity = recoilMomentum / m;
    int myNumPhotons = actualNs[i];
    int totNumPhotons = nMax * nPtcls;
    REAL3 vRecoil = myNumPhotons * recoilVelocity *
      normalize((REAL3)(kx0, ky0, kz0));

    dirs += i * nMax;
    for (int j = 0; j < myNumPhotons; ++j) {
      REAL3 dir;
      dir.x = dirs[0 * totNumPhotons + j];
      dir.y = dirs[1 * totNumPhotons + j];
      dir.z = dirs[2 * totNumPhotons + j];
      vRecoil += recoilVelocity * normalize(dir);
    }
    ax[i] += vRecoil.x / dt;
    ay[i] += vRecoil.y / dt;
    az[i] += vRecoil.z / dt;
  }
}


REAL lineShape(REAL gamma, REAL delta, REAL S)
{
  REAL halfGammaSquared = pown(0.5 * gamma, 2);
  return halfGammaSquared / 
    (halfGammaSquared * (1.0 + 2.0 * S) + pown(delta, 2));
}

REAL peakScatteringRate(REAL3 k, REAL3 v,
    REAL gamma, REAL delta, REAL S)
{
  REAL dopplerShift = -dot(k, v);
  REAL totalDetuning = delta + dopplerShift;
  return S * gamma / (2.0 * M_PI_F) * lineShape(gamma, totalDetuning, S);
}

REAL gaussianBeamProfile(REAL3 k, REAL3 r0, REAL sigma, REAL3 r)
{
  REAL3 dr = r - r0;
  REAL3 kHat = k / length(k);
  REAL rho = length(dr - dot(dr, kHat) * kHat);
  return exp(-pown(length(rho) / sigma, 2));
}

