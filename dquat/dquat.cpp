#include <cmath>
#include <iostream>
#include <random>

constexpr size_t invdim = 4;
constexpr size_t emapdim = 3;
constexpr size_t nwvec = 3;
constexpr size_t qdim = 4;

constexpr double idp_tiny_sqrt = 1.0e-90;
constexpr double zero = 0.0;
constexpr double one = 1.0;
constexpr double onehalf = 0.5;
constexpr double oneqrtr = 0.25;

#define ECMECH_NM_INDX(p, q, qDim) ((p) * (qDim) + (q))

// Simple scalar product scaling function
template <int n>
inline void vecsVxa(double *const v, double x, const double *const a) {
  for (int i = 0; i < n; ++i) {
    v[i] = x * a[i];
  }
}

// The L2 norm of a vector
template <int n> inline double vecNorm(const double *const v) {
  double retval = 0.0;
  for (int i = 0; i < n; ++i) {
    retval += v[i] * v[i];
  }

  retval = sqrt(retval);
  return retval;
}

inline void inv_to_quat(double *const quat, const double *const inv) {
  double a = inv[0] * 0.5;
  quat[0] = cos(a);
  a = sin(a);
  vecsVxa<nwvec>(&(quat[1]), a, &(inv[1]));
}

inline void emap_to_quat(double *const quat, const double *const emap) {
  double inv[invdim] = {0.0, 1.0, 0.0, 0.0};
  inv[0] = vecNorm<emapdim>(emap);
  if (inv[0] > idp_tiny_sqrt) {
    double invInv = 1.0 / inv[0];
    vecsVxa<emapdim>(&(inv[1]), invInv, emap);
  } // else, emap is effectively zero, so axis does not matter
  inv_to_quat(quat, inv);
}

/**
    ! derivative of quaternion parameters with respect to exponential map
    ! parameters
    ! This is the derivative of emap_to_quat function with a transpose applied
*/
inline void dquat_demap_T(double *const dqdeT_raw, // (emapdim, qdim)
                          const double *const emap // (emapdim)
) {
  const double theta_sm_a = 1e-9;
  const double oo48 = 1.0 / 48.0;

  double theta = vecNorm<emapdim>(emap);

  double theta_inv, sthhbyth, halfsthh, na, nb, nc;
  // When dealing with exponential maps that have very small rotations
  // along an axis the derivative term can become difficult to
  // calculate with the analytic solution
  //
  if (fabs(theta) < theta_sm_a) {
    sthhbyth =
        onehalf -
        theta * theta *
            oo48; // truncated Taylor series; probably safe to just use onehalf and be done with it
    halfsthh = theta * oneqrtr; // truncated Taylor series
    if (fabs(theta) < idp_tiny_sqrt) {
      // n is arbitrary, as theta is effectively zero
      na = one;
      nb = zero;
      nc = zero;
    } else {
      theta_inv = one / theta;
      na = emap[0] * theta_inv;
      nb = emap[1] * theta_inv;
      nc = emap[2] * theta_inv;
    }
  } else {
    halfsthh = sin(theta * onehalf);
    sthhbyth = halfsthh / theta;
    halfsthh = halfsthh * onehalf;
    theta_inv = one / theta;
    na = emap[0] * theta_inv;
    nb = emap[1] * theta_inv;
    nc = emap[2] * theta_inv;
  }
  //
  double halfcthh = cos(theta * onehalf) * onehalf;
  //
  // now have: halfsthh, sthhbyth, halfcthh, theta, na, nb, nc
  // If we made use of something like RAJA or Kokkos views or mdspans then we could
  // make some of the indexing here nicer but don't worry about that now
  // RAJA::View<double, RAJA::Layout<2> > dqdeT(dqdeT_raw, emapdim, qdim);

  dqdeT_raw[ECMECH_NM_INDX(0, 0, qdim)] = -halfsthh * na;
  dqdeT_raw[ECMECH_NM_INDX(1, 0, qdim)] = -halfsthh * nb;
  dqdeT_raw[ECMECH_NM_INDX(2, 0, qdim)] = -halfsthh * nc;

  double temp = na * na;
  dqdeT_raw[ECMECH_NM_INDX(0, 1, qdim)] =
      halfcthh * temp + sthhbyth * (one - temp);
  //
  temp = nb * nb;
  dqdeT_raw[ECMECH_NM_INDX(1, 2, qdim)] =
      halfcthh * temp + sthhbyth * (one - temp);
  //
  temp = nc * nc;
  dqdeT_raw[ECMECH_NM_INDX(2, 3, qdim)] =
      halfcthh * temp + sthhbyth * (one - temp);

  temp = halfcthh - sthhbyth;
  //
  double tempb = temp * na * nb;
  dqdeT_raw[ECMECH_NM_INDX(1, 1, qdim)] = tempb;
  dqdeT_raw[ECMECH_NM_INDX(0, 2, qdim)] = tempb;
  //
  tempb = temp * na * nc;
  dqdeT_raw[ECMECH_NM_INDX(2, 1, qdim)] = tempb;
  dqdeT_raw[ECMECH_NM_INDX(0, 3, qdim)] = tempb;
  //
  tempb = temp * nb * nc;
  dqdeT_raw[ECMECH_NM_INDX(2, 2, qdim)] = tempb;
  dqdeT_raw[ECMECH_NM_INDX(1, 3, qdim)] = tempb;
}

int main() {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(-M_PI, M_PI);

  const int num_tests = 100000;

  for (int test = 0; test < num_tests; ++test) {
    double emap[emapdim] = {};

    // Populate 'emap' with random values
    for (size_t i = 0; i < emapdim; ++i) {
      emap[i] = dis(gen);
    }

    double quat[qdim] = {};
    emap_to_quat(quat, emap);

    double dquat_dexpmap_t[emapdim * qdim] = {};
    dquat_demap_T(dquat_dexpmap_t, emap);

    // Output the results for this test
    std::cout << "Test " << test + 1 << ":" << std::endl;

    std::cout << "exponential map value:" << std::endl;
    for (size_t ie = 0; ie < emapdim; ie++) {
      std::cout << emap[ie] << " ";
    }
    std::cout << std::endl;

    std::cout << "unit quaternion value:" << std::endl;
    for (size_t iq = 0; iq < qdim; iq++) {
      std::cout << quat[iq] << " ";
    }
    std::cout << std::endl;

    std::cout << "dquat_dexpmap ^ T value:" << std::endl;
    for (size_t ie = 0; ie < emapdim; ie++) {
      for (size_t iq = 0; iq < qdim; iq++) {
        std::cout << dquat_dexpmap_t[ECMECH_NM_INDX(ie, iq, qdim)] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << "----------------------------------------" << std::endl;
  }

  return 0;
}
