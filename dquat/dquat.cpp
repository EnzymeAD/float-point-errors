#include <chrono>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>

#ifdef LOGGED
#include "fp-logger.hpp"

void thisIsNeverCalledAndJustForTheLinker() {
  enzymeLogError("", 0.0);
  enzymeLogGrad("", 0.0);
  enzymeLogValue("", 0.0, 2, nullptr);
}

int enzyme_dup;
int enzyme_dupnoneed;
int enzyme_out;
int enzyme_const;

template <typename return_type, typename... T>
return_type __enzyme_autodiff(void *, T...);
#endif

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

__attribute__((noinline)) void emap_to_quat(double *const quat,
                                            const double *const emap) {
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
__attribute__((noinline)) void
dquat_demap_T(double *const dqdeT_raw, // (emapdim, qdim)
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

int main(int argc, char *argv[]) {
#ifdef LOGGED
  initializeLogger();
#endif
  const int num_tests = 1e6;
  std::string output_path = "";
  bool save_output = false;

  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i], "--output-path") == 0) {
      if (i + 1 < argc) {
        output_path = argv[++i];
        save_output = true;
      } else {
        std::cerr << "Error: --output-path option requires a value"
                  << std::endl;
        return 1;
      }
    } else {
      std::cerr << "Unknown option: " << argv[i] << std::endl;
      std::cerr << "Usage: " << argv[0] << " [--output-path <path>]"
                << std::endl;
      return 1;
    }
  }

  std::ostream *out_stream = &std::cout;
  std::ofstream outfile;
  if (save_output) {
    outfile.open(output_path);
    if (!outfile) {
      std::cerr << "Error: Could not open output file: " << output_path
                << std::endl;
      return 1;
    }
    outfile << std::setprecision(std::numeric_limits<double>::digits10 + 1)
            << std::scientific;
    out_stream = &outfile;
  }

  std::mt19937 gen;
  gen.seed(42);

  std::uniform_real_distribution<> dis_case1(-1e-9,
                                             1e-9); // 1e-90 <= theta < 1e-9
  std::uniform_real_distribution<> dis_case2(-1e-1, 1e-1);   // theta >= 1e-9
  std::uniform_real_distribution<> dis_case3(-1e-90, 1e-90); // theta < 1e-90

  const int num_cases = 3;

  auto start_time = std::chrono::high_resolution_clock::now();

  for (int test = 0; test < num_tests; ++test) {
    double emap[emapdim] = {};

    // Determine the current case by cycling through the cases
    int current_case = test % num_cases;

    switch (current_case) {
    case 0:
      // Case 0: theta < 1e-90 (effectively zero)
      for (size_t i = 0; i < emapdim; i++) {
        emap[i] = dis_case3(gen);
      }
      break;
    case 1:
      // Case 1: 1e-90 <= theta < 1e-9 (small rotations)
      for (size_t i = 0; i < emapdim; i++) {
        emap[i] = dis_case1(gen);
      }
      break;
    case 2:
      // Case 2: theta >= 1e-9 (larger rotations)
      for (size_t i = 0; i < emapdim; i++) {
        emap[i] = dis_case2(gen);
      }
      break;
    default:
      // This should never happen
      std::cerr << "Invalid case encountered!" << std::endl;
      return 1;
    }

    // double quat[qdim] = {};
    // #ifdef LOGGED
    //     double quat_grad[qdim] = {};
    //     double emap_grad[emapdim] = {};
    //     __enzyme_autodiff<void>((void *)emap_to_quat, enzyme_dup, quat, quat_grad,
    //                             enzyme_dup, emap, emap_grad);
    // #else
    //     emap_to_quat(quat, emap);
    // #endif

    double dquat_dexpmap_t[emapdim * qdim] = {};
#ifdef LOGGED
    double dquat_dexpmap_t_grad[emapdim * qdim];
    std::fill(dquat_dexpmap_t_grad, dquat_dexpmap_t_grad + emapdim * qdim, 1.0);
    double emap_grad[emapdim];
    std::fill(emap_grad, emap_grad + emapdim, 1.0);
    __enzyme_autodiff<void>((void *)dquat_demap_T, enzyme_dup, dquat_dexpmap_t,
                            dquat_dexpmap_t_grad, enzyme_dup, emap, emap_grad);
#else
    dquat_demap_T(dquat_dexpmap_t, emap);
#endif

    // Output the results for this test
    if (save_output) {
      *out_stream << "Test " << test + 1 << ":\n";

      *out_stream << "exponential map value:\n";
      for (size_t ie = 0; ie < emapdim; ie++) {
        *out_stream << emap[ie] << " ";
      }
      *out_stream << "\n";

      //   *out_stream << "unit quaternion value:\n";
      //   for (size_t iq = 0; iq < qdim; iq++) {
      //     *out_stream << quat[iq] << " ";
      //   }
      //   *out_stream << "\n";

      *out_stream << "dquat_dexpmap ^ T value:\n";
      for (size_t ie = 0; ie < emapdim; ie++) {
        for (size_t iq = 0; iq < qdim; iq++) {
          *out_stream << dquat_dexpmap_t[ECMECH_NM_INDX(ie, iq, qdim)] << " ";
        }
        *out_stream << "\n";
      }
      *out_stream << "----------------------------------------\n";
    }
  }

  auto end_time = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed = end_time - start_time;

  std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
            << std::scientific << "Elapsed time = " << elapsed.count()
            << " (s)\n";

  if (save_output) {
    outfile.close();
    std::cout << "Results saved to: " << output_path << std::endl;
  }

#ifdef LOGGED
  printLogger();
  destroyLogger();
#endif
  return 0;
}
