#include <cmath>
#include <cstdio>
#include "sv_interpolate.hpp"

/// @brief A slightly modified Neville's interpolation algorithm
/// @see 3.2 Polynomial Interpolation and Extrapolation, Numerical Recipes,
///      3rd Edition
/// @param[in]  x The point to interpolate at
/// @param[out] y Value of interpolating polynomial at x
/// @param[out] dy Error indication for computed y value
/// @param[in] xx The x-axis data points of size array_size
/// @param[in] yy The y-axis data points (values at xx's) of size array_size
/// @param[in] array_size The size of xx and yy arrays
/// @param[in] mm Number of data points to use for the interpolation; the
///            interval will be xx[from_index,...,from_index+mm-1]
/// @param[in] cws (optional) An array of size mm; if not given it will be
///            allocated/freed during function execution
/// @param[in] dws (optional) An array of size mm; if not given it will be
///            allocated/freed during function execution
/// @return Always returns an int, following the convention:
///         0: success
///         1: not enough points to perform interpolation
///        >1: computation error
int dso::sp3::neville_interpolation(double x, double &y, double &dy,
                          const double *__restrict__ xx, const double *__restrict__ yy,
                          int array_size, int mm, int from_index, double *cws,
                          double *dws) noexcept {

  const double *xpts = xx + from_index;
  const double *ypts = yy + from_index;

  if (from_index + mm > array_size) {
    fprintf(stderr,
            "[ERROR] Not enough data points to perform interpolation "
            "(traceback: %s)\n",
            __func__);
    return 1;
  }

  // Allocate workspace
  double *c, *d;
  c = (cws == nullptr) ? new double[mm] : cws;
  d = (dws == nullptr) ? new double[mm] : dws;

  int ns = 0;
  double difx;
  double dif = std::abs(x - xpts[0]);
  // ﬁnd the index ns of the closest table entry and initialize the tableau of
  // c’s and d’s
  for (int i = 0; i < mm; i++) {
    if ((difx = std::abs(x - xpts[i])) < dif) {
      ns = i;
      dif = difx;
    }
    c[i] = d[i] = ypts[i];
  }

  // initial approximation for y
  y = ypts[ns--];

  double ho, hp, w, den;
  // For each column of the tableau we loop over the current c’s and d’s and
  // update
  for (int m = 1; m < mm; m++) {
    for (int i = 0; i < mm - m; i++) {
      ho = xpts[i] - x;
      hp = xpts[i + m] - x;
      w = c[i + 1] - d[i];
      // This error can occur only if two input xa’s are(to within roundoﬀ)
      // identical
      if ((den = ho - hp) == 0e0) {
        fprintf(
            stderr,
            "[ERROR] x-axis points too close to interpolate!(traceback: %s)\n",
            __func__);
        // do not forget to free memmory if needed ....
        if (cws == nullptr)
          delete[] c;
        if (dws == nullptr)
          delete[] d;
      }
      den = w / den;
      d[i] = hp * den;
      c[i] = ho * den;
    }
    // After each column in the tableau is completed, we decide which
    // correction, c or d, we want to add to our accumulating value of y, i.e.,
    // which path to take through the tableau -- forking up or down --. We do
    // this in such a way as to take the most "straight line" route through the
    // tableau to its apex, updating ns accordingly to keep track of where we
    // are. This route keeps the partial approximations centered (insofar as
    // possible) on the target x.The last dy added is thus the error indication.
    y += (dy = (2 * (ns + 1) < (mm - m) ? c[ns + 1] : d[ns--]));
  }
  return 0;
}