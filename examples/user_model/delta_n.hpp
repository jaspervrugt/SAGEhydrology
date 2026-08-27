#pragma once
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

namespace loss_delta {

// ------------------------ utilities ------------------------

// Median in O(n) average using nth_element (copies input)
inline double median_copy(std::vector<double> v)
{
    const size_t n = v.size();
    if (n == 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    size_t mid = n / 2;
    std::nth_element(v.begin(), v.begin() + mid, v.end());
    double med = v[mid];
    if ((n % 2) == 0) {
        auto it = std::max_element(v.begin(), v.begin() + mid);
        med = 0.5 * (med + *it);
    }
    return med;
}

// MAD (median absolute deviation) with flag=1 behavior in MATLAB: median(|x - median(x)|)
inline double mad1(const double* x, int n)
{
    std::vector<double> v((size_t)n);
    for (int i = 0; i < n; ++i) {
        v[(size_t)i] = x[i];
    }
    const double med = median_copy(v);
    for (int i = 0; i < n; ++i) {
        v[(size_t)i] = std::fabs(x[i] - med);
    }
    return median_copy(v);
}

// Approximation to norminv(p) for p in (0,1): Peter J. Acklam inverse normal CDF
inline double norminv_acklam(double p)
{
    static const double a[] = {-3.969683028665376e+01,
                               2.209460984245205e+02,
                               -2.759285104469687e+02,
                               1.383577518672690e+02,
                               -3.066479806614716e+01,
                               2.506628277459239e+00};
    static const double b[] = {-5.447609879822406e+01,
                               1.615858368580409e+02,
                               -1.556989798598866e+02,
                               6.680131188771972e+01,
                               -1.328068155288572e+01};
    static const double c[] = {-7.784894002430293e-03,
                               -3.223964580411365e-01,
                               -2.400758277161838e+00,
                               -2.549732539343734e+00,
                               4.374664141464968e+00,
                               2.938163982698783e+00};
    static const double d[] = {7.784695709041462e-03,
                               3.224671290700398e-01,
                               2.445134137142996e+00,
                               3.754408661907416e+00};

    const double plow = 0.02425;
    const double phigh = 1.0 - plow;

    if (p <= 0.0) {
        return -std::numeric_limits<double>::infinity();
    }
    if (p >= 1.0) {
        return std::numeric_limits<double>::infinity();
    }

    double q, r;
    if (p < plow) {
        q = std::sqrt(-2.0 * std::log(p));
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
               ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
    } else if (p > phigh) {
        q = std::sqrt(-2.0 * std::log(1.0 - p));
        return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
               ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
    } else {
        q = p - 0.5;
        r = q * q;
        return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
               (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0);
    }
}

static inline double sgn(double x)
{
    return (x > 0.0) - (x < 0.0);
}

// Quantile/W1 FDC loss and delta (O(n log n))
static void fdc_w1_loss_delta(const double* y,
                              const double* q,
                              int n,
                              double* loss_out,
                              std::vector<double>& delta_out)
{
    delta_out.assign((size_t)n, 0.0);

    // sort q, but keep permutation to unsort delta
    std::vector<int> P((size_t)n);
    for (int i = 0; i < n; ++i) {
        P[(size_t)i] = i;
    }

    std::sort(P.begin(), P.end(), [&](int a, int b) { return q[a] < q[b]; });

    std::vector<double> y_sorted((size_t)n);
    for (int i = 0; i < n; ++i) {
        y_sorted[(size_t)i] = y[i];
    }
    std::sort(y_sorted.begin(), y_sorted.end());

    double L = 0.0;
    const double invn = 1.0 / (double)n;

    for (int i = 0; i < n; ++i) {
        const int idx = P[(size_t)i];
        const double dq = q[idx] - y_sorted[(size_t)i];
        L += std::fabs(dq) * invn;
        // delta in sorted space
        const double dLdq_sorted = invn * sgn(dq); // d/dq_sorted
        // unsort back to original time order
        delta_out[(size_t)idx] = dLdq_sorted;
    }

    if (loss_out) {
        *loss_out = L;
    }
}

// ------------------------ API ------------------------
// Sigma_eps handling for L=2:
// - sigma_scalar != nullptr => Sigma_eps = sigma_scalar * I
// - sigma_diag   != nullptr => Sigma_eps = diag(sigma_diag) length n
// If both null => identity
inline void delta_n(int L,
                    const double* y,
                    const double* q,
                    int n,
                    std::vector<double>& delta,
                    const double* sigma_scalar = nullptr,
                    const double* sigma_diag = nullptr)
{
    delta.assign((size_t)n, 0.0);

    // residual e = y - q
    std::vector<double> e((size_t)n);
    for (int i = 0; i < n; ++i) {
        e[(size_t)i] = y[i] - q[i];
    }

    switch (L) {
    case 1: { // Sum |res|
        for (int i = 0; i < n; ++i) {
            delta[(size_t)i] = -sgn(e[(size_t)i]);
        }
    } break;

    case 2: { // GLS: delta = -Sigma^{-1} e
        if (sigma_scalar) {
            const double s = *sigma_scalar;
            const double invs = (s != 0.0) ? (1.0 / s) : 0.0;
            for (int i = 0; i < n; ++i) {
                delta[(size_t)i] = -invs * e[(size_t)i];
            }
        } else if (sigma_diag) {
            for (int i = 0; i < n; ++i) {
                const double d = sigma_diag[i];
                delta[(size_t)i] = (d != 0.0) ? -(e[(size_t)i] / d) : 0.0;
            }
        } else {
            for (int i = 0; i < n; ++i) {
                delta[(size_t)i] = -e[(size_t)i];
            }
        }
    } break;

    case 3: { // NSE: delta = -(2/SSt) * e
        double my = 0.0;
        for (int i = 0; i < n; ++i) {
            my += y[i];
        }
        my /= (double)n;

        double SSt = 0.0;
        for (int i = 0; i < n; ++i) {
            const double dy = y[i] - my;
            SSt += dy * dy;
        }
        const double c = (SSt > 0.0) ? (-2.0 / SSt) : 0.0;
        for (int i = 0; i < n; ++i) {
            delta[(size_t)i] = c * e[(size_t)i];
        }
    } break;

    case 4: { // KGE
        double my = 0.0, mq = 0.0;
        for (int i = 0; i < n; ++i) {
            my += y[i];
            mq += q[i];
        }
        my /= (double)n;
        mq /= (double)n;

        double sYY = 0.0, sQQ = 0.0, sYQ = 0.0;
        for (int i = 0; i < n; ++i) {
            const double dy = y[i] - my;
            const double dq = q[i] - mq;
            sYY += dy * dy;
            sQQ += dq * dq;
            sYQ += dy * dq;
        }
        const double denom = (n > 1) ? (double)(n - 1) : 1.0;
        const double varY = sYY / denom;
        const double varQ = sQQ / denom;
        const double covYQ = sYQ / denom;

        const double sy = (varY > 0.0) ? std::sqrt(varY) : 0.0;
        const double sq = (varQ > 0.0) ? std::sqrt(varQ) : 0.0;

        const double r = (sy > 0.0 && sq > 0.0) ? (covYQ / (sy * sq)) : 0.0;
        const double v = (sy > 0.0) ? (sq / sy) : 0.0;
        const double z = (my != 0.0) ? (mq / my) : 0.0;

        const double Lkge = std::sqrt((r - 1.0) * (r - 1.0) + (v - 1.0) * (v - 1.0) +
                                      (z - 1.0) * (z - 1.0));
        if (!(Lkge > 0.0) || !std::isfinite(Lkge)) {
            std::fill(delta.begin(), delta.end(), 0.0);
            break;
        }

        const double invn = 1.0 / (double)n;
        const double invnm1 = (n > 1) ? (1.0 / (double)(n - 1)) : 0.0;
        const double invsy = (sy > 0.0) ? (1.0 / sy) : 0.0;
        const double invsq = (sq > 0.0) ? (1.0 / sq) : 0.0;
        const double invsq2 = (sq > 0.0) ? (1.0 / (sq * sq)) : 0.0;

        for (int i = 0; i < n; ++i) {
            const double qcen = q[i] - mq;
            const double ycen = y[i] - my;

            const double dmqdqt = invn;
            const double dsqdqt = (sq > 0.0) ? (qcen * invnm1 * invsq) : 0.0;

            const double drdqt =
                (sy > 0.0 && sq > 0.0)
                    ? (ycen * invnm1 * invsq * invsy)
                    : 0.0 - r * ((sq > 0.0) ? (qcen * invnm1 * invsq2) : 0.0);

            const double dvdqt = dsqdqt * invsy;
            const double dzdqt = (my != 0.0) ? (dmqdqt / my) : 0.0;

            delta[(size_t)i] =
                ((r - 1.0) * drdqt + (v - 1.0) * dvdqt + (z - 1.0) * dzdqt) / Lkge;
        }
    } break;

    case 5: { // Huber
        const double kappa = 1.0 / norminv_acklam(0.75);
        const double S = kappa * mad1(y, n);
        const double invS = (S > 0.0) ? (1.0 / S) : 0.0;
        const double c = 1.345;

        for (int i = 0; i < n; ++i) {
            const double u = e[(size_t)i] * invS;
            double psi = (std::fabs(u) <= c) ? u : (c * sgn(u));
            delta[(size_t)i] = -invS * psi;
        }
    } break;

    case 6: { // FDC pairwise O(n^2) -- WARNING for large n
        const double invn2 = 1.0 / ((double)n * (double)n);
        for (int i = 0; i < n; ++i) {
            double t1 = 0.0, t2 = 0.0;
            const double qi = q[i];
            for (int j = 0; j < n; ++j) {
                t1 += sgn(qi - y[j]);
                t2 += sgn(qi - q[j]);
            }
            delta[(size_t)i] = (t1 - t2) * invn2;
        }
    } break;

    default:
        std::fill(delta.begin(), delta.end(), 0.0);
        break;
    }
}

} // namespace loss_delta