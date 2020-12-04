#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#define MATHLIB_STANDALONE 1
#define xmalloc(nbytes) malloc_or_exit(nbytes, __FILE__, __LINE__)

void *malloc_or_exit(size_t nbytes, const char * file, int line) {
  void *p = malloc(nbytes);
  if (p == NULL) {
    fprintf(stderr, "%s, line %d:"
      "unable to allocate %lu bytes, calling exit()\n", file, line,
      ((unsigned long) nbytes));
  }
  return p;
}

// define vector
double *dvector(int n) {
  return xmalloc(n * sizeof(double));
}

// define function to compute unbiased estimate of tr(C_suhk C^T_svhl)
// the unbiased estimator is defined in the text below paragraph 5 of
// Zhong, Li, Santo (2019)
double compute_tr_chat(double *Yvec, int n, int p, int T, int s, int s1, int h, int h1) {
  int i, j, k, l, d;
  double temp, temp1, temp2, temp3, temp4, Yjh1, Yks1, Yis1, Ykh1, Yls1;
  double Ys = 0.0, Yh = 0.0, Ys1 = 0.0, Yh1 = 0.0;
  double tr_cov = 0.0, tr_cov0 = 0.0, tr_cov1 = 0.0, tr_cov2 = 0.0, tr_cov3 = 0.0;

  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if (i != j) {
        temp = 0.0;
        temp1 = 0.0;
        for (d = 0; d < p; d++) {
          Ys = Yvec[(unsigned long)(s - 1) * p * n + i * p + d];
          Yh = Yvec[(unsigned long)(h - 1) * p * n + j * p + d];
          Ys1 = Yvec[(unsigned long)(s1 - 1) * p * n + i * p + d];
          Yh1 = Yvec[(unsigned long)(h1 - 1) * p * n + j * p + d];
          temp += Ys * Yh;
          temp1 += Ys1 * Yh1;
        }
				// term 1
        tr_cov += temp * temp1 / (((double) n) * (n - 1.0));
        for (k = 0; k < n; k++) {
          if (k != i && k != j) {
            temp2 = 0.0;
            temp3 = 0.0;
            for (d = 0; d < p; d++) {
              Yjh1 = Yvec[(unsigned long)(h1 - 1) * p * n + j * p + d];
              Yks1 = Yvec[(unsigned long)(s1 - 1) * p * n + k * p + d];
              temp2 += Yjh1 * Yks1;
              Yis1 = Yvec[(unsigned long)(s1 - 1) * p * n + i * p + d];
              Ykh1 = Yvec[(unsigned long)(h1 - 1) * p * n + k * p + d];
              temp3 += Yis1 * Ykh1;
            }
						// terms 2 and 3
            tr_cov1 += temp * temp2 / ((double) n * (n - 1.0) * (n - 2.0));
            tr_cov2 += temp * temp3 / ((double) n * (n - 1.0) * (n - 2.0));
            for (l = 0; l < n; l++) {
              if (l != i && l != j && l != k) {
                temp4 = 0.0;
                for (d = 0; d < p; d++) {
                  Ykh1 = Yvec[(unsigned long)(h1 - 1) * p * n + k * p + d];
                  Yls1 = Yvec[(unsigned long)(s1 - 1) * p * n + l * p + d];
                  temp4 += Ykh1 * Yls1;
                }
								// term 4
                tr_cov3 += temp * temp4 / ((double) n * (n - 1.0) * (n - 2.0) * (n - 3.0));
              }
            }
          }
        }
      }
    }
  }

  tr_cov0 = tr_cov - tr_cov1 - tr_cov2 + tr_cov3;
  return (tr_cov0);
}

// define function to compute the variance for a given k1 and k2
// given k1 and k2, this function computes Q_hat defined below Theorem 3
// of Zhong, Li, Santo (2019)
double compute_variance(double *tr_covvec, int n, int T, int k1, int k2) {

  int s, s1, h, h1;
  double tr0, tr1, tr2, tr3, tr4, tr5, tr6, tr7, tr8, tr9, varvec = 0.0;

  for (s = 1; s < (k1 + 1); s++) {
    for (s1 = 1; s1 < (k2 + 1); s1++) {
      tr0 = tr_covvec[(s - 1) + T * (s1 - 1) + (unsigned long) T * T * (s - 1) + (unsigned long) T * T * T * (s1 - 1)];

      for (h = (k1 + 1); h < (T + 1); h++) {
        tr2 = tr_covvec[(h - 1) + T * (s1 - 1) + (unsigned long) T * T * (h - 1) + (unsigned long) T * T * T * (s1 - 1)];
        tr6 = tr_covvec[(s - 1) + T * (s1 - 1) + (unsigned long) T * T * (h - 1) + (unsigned long) T * T * T * (s1 - 1)];

        for (h1 = (k2 + 1); h1 < (T + 1); h1++) {
          tr1 = tr_covvec[(s - 1) + T * (h1 - 1) + (unsigned long) T * T * (s - 1) + (unsigned long) T * T * T * (h1 - 1)];
          tr3 = tr_covvec[(h - 1) + T * (h1 - 1) + (unsigned long) T * T * (h - 1) + (unsigned long) T * T * T * (h1 - 1)];
          tr4 = tr_covvec[(s - 1) + T * (s1 - 1) + (unsigned long) T * T * (s - 1) + (unsigned long) T * T * T * (h1 - 1)];
          tr5 = tr_covvec[(h - 1) + T * (s1 - 1) + (unsigned long) T * T * (h - 1) + (unsigned long) T * T * T * (h1 - 1)];
          tr7 = tr_covvec[(s - 1) + T * (h1 - 1) + (unsigned long) T * T * (h - 1) + (unsigned long) T * T * T * (h1 - 1)];
          tr8 = tr_covvec[(s - 1) + T * (s1 - 1) + (unsigned long) T * T * (h - 1) + (unsigned long) T * T * T * (h1 - 1)];
          tr9 = tr_covvec[(s - 1) + T * (h1 - 1) + (unsigned long) T * T * (h - 1) + (unsigned long) T * T * T * (s1 - 1)];

          varvec += (tr0 * tr0 + tr1 * tr1 + tr2 * tr2 + tr3 * tr3 - 2 * tr4
						* tr4 - 2 * tr5 * tr5 - 2 * tr6 * tr6 - 2 * tr7 * tr7 + 2 * tr8
						* tr8 + 2 * tr9 * tr9);
        }
      }
    }
  }

  varvec = 4 * varvec / ((double) n * (n - 1.0) * k1 * k2 * (T - k1) * (T - k2));
  return (varvec);
}

// define function to compute dk hat
// quantity d-hat is defined below expression 4 in Zhong, Li, Santo (2019),
// here d-hat is computed in a recursive manner to improve the computation
// time
int compute_dk_hat(double *Yvec, int n, int p, int T, double *dk, int *khat) {
  int k, j, t, s, s0, t0, k0, h, h1;
  double *tr_sigk = dvector(T), temp_ck, temp_ck1, temp_ck2;
	double ck, ck1, ck2, mk, mTk, largest;

  for (k = 0; k < T; k++) {
    j = k + 1;
    tr_sigk[k] = compute_tr_chat(Yvec, n, p, T, j, j, j, j);
  }

  for (k = 0; k < 1; k++) {
    ck = 0;
    k0 = k + 1;
    for (t = 0; t < k + 1; t++) {
      for (s = 0; s < T - k0; s++) {
        s0 = s + k0 + 1;
        t0 = t + 1;
        temp_ck = compute_tr_chat(Yvec, n, p, T, t0, t0, s0, s0);
        ck += temp_ck;
      }
    }
    mk = 0.0;
    for (h = 0; h < k + 1; h++) {
      mk += tr_sigk[h] / (k + 1.0);
    }
    mTk = 0.0;
    for (h1 = 0; h1 < T - k0; h1++) {
      mTk += tr_sigk[h1 + k0] / (T - k0);
    }
    dk[k] = mk + mTk - 2 * ck / (k0 * (T - k0) + 0.0);
  }

  for (k = 1; k < T - 1; k++) {
    ck1 = 0;
    ck2 = 0;
    k0 = k + 1;

    mk = 0.0;
    for (h = 0; h < k; h++) {
      mk += tr_sigk[h];
      temp_ck1 = compute_tr_chat(Yvec, n, p, T, k + 1, k + 1, h + 1, h + 1);
      ck1 += temp_ck1;
    }
    mTk = 0.0;
    for (h1 = 0; h1 < T - k0; h1++) {
      mTk += tr_sigk[h1 + k0];
      temp_ck2 = compute_tr_chat(Yvec, n, p, T, h1 + k0 + 1, h1 + k0 + 1, k + 1, k + 1);
      ck2 += temp_ck2;
    }

    dk[k] = ((k0 - 1) * (T - (k0 - 1)) * dk[k - 1] + (T - 2 * k0 + 1)
		* tr_sigk[k] - mk + mTk + 2 * ck1 - 2 * ck2) / (k0 * (T - k0));
  }

  khat[0] = 1;
  largest = dk[0];
  if (T > 2) {
    for (k = 0; k < T - 2; k++) {
      k0 = k + 1;
      if (dk[k0] > largest) {
        khat[0] = k0 + 1;
        largest = dk[k0];
      }
    }
  }

  free(tr_sigk);
  return 0;
}

// perform steps for testing and change point identification
int compute_test_cp(double *Yvec, int n, int p, int T, int *khat, double *std_dk, double *max_std_dk, double *corr_mat) {
  int T0 = T - 1, k1, k2, nk2, s, s1, h, h1;
  double *tr_covvec = dvector((unsigned long) T * T * T * T);
	double *dk = dvector(T - 1);
	double *cov_mat = dvector(T0 * (T0 + 1) / 2), varvec, tr_vecval;

  compute_dk_hat(Yvec, n, p, T, dk, khat);

	// compute tr(C_suhk C^T_svhl) for all possible time points
	// store in a vector to make computation of Q-hat fast in the next step
  for (s = 1; s < (T + 1); s++) {
    for (s1 = 1; s1 < (T + 1); s1++) {
      for (h = 1; h < (T + 1); h++) {
        for (h1 = 1; h1 < (T + 1); h1++) {
          tr_vecval = compute_tr_chat(Yvec, n, p, T, s, s1, h, h1);
          tr_covvec[(h1 - 1) + T * (h - 1) + (unsigned long) T * T * (s1 - 1) + (unsigned long) T * T * T * (s - 1)] = tr_vecval;
        }
      }
    }
  }

	// compute lower triangle of Q-hat
  for (k1 = 0; k1 < T0; k1++) {
    for (k2 = 0; k2 < T0 - k1; k2++) {
      nk2 = k2 + k1;
      varvec = compute_variance(tr_covvec, n, T, k1 + 1, nk2 + 1);
      cov_mat[(2 * T0 - k1 + 1) * k1 / 2 + k2] = varvec;
    }
  }

	// compute Mn as defined in equation 9 of Zhong, Li, Santo (2019)
	// maximum standardized value of D-hat across T - 1 time points
  max_std_dk[0] = dk[0] / sqrt(cov_mat[0]);
  for (k1 = 0; k1 < T0; k1++) {
    std_dk[k1] = dk[k1] / sqrt(cov_mat[(2 * T0 - k1 + 1) * k1 / 2]);
    if (std_dk[k1] > max_std_dk[0]) {
      max_std_dk[0] = std_dk[k1];
    }
  }

	// compute correlation matrix, defined as V-hat below Theorem 3 in
	// Zhong, Li, Santo (2019)
  for (k1 = 0; k1 < T0; k1++) {
    for (k2 = 0; k2 < T0 - k1; k2++) {
      corr_mat[(2 * T0 - k1 + 1) * k1 / 2 + k2] = cov_mat[(2 * T0 - k1 + 1) * k1 / 2 + k2] /
        sqrt(cov_mat[(2 * T0 - k1 + 1) * k1 / 2] * cov_mat[(2 * T0 - k1 - k2 + 1) * (k1 + k2) / 2]);
    }
  }

  free(dk);
  free(tr_covvec);
  free(cov_mat);
  return 0;
}

// function for linking R and C code
void compute_test_cp_r(double *Yvec, int *n, int *p, int *T, int *khat, double *std_dk, double *max_std_dk, double *corr_mat) {
  compute_test_cp(Yvec, n[0], p[0], T[0], khat, std_dk, max_std_dk, corr_mat);
}
