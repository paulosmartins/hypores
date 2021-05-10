#include "HyPoRes.h"
#include <NTL/mat_ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/LLL.h>
#include <assert.h>

/* precomputations */
std::array<value_t, h2> b2i_over_B2;
std::array<std::array<value_t, n>, h1> m_1;
std::array<std::array<value_t, n>, h1> beta_m_1;
std::array<std::array<value_t, h1>, h2> B1_over_b1j_mod_b2i;
std::array<value_t, h1> B1_over_b1i_sk;
std::array<value_t, h2> B2_over_b2i_B1;
std::array<std::array<value_t, n>, h2> m_2;
std::array<std::array<value_t, n>, h2> beta_m_2;
std::array<value_t, h2> b2i_inv_sk;
value_t B1_inv_sk;
value_t B2_inv_sk;
std::array<value_t, n> m_sk;
std::array<value_t, n> beta_m_sk;
std::array<std::array<value_t, h2>, h1> B2_over_b2i_mod_b1i;
std::array<value_t, h1> B2_mod_b1i;

/* channel arithmetic */
std::array<value_t, h1> c1;
std::array<value_t, h2> c2;
value_t mask;
value_t mask_sk;
  
value_t reduce (greater_value_t x, value_t c, value_t p);
value_t addmod (greater_value_t a, value_t b, value_t p);
value_t submod (value_t a, value_t b, value_t p);
value_t mulmod (greater_value_t a, value_t b, value_t c, value_t p);
value_t mulbeta (greater_value_t a, value_t c, value_t p);
value_t mulbeta_sk (value_t a);

template<size_t h>
std::array<std::array<value_t, n>, h>
rns_pol (const std::array<NTL::ZZ, n> &x,
	 const std::array<value_t, h> &ps)
{
  std::array<std::array<value_t, n>, h> x_red;

  for (size_t i = 0; i < h; i++) {
    NTL::ZZ pi (ps[i]);
    for (size_t j = 0; j < n; j++) {
      NTL::ZZ x_red_ij = x[j] % pi;
      x_red[i][j] = NTL::conv<value_t> (x_red_ij);
    }
  }

  return x_red;
}

template<size_t h>
std::array<NTL::ZZ, n>
rns_pol_inv (const std::array<std::array<value_t, n>, h> &x1,
	     const std::array<value_t, h> &ps)
{
  std::array<NTL::ZZ, n> X;
  std::array<NTL::ZZ, n> P;

  for (size_t i = 0; i < n; i++) {
    if (x1[0][i] > (ps[0]/2)) {
      X[i] = NTL::ZZ (x1[0][i]) - NTL::ZZ (ps[0]);
    } else {
      X[i] = NTL::ZZ (x1[0][i]);
    }
    P[i] = ps[0];
  }

  for (size_t i = 1; i < h; i++) {
    for (size_t j = 0; j < n; j++) {
      NTL::CRT (X[j], P[j], NTL::ZZ (x1[i][j]), NTL::ZZ (ps[i]));
    }
  }

  return X;
}

HPolRNS::HPolRNS ()
{
  /* channel arithmetic */
  for (size_t i = 0; i < h1; i++) {
    c1[i] = (((greater_value_t)1)<<logw) - b1[i];
  }
  for (size_t i = 0; i < h2; i++) {
    c2[i] = (((greater_value_t)1)<<logw) - b2[i];
  }

  mask = (((greater_value_t)1)<<logw) - 1;
  mask_sk = (((greater_value_t)1)<<log_bsk) - 1;

  /* precomputations */
  NTL::ZZ psk (1);
  psk <<= log_bsk;

  NTL::ZZ B2 (1);
  for (size_t i = 0; i < h2; i++) {
    B2 *= NTL::ZZ (b2[i]);
  }
  NTL::ZZ B1 (1);
  for (size_t i = 0; i < h1; i++) {
    B1 *= NTL::ZZ (b1[i]);
  }
  
  for (size_t i = 0; i < h2; i++) {
    NTL::ZZ B2_over_b2i (B2 / NTL::ZZ (b2[i]));
    b2i_over_B2[i] = NTL::conv<value_t>
      (NTL::InvMod (B2_over_b2i % NTL::ZZ (b2[i]), NTL::ZZ (b2[i])));
    value_t B1inv = NTL::conv<value_t>
      (NTL::InvMod (B1 % NTL::ZZ (b2[i]), NTL::ZZ (b2[i])));
    B2_over_b2i_B1[i] = NTL::conv<value_t>
      ((B2_over_b2i * NTL::ZZ (B1inv)) % NTL::ZZ (b2[i]));

    for (size_t j = 0; j < n; j++) {
      m_2[i][j] = mulmod (rns_m_2[i][j], b2i_over_B2[i], c2[i], b2[i]);
      m_2[i][j] = mulmod (m_2[i][j], B1inv, c2[i], b2[i]);
      beta_m_2[i][j] = mulbeta (m_2[i][j], c2[i], b2[i]);
    }

    for (size_t j = 0; j < h1; j++) {
      B2_over_b2i_mod_b1i[j][i] = NTL::conv<value_t>
	(B2_over_b2i % NTL::ZZ (b1[j]));
    }
  }
  
  for (size_t i = 0; i < h1; i++) {
    NTL::ZZ B1_over_b1i (B1 / NTL::ZZ (b1[i]));
    value_t b1i_over_B1 = NTL::conv<value_t>
      (NTL::InvMod (B1_over_b1i % NTL::ZZ (b1[i]), NTL::ZZ (b1[i])));
    
    for (size_t j = 0; j < n; j++) {
      m_1[i][j] = mulmod (rns_minv_m_1[i][j], b1i_over_B1, c1[i], b1[i]);
      beta_m_1[i][j] = mulbeta (m_1[i][j], c1[i], b1[i]);
    }

    for (size_t j = 0; j < h2; j++) {
      B1_over_b1j_mod_b2i[j][i] = NTL::conv<value_t>
	(B1_over_b1i % NTL::ZZ (b2[j]));
    }

    B1_over_b1i_sk[i] = NTL::conv<value_t>
      (B1_over_b1i % psk);

    B2_mod_b1i[i] = NTL::conv<value_t>
      (B2 % NTL::ZZ (b1[i]));
  }

  for (size_t i = 0; i < h2; i++) {
    b2i_inv_sk[i] = NTL::conv<value_t>
      (NTL::InvMod (NTL::ZZ (b2[i]) % psk, psk));
  }
  B1_inv_sk = NTL::conv<value_t>
    (NTL::InvMod (B1 % psk, psk));
  B2_inv_sk = NTL::conv<value_t>
    (NTL::InvMod (B2 % psk, psk));
  for (size_t i = 0; i < n; i++) {
    m_sk[i] = rns_m_sk[i] * B1_inv_sk;
    m_sk[i] &= mask_sk;
    beta_m_sk[i] = mulbeta_sk (m_sk[i]);
  }
}

value_t
reduce (greater_value_t x, value_t c, value_t p)
{
  const greater_value_t bound = ((greater_value_t)1)<<logw;
  
  while (x >= bound) {
    greater_value_t x0 = x & mask;
    greater_value_t x1 = x >> logw;
    x = x0 + c * x1;
  }

  if (x >= p) x -= p;

  return x;
}

value_t
addmod (greater_value_t a, value_t b, value_t p)
{
  greater_value_t c = a + b;
  if (c >= p) c -= p;
  return c;
}

value_t
submod (value_t a, value_t b, value_t p)
{
  return (b == 0 ? a : addmod (a, p-b, p));
}

value_t
mulmod (greater_value_t a, value_t b, value_t c, value_t p)
{
  greater_value_t d = a*b;
  return reduce (d, c, p);
}

template<typename T, size_t beta1, size_t bit, size_t i>
struct _mul_beta { };

template<typename T, size_t i>
struct _mul_beta<T, 0, 0, i>
{
  T apply (const T &x)
  {
    return 0;
  }
};

template<typename T, size_t beta1, size_t i>
struct _mul_beta<T, beta1, 0, i>
{
  T apply (const T &x)
  {
    return _mul_beta<T, (beta1>>1), (beta1&1), i+1>{}.apply (x);
  }
};

template<typename T, size_t beta1, size_t i>
struct _mul_beta<T, beta1, 1, i>
{
  T apply (const T &x)
  {
    return (x<<i) +
      _mul_beta<T, (beta1>>1), (beta1&1), i+1>{}.apply (x);
  }
};

value_t
mulbeta (greater_value_t a, value_t c, value_t p)
{
  greater_value_t d = _mul_beta<greater_value_t, (beta>>1), beta&1, 0>{}.apply (a);
  
  return reduce (d, c, p);
}

value_t
mulbeta_sk (value_t a)
{
  value_t d = _mul_beta<value_t, (beta>>1), beta&1, 0>{}.apply (a);
  
  return d & mask_sk;
}

HPolRNS_ZZq
HPolRNS::to_ZZq (const NTL::ZZ &a)
{
  NTL::mat_ZZ Bm;
  Bm.SetDims (n+1, n+1);

  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      if (j >= i) {
	Bm[i][j] = m[j-i];
      } else {
	Bm[i][j] = NTL::ZZ (beta) * m[n+j-i];
      }
    }
  }

  Bm[n][0] = a;
  Bm[n][n] = P;

  NTL::BKZ_RR (Bm);

  std::array<NTL::ZZ, n> x;
  for (size_t i = 0; i < n; i++) {
    x[i] = Bm[n][i];
  }

  auto a_1 = rns_pol (x, b1);
  auto a_2 = rns_pol (x, b2);

  for (size_t i = 0; i < h2; i++) {
    for (size_t j = 0; j < n; j++) {
      a_2[i][j] = mulmod (a_2[i][j], b2i_over_B2[i], c2[i], b2[i]);
    }
  }

  std::array<value_t, n> a_sk;
  NTL::ZZ psk (1);
  psk <<= log_bsk;
  for (size_t i = 0; i < n; i++) {
    a_sk[i] = NTL::conv<value_t> (x[i] % psk);
  }

  return HPolRNS_ZZq
    {
      a_1, a_2, a_sk
    };
}

NTL::ZZ
HPolRNS::to_ZZ (const HPolRNS_ZZq &A)
{
  auto as = rns_pol_inv (A.a_1, b1);
  NTL::ZZ x (0);
  NTL::ZZ gamma_i (1);
  for (size_t i = 0; i < n; i++) {
    x = (x + as[i] * gamma_i) % P;
    gamma_i = (gamma_i * hpr_gamma) % P;
  }

  return x;
}

template<size_t h>
void rns_pol_mul
(std::array<std::array<value_t, n>, h> &zs,
 const std::array<std::array<value_t, n>, h> &xs,
 const std::array<std::array<value_t, n>, h> &ys,
 const std::array<std::array<value_t, n>, h> &ys_beta,
 const std::array<value_t, h> &cs,
 const std::array<value_t, h> &ps)
{
  for (size_t i = 0; i < h; i++) {
    for (size_t j = 0; j < n; j++) {
      zs[i][j] = 0;
      
      for (size_t k = 0; k <= j; k++) {
	value_t prod = mulmod (xs[i][k], ys[i][j-k], cs[i], ps[i]);
	zs[i][j] = addmod (zs[i][j], prod, ps[i]);
      }

      for (size_t k = j+1; k < n; k++) {
	value_t prod = mulmod (xs[i][k], ys_beta[i][n+j-k], cs[i], ps[i]);
	zs[i][j] = addmod (zs[i][j], prod, ps[i]);
      }
    }
  }
}

void rns_pol_mul_sk
(std::array<value_t, n> &zs,
 const std::array<value_t, n> &xs,
 const std::array<value_t, n> &ys,
 const std::array<value_t, n> &ys_beta)
{
  for (size_t j = 0; j < n; j++) {
    zs[j] = 0;
      
    for (size_t k = 0; k <= j; k++) {
      value_t prod = xs[k] * ys[j-k];
      zs[j] += prod;
    }

    for (size_t k = j+1; k < n; k++) {
      value_t prod = xs[k] * ys_beta[n+j-k];
      zs[j] += prod;
    }

    zs[j] &= mask_sk;
  }
}


template<size_t h>
void
rns_pol_mul_beta (std::array<std::array<value_t, n>, h> &zs,
		  const std::array<std::array<value_t, n>, h> &xs,
		  const std::array<std::array<value_t, n>, h> &ys,
		  const std::array<value_t, h> &cs,
		  const std::array<value_t, h> &ps)
{
  for (size_t i = 0; i < h; i++) {
    for (size_t j = 0; j < n; j++) {
      zs[i][j] = 0;

      for (size_t k = j+1; k < n; k++) {
	value_t prod = mulmod (xs[i][k], ys[i][n+j-k], cs[i], ps[i]);
	zs[i][j] = addmod (zs[i][j], prod, ps[i]);
      }
      zs[i][j] = mulbeta (zs[i][j], cs[i], ps[i]);

      for (size_t k = 0; k <= j; k++) {
	value_t prod = mulmod (xs[i][k], ys[i][j-k], cs[i], ps[i]);
	zs[i][j] = addmod (zs[i][j], prod, ps[i]);
      }
    }
  }
}

void
rns_pol_mul_sk_beta (std::array<value_t, n> &zs,
		     const std::array<value_t, n> &xs,
		     const std::array<value_t, n> &ys)
{
  for (size_t j = 0; j < n; j++) {
    zs[j] = 0;

    for (size_t k = j+1; k < n; k++) {
      value_t prod = xs[k] * ys[n+j-k];
      zs[j] += prod;
    }
    zs[j] = mulbeta_sk (zs[j]);

    for (size_t k = 0; k <= j; k++) {
      value_t prod = xs[k] * ys[j-k];
      zs[j] += prod;
    }

    zs[j] &= mask_sk;
  }
}

template<size_t h1, size_t h2>
void
rns_pol_basis_extension (std::array<std::array<value_t, n>, h2> &zs,
			 const std::array<std::array<value_t, n>, h1> &xs,
			 const std::array<std::array<value_t, h1>, h2> &bs,
			 const std::array<value_t, h2> cs2,
			 const std::array<value_t, h2> ps2,
			 const std::array<value_t, h1> ps1)
{
  for (size_t i = 0; i < h2; i++) {
    for (size_t j = 0; j < n; j++) {
      zs[i][j] = 0;
      
      for (size_t k = 0; k < h1; k++) {
	value_t pdiff = (ps1[k] > ps2[i] ? ps2[i]<<1 : ps2[i]) - ps1[k];
	value_t xskj = xs[k][j];
	bool negative = false;
	if (xskj > ps1[k]/2) {
	  negative = true;
	}
	if (xskj > ps2[i]) {
	  xskj -= ps2[i];
	}
	if (negative) {
	  xskj = addmod (xskj, pdiff, ps2[i]);
	}
	value_t prod = mulmod (xskj, bs[i][k], cs2[i], ps2[i]);
	zs[i][j] = addmod (zs[i][j], prod, ps2[i]);
      }
    }
  }
}
			 
template<size_t h1>
void
rns_pol_basis_extension_sk (std::array<value_t, n> &zs,
			    const std::array<std::array<value_t, n>, h1> &xs,
			    const std::array<value_t, h1> &bs,
			    const std::array<value_t, h1> &ps1)
{
  for (size_t j = 0; j < n; j++) {
    zs[j] = 0;
      
    for (size_t k = 0; k < h1; k++) {
      value_t xskj = xs[k][j];
      if (xskj > ps1[k]/2) {
	xskj -= ps1[k];
      }
      value_t prod = xskj * bs[k];
      zs[j] += prod;
    }

    zs[j] &= mask_sk;
  }
}

template<size_t h>
void
rns_pol_scale (std::array<std::array<value_t, n>, h> &zs,
	       const std::array<std::array<value_t, n>, h> &xs,
	       const std::array<value_t, h> &ys,
	       const std::array<value_t, h> &cs,
	       const std::array<value_t, h> &ps)
{
  for (size_t i = 0; i < h; i++) {
    for (size_t j = 0; j < n; j++) {
      zs[i][j] = mulmod (xs[i][j], ys[i], cs[i], ps[i]);
    }
  }
}

void
rns_pol_scale_sk (std::array<value_t, n> &zs,
		  const std::array<value_t, n> &xs,
		  value_t y)
{
  for (size_t i = 0; i < n; i++) {
    zs[i] = xs[i] * y;
    zs[i] &= mask_sk;
  }
}

template<size_t h>
void
rns_pol_add (std::array<std::array<value_t, n>, h> &zs,
	     const std::array<std::array<value_t, n>, h> &xs,
	     const std::array<std::array<value_t, n>, h> &ys,
	     const std::array<value_t, h> &ps)
{
  for (size_t i = 0; i < h; i++) {
    for (size_t j = 0; j < n; j++) {
      zs[i][j] = addmod (xs[i][j], ys[i][j], ps[i]);
    }
  }
}

void
rns_pol_add_sk (std::array<value_t, n> &zs,
		const std::array<value_t, n> &xs,
		const std::array<value_t, n> &ys)
{
  for (size_t j = 0; j < n; j++) {
    zs[j] = xs[j] + ys[j];
    zs[j] &= mask_sk;
  }
}

void
rns_pol_sub_sk (std::array<value_t, n> &zs,
		const std::array<value_t, n> &xs,
		const std::array<value_t, n> &ys)
{
  for (size_t j = 0; j < n; j++) {
    zs[j] = xs[j] - ys[j];
    zs[j] &= mask_sk;
  }
}

void
rns_pol_correct (std::array<std::array<value_t, n>, h1> &a_1,
		 const std::array<value_t, n> &alpha,
		 const std::array<value_t, h1> B2_mod_b1i)
{
  for (size_t i  = 0; i < h1; i++) {
    for (size_t j = 0; j < n; j++) {
      value_t alpha_j = alpha[j];
      bool negative = false;
      if (alpha_j > (((value_t)1)<<(log_bsk-1))) {
	negative = true;
	alpha_j = -alpha_j;
	alpha_j &= mask_sk;
      }
      if (alpha_j > b1[i])
	alpha_j -= b1[i];
      value_t prod = mulmod (alpha_j, B2_mod_b1i[i], c1[i], b1[i]);
      if (negative) {
	prod = b1[i] - prod;
      }
      a_1[i][j] = submod (a_1[i][j], prod, b1[i]);
    }
  }
}

void
HPolRNS::mul (HPolRNS_ZZq &c, const HPolRNS_ZZq &a, const HPolRNS_ZZq &b)
{
  std::array<std::array<value_t, n>, h1> d_1;
  rns_pol_mul_beta (d_1, a.a_1, b.a_1, c1, b1);
  std::array<value_t, n> d_sk;
  rns_pol_mul_sk_beta (d_sk, a.a_sk, b.a_sk);
  std::array<std::array<value_t, n>, h1> xi_1;
  rns_pol_mul (xi_1, d_1, m_1, beta_m_1, c1, b1);
  std::array<std::array<value_t, n>, h2> q_2;
  rns_pol_basis_extension (q_2, xi_1, B1_over_b1j_mod_b2i, c2, b2, b1);
  std::array<value_t, n> q_sk;
  rns_pol_basis_extension_sk (q_sk, xi_1, B1_over_b1i_sk, b1);
  rns_pol_mul_beta (c.xi_2, a.xi_2, b.xi_2, c2, b2);
  rns_pol_scale (c.xi_2, c.xi_2, B2_over_b2i_B1, c2, b2);
  std::array<std::array<value_t, n>, h2> qm_2;
  rns_pol_mul (qm_2, q_2, m_2, beta_m_2, c2, b2);
  rns_pol_add (c.xi_2, c.xi_2, qm_2, b2);
  rns_pol_scale_sk (d_sk, d_sk, B1_inv_sk);
  rns_pol_mul_sk (c.a_sk, q_sk, m_sk, beta_m_sk);
  rns_pol_add_sk (c.a_sk, c.a_sk, d_sk);
  std::array<value_t, n> r_over_B2;
  rns_pol_scale_sk (r_over_B2, c.a_sk, B2_inv_sk);
  std::array<value_t, n> alpha_sk;
  rns_pol_basis_extension_sk (alpha_sk, c.xi_2, b2i_inv_sk, b2);
  rns_pol_sub_sk (alpha_sk, alpha_sk, r_over_B2);
  rns_pol_basis_extension (c.a_1, c.xi_2, B2_over_b2i_mod_b1i, c1, b1, b2);
  rns_pol_correct (c.a_1, alpha_sk, B2_mod_b1i);
}

HPolRNS_ZZq
HPolRNS::mul (const HPolRNS_ZZq &a, const HPolRNS_ZZq &b)
{
  HPolRNS_ZZq c;
  mul (c, a, b);
  return c;
}
