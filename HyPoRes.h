#include "Params.h"

struct HPolRNS_ZZq {
  std::array<std::array<value_t, n>, h1> a_1;
  std::array<std::array<value_t, n>, h2> xi_2;
  std::array<value_t, n> a_sk;
};

struct HPolRNS {
  HPolRNS ();
  HPolRNS_ZZq to_ZZq (const NTL::ZZ &a);
  NTL::ZZ to_ZZ (const HPolRNS_ZZq &A);
  void mul (HPolRNS_ZZq &c, const HPolRNS_ZZq &a, const HPolRNS_ZZq &b);
  HPolRNS_ZZq mul (const HPolRNS_ZZq &a, const HPolRNS_ZZq &b);
};
