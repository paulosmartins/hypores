#include "HyPoRes.h"
#include <vector>
#include <chrono>
#include <assert.h>

int
main (int argc, char *argv[])
{
  size_t numtests = 10000;
  NTL::SetSeed (NTL::ZZ (0));
  HPolRNS hpolrns;

  NTL::ZZ B1 (1);
  for (size_t i = 0; i < h1; i++) {
    B1 *= NTL::ZZ (b1[i]);
  }
  NTL::ZZ B1inv = NTL::InvMod (B1 % P, P);

  std::vector<HPolRNS_ZZq> xs(numtests);
  std::vector<HPolRNS_ZZq> ys(numtests);
  std::vector<HPolRNS_ZZq> zs(numtests);
  std::vector<NTL::ZZ> exp_zs(numtests);

  for (size_t i = 0; i < numtests; i++) {
    NTL::ZZ x (RandomBnd (P));
    NTL::ZZ y (RandomBnd (P));
    NTL::ZZ exp_z = (x*y*B1inv) % P;

    xs[i] = hpolrns.to_ZZq (x);
    ys[i] = hpolrns.to_ZZq (y);
    exp_zs[i] = exp_z;
  }

  //warmup
  for (size_t i = 0; i < numtests; i++) {
    hpolrns.mul (zs[i], xs[i], ys[i]);
  }

  auto start = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < numtests; i++) {
    hpolrns.mul (zs[i], xs[i], ys[i]);
  }
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;

  for (size_t i = 0; i < numtests; i++) {
    assert (hpolrns.to_ZZ (zs[i]) == exp_zs[i]);
  }

  std::cout << "average execution time = " << elapsed_seconds.count () / numtests<< "\n";
}
