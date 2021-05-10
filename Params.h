#include <NTL/ZZ.h>
#include <stddef.h>
#include <stdint.h>
#include <array>

#define E512

#ifdef E512
typedef uint32_t value_t;
typedef uint64_t greater_value_t;
  
constexpr size_t n = 3;
constexpr size_t h1 = 6;
constexpr size_t h2 = 5;
constexpr size_t k = 10;
constexpr size_t beta = 2;
constexpr size_t logw = 32;
  
const NTL::ZZ P (NTL::INIT_VAL, "6864797660130609714981900799081393217269435300143305409394463459185543183397656052122559640661454554977296311391480858037121987999716643812574028291115057151");
const NTL::ZZ hpr_gamma (NTL::INIT_VAL, "23945242826029513411849172299223580994042798784118784");
const std::array<NTL::ZZ, n> m = {NTL::ZZ (NTL::INIT_VAL, "-1"),
				  NTL::ZZ (NTL::INIT_VAL, "0"),
				  NTL::ZZ (NTL::INIT_VAL, "11972621413014756705924586149611790497021399392059392")};
const NTL::ZZ rho (NTL::INIT_VAL, "1197262141301475670592458614961179049702139939205939201");

const std::array<value_t, h1> b1 = {4294967197, 4294967161, 4294967029, 4294966981, 4294966927, 4294966813};
const std::array<value_t, h2> b2 = {4294967295, 4294967293, 4294967291, 4294967287, 4294967281};
constexpr size_t log_bsk = 32;

const std::array<std::array<value_t, n>, h1> rns_minv_m_1 =
  {
    std::array<value_t, n>{972066260, 1566747391, 317381676},
    std::array<value_t, n>{1165645764, 1503422105, 2787292743},
    std::array<value_t, n>{3092222663, 2127761643, 1408890592},
    std::array<value_t, n>{2362903563, 216117340, 2298552507},
    std::array<value_t, n>{715306440, 1322665808, 4275389740},
    std::array<value_t, n>{3197495866, 4124011057, 1554570859}
  };

const std::array<value_t, h2> rns_inv_B1_2 = {1290750217, 1245135352, 1311057085, 2931660051, 3639238195};
const value_t rns_inv_B1_bsk = 2834903863;

const std::array<std::array<value_t, n>, h2> rns_m_2 =
  {
    std::array<value_t, n>{4294967294, 0, 8192},
    std::array<value_t, n>{4294967292, 0, 1990656},
    std::array<value_t, n>{4294967290, 0, 25600000},
    std::array<value_t, n>{4294967286, 0, 483729408},
    std::array<value_t, n>{4294967280, 0, 1925832719}
  };

const std::array<value_t, n> rns_m_sk = {4294967295, 0, 0};

/*
  num_ops = 474
*/
#endif

#ifdef Ed448
typedef uint32_t value_t;
typedef uint64_t greater_value_t;
  
constexpr size_t n = 3;
constexpr size_t h1 = 6;
constexpr size_t h2 = 4;
constexpr size_t k = 10;
constexpr size_t beta = 2;
constexpr size_t logw = 32;
  
const NTL::ZZ P (NTL::INIT_VAL, "726838724295606890549323807888004534353641360687318060281490199180612328166730772686396383698676545930088884461843637361053498018365439");
const NTL::ZZ hpr_gamma (NTL::INIT_VAL, "617037013236874150081232979704524838882978450634686916049074439753925120009974044466434996308660572821900599056417811283480361754355140");
const std::array<NTL::ZZ, n> m = {NTL::ZZ (NTL::INIT_VAL, "-357400072521332008886097893349417451849548435"),
				  NTL::ZZ (NTL::INIT_VAL, "458411098151774579540799965020532854365487420"),
				  NTL::ZZ (NTL::INIT_VAL, "-469045874351789005827208204177461134721497131")};
const NTL::ZZ rho (NTL::INIT_VAL, "234522937175894502913604102088730567360748565501");

const std::array<value_t, h1> b1 = {4294967197, 4294967161, 4294967029, 4294966981, 4294966927, 4294966813};
const std::array<value_t, h2> b2 = {4294967295, 4294967293, 4294967291, 4294967287};
constexpr size_t log_bsk = 32;

const std::array<std::array<value_t, n>, h1> rns_minv_m_1 =
  {
    std::array<value_t, n>{2385268492, 2011451587, 3958203258},
    std::array<value_t, n>{1512448157, 4139946623, 3096119837},
    std::array<value_t, n>{2316157701, 1104920401, 1207954308},
    std::array<value_t, n>{2827138429, 3185320419, 1891021641},
    std::array<value_t, n>{371751057, 3852556880, 2887819653},
    std::array<value_t, n>{2795903457, 850984278, 1134265698}
  };

const std::array<value_t, h2> rns_inv_B1_2 = {1290750217, 1245135352, 1311057085, 2931660051};
const value_t rns_inv_B1_bsk = 2834903863;

const std::array<std::array<value_t, n>, h2> rns_m_2 =
  {
    std::array<value_t, n>{4180971720, 13650445, 388273169},
    std::array<value_t, n>{795281362, 1862036922, 2255994296},
    std::array<value_t, n>{869907867, 1784143040, 2633384481},
    std::array<value_t, n>{75804916, 467436066, 1553756707}
  };

const std::array<value_t, n> rns_m_sk = {263945581, 1476526396, 2664220629};

/*
  num_ops = 420
*/
#endif

#ifdef M383
typedef uint32_t value_t;
typedef uint64_t greater_value_t;
  
constexpr size_t n = 2;
constexpr size_t h1 = 7;
constexpr size_t h2 = 6;
constexpr size_t k = 10;
constexpr size_t beta = 3;
constexpr size_t logw = 32;
  
const NTL::ZZ P (NTL::INIT_VAL, "19701003098197239606139520050071806902539869635232723333974146702122860885748605305707133127442457820403313995153221");
const NTL::ZZ hpr_gamma (NTL::INIT_VAL, "15751658786517026077004411653439046205381357236872584912561811878014957218639809853408982590108608513434801029743689");
const std::array<NTL::ZZ, n> m = {NTL::ZZ (NTL::INIT_VAL, "-991151885317490685877537319994551485852631552512982631751"),
				  NTL::ZZ (NTL::INIT_VAL, "3668986611048655554039381508783546278051675539955460700691")};
const NTL::ZZ rho (NTL::INIT_VAL, "366898661104865555403938150878354627805167553995546070069101");

const std::array<value_t, h1> b1 = {4294967291, 4294967279, 4294967231, 4294967197, 4294967161, 4294967111, 4294967087};
const std::array<value_t, h2> b2 = {4294967295, 4294967293, 4294967287, 4294967281, 4294967273, 4294967269};
constexpr size_t log_bsk = 32;

const std::array<std::array<value_t, n>, h1> rns_minv_m_1 =
  {
    std::array<value_t, n>{829709425, 1051124881},
    std::array<value_t, n>{334487361, 3750368229},
    std::array<value_t, n>{2368309248, 2568970184},
    std::array<value_t, n>{1459942267, 3271723987},
    std::array<value_t, n>{335872509, 2381612497},
    std::array<value_t, n>{3258111871, 326558630},
    std::array<value_t, n>{3639275640, 2084938148}
  };

const std::array<value_t, h2> rns_inv_B1_2 = {4250981396, 2606369501, 2965194321, 1545768811, 1434137311, 3012873487};
const value_t rns_inv_B1_bsk = 3964892575;

const std::array<std::array<value_t, n>, h2> rns_m_2 =
  {
    std::array<value_t, n>{1215838799, 2023366281},
    std::array<value_t, n>{1513834288, 1236202703},
    std::array<value_t, n>{3756854056, 214652324},
    std::array<value_t, n>{3055935930, 648874192},
    std::array<value_t, n>{360785058, 800875203},
    std::array<value_t, n>{2111665454, 3635942859}
  };

const std::array<value_t, n> rns_m_sk = {4055513785, 3112008211};

/*
  num_ops = 340
*/
#endif
