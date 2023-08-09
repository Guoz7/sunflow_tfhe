//#include <stdio.h>
#include <stdlib.h>
//#include <math.h>
//#include <fftw3.h>



#include "lagrangehalfc_impl.h"
/*
#include "ni.h"
#include "ni.c"
#include "sdram.h"
*/



//#include "tfhe_io.h"

#define torusPolynomialMulR torusPolynomialMultFFT
#define torusPolynomialAddMulR torusPolynomialAddMulRFFT
#define torusPolynomialSubMulR torusPolynomialSubMulRFFT

static const int64_t _two32 = INT64_C(1) << 32; // 2^32























struct TorusPolynomial {
    int32_t N;
   Torus32* coefsT;

#ifdef __cplusplus   
   TorusPolynomial(const int32_t N);
   ~TorusPolynomial();
   TorusPolynomial(const TorusPolynomial&) = delete; //forbidden 
   TorusPolynomial* operator=(const TorusPolynomial&) = delete; //forbidden
#endif
};





struct LweParams {
	const int32_t n;
	const double alpha_min;//le plus petit bruit tq sur
	const double alpha_max;//le plus gd bruit qui permet le déchiffrement
//since all members are declared constant, a constructor is 
//required in the structure.
#ifdef __cplusplus
	LweParams(int32_t n, double alpha_min, double alpha_max);
	~LweParams();
	LweParams(const LweParams&) = delete; //forbidden
	LweParams& operator=(const LweParams& ) = delete; //forbidden
#endif
};

struct LweSample {
	Torus32* a; //-- the n coefs of the mask
    Torus32 b;  //
   	double current_variance; //-- average noise of the sample

#ifdef __cplusplus
   LweSample(const LweParams* params);
   ~LweSample();
   LweSample(const LweSample&)=delete;
   LweSample& operator=(const LweSample&)=delete;
#endif
};

struct TLweParams {
    const int32_t N; ///< a power of 2: degree of the polynomials
    const int32_t k; ///< number of polynomials in the mask
    const double alpha_min; ///< minimal noise s.t. the sample is secure
    const double alpha_max; ///< maximal noise s.t. we can decrypt
    struct LweParams extracted_lweparams; ///< lwe params if one extracts

#ifdef __cplusplus

    TLweParams(int32_t N, int32_t k, double alpha_min, double alpha_max);

    ~TLweParams();

    TLweParams(const TLweParams &) = delete;

    void operator=(const TLweParams &) = delete;

#endif
};

struct TLweSample {
    struct TorusPolynomial *a; ///< array of length k+1: mask + right term
    struct TorusPolynomial *b; ///< alias of a[k] to get the right term
    double current_variance; ///< avg variance of the sample
    int32_t k;
#ifdef __cplusplus

    TLweSample(const TLweParams *params);

    ~TLweSample();

    TLweSample(const TLweSample &) = delete;

    void operator=(const TLweSample &) = delete;

#endif
};

struct IntPolynomial {
    int32_t N;
   int32_t* coefs;

#ifdef __cplusplus   
   IntPolynomial(const int32_t N);
   ~IntPolynomial();
   IntPolynomial(const IntPolynomial&) = delete; //forbidden 
   IntPolynomial* operator=(const IntPolynomial&) = delete; //forbidden
#endif
};



struct TGswParams {
     int32_t l; ///< decomp length
     int32_t Bgbit;///< log_2(Bg)
     int32_t Bg;///< decomposition base (must be a power of 2)
     int32_t halfBg; ///< Bg/2
     uint32_t maskMod; ///< Bg-1
     struct TLweParams *tlwe_params; ///< Params of each row
     int32_t kpl; ///< number of rows = (k+1)*l
    Torus32 *h; ///< powers of Bgbit
    uint32_t offset; ///< offset = Bg/2 * (2^(32-Bgbit) + 2^(32-2*Bgbit) + ... + 2^(32-l*Bgbit))

#ifdef __cplusplus

    TGswParams(int32_t l, int32_t Bgbit, const TLweParams *tlwe_params);

    ~TGswParams();

    TGswParams(const TGswParams &) = delete;

    void operator=(const TGswParams &) = delete;

#endif
};

struct TGswSample {
    struct TLweSample *all_sample; ///< TLweSample* all_sample; (k+1)l TLwe Sample
    struct TLweSample **bloc_sample;///< accès optionnel aux différents blocs de taille l.
    // double current_variance;
     int32_t k;
     int32_t l;

#ifdef __cplusplus

    inline TGswSample(TLweSample *all_sample, TLweSample **bloc_sample, const int32_t k, const int32_t l) :
            all_sample(all_sample),
            bloc_sample(bloc_sample),
            k(k), l(l) {}

    inline ~TGswSample() {}

    TGswSample(const TGswSample &) = delete;

    void operator=(const TGswSample &) = delete;

#endif
};


struct LweKeySwitchKey {
    int32_t n; ///< length of the input key: s'
    int32_t t; ///< decomposition length
    int32_t basebit; ///< log_2(base)
    int32_t base; ///< decomposition base: a power of 2 
    const struct  LweParams* out_params; ///< params of the output key s 
    struct LweSample* ks0_raw; //tableau qui contient tout les Lwe samples de taille nlbase
    struct LweSample** ks1_raw;// de taille nl  pointe vers un tableau ks0_raw dont les cases sont espaceés de base positions
    struct LweSample*** ks; ///< the keyswitch elements: a n.l.base matrix
    // de taille n pointe vers ks1 un tableau dont les cases sont espaceés de ell positions

#ifdef __cplusplus
    LweKeySwitchKey(int32_t n, int32_t t, int32_t basebit, const LweParams* out_params, LweSample* ks0_raw);
    ~LweKeySwitchKey();
    LweKeySwitchKey(const LweKeySwitchKey&) = delete;
    void operator=(const LweKeySwitchKey&) = delete;
#endif
};

struct LweBootstrappingKey{
    const struct LweParams* in_out_params; ///< paramètre de l'input et de l'output. key: s
    const  struct TGswParams* bk_params; ///< params of the Gsw elems in bk. key: s"
    const struct TLweParams* accum_params; ///< params of the accum variable key: s"
    const struct LweParams* extract_params; ///< params after extraction: key: s' 
    struct TGswSample* bk; ///< the bootstrapping key (s->s")
    struct LweKeySwitchKey* ks; ///< the keyswitch key (s'->s)


#ifdef __cplusplus
   LweBootstrappingKey(const LweParams* in_out_params, 
    const TGswParams* bk_params,
    const TLweParams* accum_params,
    const LweParams* extract_params,
    TGswSample* bk,
    LweKeySwitchKey* ks);
    ~LweBootstrappingKey();
    LweBootstrappingKey(const LweBootstrappingKey&) = delete;
    void operator=(const LweBootstrappingKey&) = delete;
  
#endif
};


struct TFheGateBootstrappingParameterSet {
    const int32_t ks_t;
    const int32_t ks_basebit;
    const struct LweParams *const in_out_params;
    const struct TGswParams *const tgsw_params;
#ifdef __cplusplus

    TFheGateBootstrappingParameterSet(const int32_t ks_t, const int32_t ks_basebit, const LweParams *const in_out_params,
                                      const TGswParams *const tgsw_params);

    TFheGateBootstrappingParameterSet(const TFheGateBootstrappingParameterSet &) = delete;

    void operator=(const TFheGateBootstrappingParameterSet &)= delete;

#endif
};

struct LweKey {
   const struct LweParams* params;
   int32_t* key;

#ifdef __cplusplus   
   LweKey(const LweParams* params);
   ~LweKey();
   LweKey(const LweKey&) = delete; //forbidden 
   LweKey* operator=(const LweKey&) = delete; //forbidden
#endif
};
struct TLweKey {
    const struct TLweParams *params; ///< the parameters of the key
    struct IntPolynomial *key; ///< the key (i.e k binary polynomials)
#ifdef __cplusplus

    TLweKey(const TLweParams *params);

    ~TLweKey();

    TLweKey(const TLweKey &) = delete;

    void operator=(const TLweKey &) = delete;

#endif
};

struct TGswKey {
    const struct TGswParams *params; ///< the parameters
    const struct TLweParams *tlwe_params; ///< the tlwe params of each rows
    struct IntPolynomial *key; ///< the key (array of k polynomials)
    struct TLweKey tlwe_key;

#ifdef __cplusplus

    TGswKey(const TGswParams *params);

    ~TGswKey();

    TGswKey(const TGswKey &) = delete;

    void operator=(const TGswKey &) = delete;

#endif
};
struct TLweSampleFFT {
    struct LagrangeHalfCPolynomial *a; ///< array of length k+1: mask + right term
    struct LagrangeHalfCPolynomial *b; ///< alias of a[k] to get the right term
    double current_variance; ///< avg variance of the sample
     int32_t k; //required during the destructor call...
#ifdef __cplusplus

    TLweSampleFFT(const TLweParams *params, LagrangeHalfCPolynomial *a, double current_variance);

    ~TLweSampleFFT();

    TLweSampleFFT(const TLweSampleFFT &) = delete;

    void operator=(const TLweSampleFFT &) = delete;

#endif
};
struct TGswSampleFFT {
    struct TLweSampleFFT *all_samples; ///< TLweSample* all_sample; (k+1)l TLwe Sample
    struct TLweSampleFFT **sample; ///< accès optionnel aux différents blocs de taille l.
    //double current_variance;
    const int32_t k;
    const int32_t l;

#ifdef __cplusplus

    TGswSampleFFT(const TGswParams *params, TLweSampleFFT *all_samples);

    ~TGswSampleFFT();

    TGswSampleFFT(const TGswSampleFFT &) = delete;

    void operator=(const TGswSampleFFT &) = delete;

#endif
};
struct LweBootstrappingKeyFFT {
    const struct LweParams* in_out_params; ///< paramètre de l'input et de l'output. key: s
    const struct TGswParams* bk_params; ///< params of the Gsw elems in bk. key: s"
    const struct TLweParams* accum_params; ///< params of the accum variable key: s"
    const struct LweParams* extract_params; ///< params after extraction: key: s' 
    const struct TGswSampleFFT* bkFFT; ///< the bootstrapping key (s->s")
    const struct LweKeySwitchKey* ks; ///< the keyswitch key (s'->s)


#ifdef __cplusplus
   LweBootstrappingKeyFFT(const LweParams* in_out_params, 
    const TGswParams* bk_params,
    const TLweParams* accum_params,
    const LweParams* extract_params, 
    const TGswSampleFFT* bkFFT,
    const LweKeySwitchKey* ks);
    ~LweBootstrappingKeyFFT();
    LweBootstrappingKeyFFT(const LweBootstrappingKeyFFT&) = delete;
    void operator=(const LweBootstrappingKeyFFT&) = delete;
  
#endif


};

struct TFheGateBootstrappingCloudKeySet {
    const struct TFheGateBootstrappingParameterSet *const params;
    const struct LweBootstrappingKey *const bk;
    //const LweBootstrappingKeyFFT *const bkFFT;
#ifdef __cplusplus

    TFheGateBootstrappingCloudKeySet(
            const TFheGateBootstrappingParameterSet *const params,
            const LweBootstrappingKey *const bk,
            const LweBootstrappingKeyFFT *const bkFFT);

    TFheGateBootstrappingCloudKeySet(const TFheGateBootstrappingCloudKeySet &) = delete;

    void operator=(const TFheGateBootstrappingCloudKeySet &)= delete;

#endif
};

struct TFheGateBootstrappingSecretKeySet {
    const struct TFheGateBootstrappingParameterSet *params;
    const struct LweKey *lwe_key;
    const struct TGswKey *tgsw_key;
    const struct TFheGateBootstrappingCloudKeySet cloud;
#ifdef __cplusplus

    TFheGateBootstrappingSecretKeySet(
            const TFheGateBootstrappingParameterSet *const params,
            const LweBootstrappingKey *const bk,
            const LweBootstrappingKeyFFT *const bkFFT,
            const LweKey *lwe_key,
            const TGswKey *tgsw_key);

    TFheGateBootstrappingSecretKeySet(const TFheGateBootstrappingSecretKeySet &) = delete;

    void operator=(const TFheGateBootstrappingSecretKeySet &)= delete;

#endif
};
/** This structure represents an torus polynomial modulo X^N+1 */


/** 
 * This structure is used for FFT operations, and is a representation
 * over C of a polynomial in R[X]/X^N+1
 * This type is meant to be specialized, and all implementations of the structure must be compatible 
 * (reinterpret_cast) with this one. Namely, they should contain at most 2 pointers 
 */
struct LagrangeHalfCPolynomial
{
   void* data;
   void* precomp;
};



// Idea:
// we may want to represent an element x of the real torus by 
// the integer rint(2^32.x) modulo 2^32
//  -- addition, subtraction and integer combinations are native operation
//  -- modulo 1 is mapped to mod 2^32, which is also native!
// This looks much better than using float/doubles, where modulo 1 is not
// natural at all.

//typedef int64_t Torus64; //avant uint64_t
struct LweParams;
struct LweKey;
struct LweSample;
struct LweKeySwitchKey;
struct TLweParams;
struct TLweKey;
struct TLweSample;
struct TLweSampleFFT;
struct TGswParams;
struct TGswKey;
struct TGswSample;
struct TGswSampleFFT;
struct LweBootstrappingKey;
struct LweBootstrappingKeyFFT;
struct IntPolynomial;
struct TorusPolynomial;
struct LagrangeHalfCPolynomial;
struct TFheGateBootstrappingParameterSet;
struct TFheGateBootstrappingCloudKeySet;
struct TFheGateBootstrappingSecretKeySet;

//this is for compatibility with C code, to be able to use
//"LweParams" as a type and not "struct LweParams"
typedef struct LweParams           LweParams;
typedef struct LweKey              LweKey;
typedef struct LweSample           LweSample;
typedef struct LweKeySwitchKey     LweKeySwitchKey;
typedef struct TLweParams       TLweParams;
typedef struct TLweKey          TLweKey;
typedef struct TLweSample       TLweSample;
typedef struct TLweSampleFFT       TLweSampleFFT;
typedef struct TGswParams       TGswParams;
typedef struct TGswKey          TGswKey;
typedef struct TGswSample       TGswSample;
typedef struct TGswSampleFFT       TGswSampleFFT;
typedef struct LweBootstrappingKey LweBootstrappingKey;
typedef struct LweBootstrappingKeyFFT LweBootstrappingKeyFFT;
typedef struct IntPolynomial	   IntPolynomial;
typedef struct TorusPolynomial	   TorusPolynomial;
typedef struct LagrangeHalfCPolynomial	   LagrangeHalfCPolynomial;
typedef struct TFheGateBootstrappingParameterSet TFheGateBootstrappingParameterSet;
typedef struct TFheGateBootstrappingCloudKeySet TFheGateBootstrappingCloudKeySet;
typedef struct TFheGateBootstrappingSecretKeySet TFheGateBootstrappingSecretKeySet;

// TorusPolynomial的初始化函数
void TorusPolynomial_init(TorusPolynomial* poly,  int32_t N) {
    poly->N = N;
    poly->coefsT = (Torus32*)malloc(N * sizeof(Torus32));
}

// TorusPolynomial的清理函数（类似于析构函数）
void TorusPolynomial_destroy(struct TorusPolynomial* poly) {
    free(poly->coefsT);
}








//#include "tfhe.h"
#include "tfhe_io.h"


//  wle functionfnew_TorusPolynomial_array

EXPORT void lweNoiselessTrivial(LweSample* result, Torus32 mu, const LweParams* params){
    const int32_t n = params->n;

    for (int32_t i = 0; i < n; ++i) result->a[i] = 0;
    result->b = mu;
    result->current_variance = 0.;
}


/** result = result - sample */
EXPORT void lweSubTo(LweSample* result, const LweSample* sample, const LweParams* params){
    const int32_t n = params->n;
    const Torus32* __restrict sa = sample->a;
    Torus32* __restrict ra = result->a;

#ifdef __AVX2__
    intVecSubTo_avx(ra,sa,n);
#else
    for (int32_t i = 0; i < n; ++i) ra[i] -= sa[i];
#endif
    result->b -= sample->b;
    result->current_variance += sample->current_variance; 
}



EXPORT void lweAddTo(LweSample* result, const LweSample* sample, const LweParams* params){
    const int32_t n = params->n;

    for (int32_t i = 0; i < n; ++i) result->a[i] += sample->a[i];
    result->b += sample->b;
    result->current_variance += sample->current_variance; 
}


//typedef unsigned long uint64_t1;



EXPORT int32_t modSwitchFromTorus32(Torus32 phase1, int32_t Msize){
    uint64_t interv = ((UINT64_C(1)<<63)/Msize)*2; // width of each intervall
    uint64_t half_interval = interv/2; // begin of the first intervall
    uint64_t phase64 = (uint64_t)phase1;
    phase64 = (phase64 << 32) + half_interval;
   // << 32) + half_interval;
    //floor to the nearest multiples of interv
    return phase64/interv;
}


EXPORT TorusPolynomial* new_TorusPolynomial(const int32_t N) {
    TorusPolynomial* tmp;
    TorusPolynomial_init(tmp,N);
    return tmp;
}

// TorusPolynomial::TorusPolynomial(const int32_t N): N(N)
// {
//     this->coefsT = new Torus32[N]; 
// }

EXPORT TorusPolynomial* alloc_TorusPolynomial_array(int32_t nbelts) {
    return (TorusPolynomial*) malloc(nbelts*sizeof(TorusPolynomial));
}

void init_TorusPolynomial_array(int32_t nbelts, struct TorusPolynomial* obj, const int32_t N) {
    for (int32_t i = 0; i < nbelts; i++) {
        TorusPolynomial_init(&obj[i], N);
    }
}


EXPORT TorusPolynomial* new_TorusPolynomial_array(int32_t nbelts, const int32_t N) {
    TorusPolynomial* obj = alloc_TorusPolynomial_array(nbelts);
    init_TorusPolynomial_array(nbelts,obj,N);
    return obj;
}


TLweSample *new_TLweSample(const TLweParams* params) {
    const int32_t k = params->k;
    TLweSample* result = (TLweSample*) malloc(sizeof(TLweSample));
    result->a = new_TorusPolynomial_array(k+1, params->N);
    result->b = result->a + k;
    result->current_variance = 0.;
    return result;
}
struct TLweSample *new_TLweSample_array(int32_t nbelts, const struct TLweParams *params) {
    struct TLweSample *obj = (struct TLweSample *)malloc(nbelts * sizeof(struct TLweSample));
    for (int32_t i = 0; i < nbelts; i++) {
        obj[i] = *new_TLweSample(params);
    }
    return obj;
}

// TorusPolynomial = 0
EXPORT void torusPolynomialClear(TorusPolynomial *result) {
    const int32_t N = result->N;

    for (int32_t i = 0; i < N; ++i) result->coefsT[i] = 0;
}

EXPORT void tfhe_blindRotate(TLweSample* accum, const TGswSample* bk, const int32_t* bara, const int32_t n, const TGswParams* bk_params);
EXPORT void tfhe_blindRotateAndExtract(LweSample* result, const TorusPolynomial* v, const TGswSample* bk, const int32_t barb, const int32_t* bara, const int32_t n, const TGswParams* bk_params);
EXPORT void tfhe_bootstrap_woKS(LweSample* result, const LweBootstrappingKey* bk, Torus32 mu, const LweSample* x);
EXPORT void tfhe_bootstrap(LweSample* result, const LweBootstrappingKey* bk, Torus32 mu, const LweSample* x);
EXPORT void tfhe_createLweBootstrappingKey(LweBootstrappingKey* bk, const LweKey* key_in, const TGswKey* rgsw_key);




EXPORT void die_dramatically(const char* message);

// // //torusPolynomialAddMulRFFT
// class FFT_Processor_fftw {
//     public:
//     const int32_t _2N;
//     const int32_t N;    
//     const int32_t Ns2;
//     private:
//     double* rev_in;
//     fftw_complex* rev_out;
//     fftw_complex* in;
//     double* out;
//     fftw_plan p;
//     fftw_plan rev_p;
//     void plan_fftw();
//     public:
//     cplx* omegaxminus1;

//     FFT_Processor_fftw(const int32_t N);
//     void execute_reverse_int(cplx* res, const int32_t* a);
//     void execute_reverse_torus32(cplx* res, const Torus32* a);
//     void execute_direct_Torus32(Torus32* res, const cplx* a);
//     ~FFT_Processor_fftw();
// };

// struct cplx {
//     double real; // 实部
//     double imag; // 虚部
// };

// // 定义FFT_Processor_fftw结构体
// struct FFT_Processor_fftw {
//     const int32_t _2N;
//     const int32_t N;    
//     const int32_t Ns2;
//     double* rev_in;
//     struct cplx* rev_out;
//     struct cplx* in;
//     double* out;
//     void (*plan_fftw)();
//     struct cplx* omegaxminus1;
// };

// thread_local struct FFT_Processor_fftw fp1024_fftw;

/////  use the new method to represent the fushu

//typedef std::complex<double> cplx;
// struct LagrangeHalfCPolynomial_IMPL
// {
//    double* coefsC;
//    struct FFT_Processor_fftw* proc;

//    LagrangeHalfCPolynomial_IMPL(int32_t N);
//    ~LagrangeHalfCPolynomial_IMPL();
// };



// extern struct FFT_Processor_fftw ;
// struct FFt_Processor_fftw fp1024_fftw;



// EXPORT void IntPolynomial_ifft(LagrangeHalfCPolynomial* result, const IntPolynomial* p) {
//     fp1024_fftw.execute_reverse_int = execute_reverse_int;
//     fp1024_fftw.execute_reverse_int(((struct LagrangeHalfCPolynomial_IMPL*)result)->coefsC, p->coefs);
// }
// EXPORT void TorusPolynomial_ifft(LagrangeHalfCPolynomial* result, const TorusPolynomial* p) {
//     fp1024_fftw.execute_reverse_torus32 = execute_reverse_torus32;  
//     fp1024_fftw.execute_reverse_torus32(((struct LagrangeHalfCPolynomial_IMPL*)result)->coefsC, p->coefsT);
// }
// EXPORT void TorusPolynomial_fft(TorusPolynomial* result, const LagrangeHalfCPolynomial* p) {
//     fp1024_fftw.execute_direct_Torus32 = execute_direct_Torus32;
//     fp1024_fftw.execute_direct_Torus32(result->coefsT, ((struct LagrangeHalfCPolynomial_IMPL*)p)->coefsC);
// }


EXPORT LagrangeHalfCPolynomial* alloc_LagrangeHalfCPolynomial_array(int32_t nbelts) {
    return (LagrangeHalfCPolynomial*) malloc(nbelts*sizeof(LagrangeHalfCPolynomial));
}



// void LagrangeHalfCPolynomial_IMPL_init(struct LagrangeHalfCPolynomial_IMPL* obj, const int32_t N) {
//     ////assert(N == 1024);
//     //obj->N = N;

//     // 分配内存并初始化 coefsC 数组
//     obj->coefsC = (struct cplx*)malloc((N / 2) * sizeof(struct cplx));
//     // 这里可以根据需要对 coefsC 数组进行初始化
//     // ...
//     // 假设 fp1024_fftw 是一个 FFT_Processor_fftw 结构体的全局变量，
//     // 这里将 obj 的 proc 成员指向 fp1024_fftw 结构体的地址
//     // 如果需要动态分配 proc 成员，也可以在这里进行分配
//     obj->proc = &fp1024_fftw;
// }

// void init_LagrangeHalfCPolynomial(struct LagrangeHalfCPolynomial* obj,const int32_t N){
//     obj.coefsC = struct cplx[N/2];
//     obj.proc = &fp1024_fftw;
// }


EXPORT void init_LagrangeHalfCPolynomial_array(int32_t nbelts, LagrangeHalfCPolynomial* obj, const int32_t N) {
    for (int32_t i=0; i<nbelts; i++) {
	//  LagrangeHalfCPolynomial_IMPL_init(&(obj[i]),N);
    }
    }
EXPORT LagrangeHalfCPolynomial* new_LagrangeHalfCPolynomial_array(int32_t nbelts, const int32_t N) {
    LagrangeHalfCPolynomial* obj = alloc_LagrangeHalfCPolynomial_array(nbelts);
    init_LagrangeHalfCPolynomial_array(nbelts,obj,N);
    return obj;
}

// TorusPolynomial += TorusPolynomial
EXPORT void torusPolynomialAddTo(TorusPolynomial *result, const TorusPolynomial *poly2) {
    const int32_t N = poly2->N;
    Torus32 *r = result->coefsT;
    const Torus32 *b = poly2->coefsT;

    for (int32_t i = 0; i < N; ++i)
        r[i] += b[i];
}


EXPORT void torusPolynomialAddMulRFFT(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2) {
    return;
    const int32_t N = poly1->N;
    struct LagrangeHalfCPolynomial* tmp = new_LagrangeHalfCPolynomial_array(3,N);
    struct TorusPolynomial* tmpr = new_TorusPolynomial(N);
    IntPolynomial_ifft(tmp+0,poly1);
    TorusPolynomial_ifft(tmp+1,poly2);
    LagrangeHalfCPolynomialMul(tmp+2,tmp+0,tmp+1);
    TorusPolynomial_fft(tmpr, tmp+2);
    torusPolynomialAddTo(result, tmpr);
    delete_TorusPolynomial(tmpr);
    delete_LagrangeHalfCPolynomial_array(3,tmp);
}
//tGswExternMulToTLwe

EXPORT IntPolynomial* alloc_IntPolynomial_array(int32_t nbelts) {
    return (IntPolynomial*) malloc(nbelts*sizeof(IntPolynomial));
}

// IntPolynomial::IntPolynomial(const int32_t N): N(N)
// {
//     this->coefs = new int32_t[N]; 
// }

// EXPORT void init_IntPolynomial_array(int32_t nbelts, IntPolynomial* obj, const int32_t N) {
//     for (int32_t i=0; i<nbelts; i++) {
// 	new(obj+i) IntPolynomial(N);
//     }
// }

void init_IntPolynomial_array(int32_t nbelts, struct IntPolynomial* obj, const int32_t N) {
    for (int32_t i = 0; i < nbelts; i++) {
        // 为每个元素分配内存
        obj[i].coefs = (int32_t*)malloc(N * sizeof(int32_t));

        // 初始化 N
        obj[i].N = N;

        // 在这里可以根据需要对 coefs 数组进行初始化
        // 例如，可以将 coefs 数组的元素全部初始化为0
        for (int32_t j = 0; j < N; j++) {
            obj[i].coefs[j] = 0;
        }
    }
}

EXPORT IntPolynomial* new_IntPolynomial_array(int32_t nbelts, const int32_t N) {
    IntPolynomial* obj = alloc_IntPolynomial_array(nbelts);
    init_IntPolynomial_array(nbelts,obj,N);
    return obj;
}



#define INCLUDE_TGSW_TORUS32POLYNOMIAL_DECOMP_H 


#if defined INCLUDE_ALL || defined INCLUDE_TGSW_TORUS32POLYNOMIAL_DECOMP_H
#undef INCLUDE_TGSW_TORUS32POLYNOMIAL_DECOMP_H
EXPORT void
tGswTorus32PolynomialDecompH(IntPolynomial *result, const TorusPolynomial *sample, const TGswParams *params) {
    const int32_t N = params->tlwe_params->N;
    const int32_t l = params->l;
    const int32_t Bgbit = params->Bgbit;
    uint32_t *buf = (uint32_t *) sample->coefsT;
//#define __AVX2__ //(to test)
#ifndef __AVX2__
    const uint32_t maskMod = params->maskMod;
    const int32_t halfBg = params->halfBg;
    const uint32_t offset = params->offset;
#else
    const uint32_t* maskMod_addr = &params->maskMod;
    const int32_t* halfBg_addr = &params->halfBg;
    const uint32_t* offset_addr = &params->offset;
    //const uint32_t offset = params->offset;
    //const uint32_t maskMod = params->maskMod;
    //const int32_t halfBg = params->halfBg;
#endif

    //First, add offset to everyone
#ifndef __AVX2__
    for (int32_t j = 0; j < N; ++j) buf[j] += offset;
#else
    {
    const uint32_t* sit = buf;
    const uint32_t* send = buf+N;
    __asm__ __volatile__ (
        "vpbroadcastd (%2),%%ymm0\n"
        "0:\n"
        "vmovdqu (%0),%%ymm3\n"
        "vpaddd %%ymm0,%%ymm3,%%ymm3\n" // add offset
        "vmovdqu %%ymm3,(%0)\n"
        "addq $32,%0\n"
        "cmpq %1,%0\n"
        "jb 0b\n"
        : "=r"(sit),"=r"(send),"=r"(offset_addr)
        :  "0"(sit), "1"(send), "2"(offset_addr)
        : "%ymm0","%ymm3","memory"
        );
    }
#endif

    //then, do the decomposition (in parallel)
    for (int32_t p = 0; p < l; ++p) {
        const int32_t decal = (32 - (p + 1) * Bgbit);
#ifndef __AVX2__
        int32_t *res_p = result[p].coefs;
        for (int32_t j = 0; j < N; ++j) {
            uint32_t temp1 = (buf[j] >> decal) & maskMod;
            res_p[j] = temp1 - halfBg;
        }
#else
        int32_t* dst = result[p].coefs;
        const uint32_t* sit = buf;
        const uint32_t* send = buf+N;
        const int32_t* decal_addr = &decal;
        __asm__ __volatile__ (
            "vpbroadcastd (%4),%%ymm0\n"
            "vpbroadcastd (%5),%%ymm1\n"
            "vmovd (%3),%%xmm2\n"
            "1:\n"
            "vmovdqu (%1),%%ymm3\n"
            "VPSRLD %%xmm2,%%ymm3,%%ymm3\n" // shift by decal
            "VPAND %%ymm1,%%ymm3,%%ymm3\n"  // and maskMod
            "VPSUBD %%ymm0,%%ymm3,%%ymm3\n" // sub halfBg
            "vmovdqu %%ymm3,(%0)\n"
            "addq $32,%0\n"
            "addq $32,%1\n"
            "cmpq %2,%1\n"
            "jb 1b\n"
            : "=r"(dst),"=r"(sit),"=r"(send),"=r"(decal_addr),"=r"(halfBg_addr),"=r"(maskMod_addr)
            :  "0"(dst), "1"(sit), "2"(send), "3"(decal_addr), "4"(halfBg_addr) ,"5"(maskMod_addr)
            : "%ymm0","%ymm1","%ymm2","%ymm3","memory"
            );
        /* // verify that the assembly block was ok
        int32_t* res_p = result[p].coefs;
        for (int32_t j = 0; j < N; ++j)
        {
            uint32_t temp1 = (buf[j] >> decal) & maskMod;
            if (res_p[j] != int32_t(temp1 - halfBg)) {
            fprintf(stderr, "j=%d,buf[j]=%u,decal=%u,mask=%u,halfbg=%d,res_p[j]=%d\n",j,buf[j],decal,maskMod,halfBg,res_p[j]);
            abort();
            }
        }*/

#endif
    }

    //finally, remove offset to everyone
#ifndef __AVX2__
    for (int32_t j = 0; j < N; ++j) buf[j] -= offset;
#else
    {
    const uint32_t* sit = buf;
    const uint32_t* send = buf+N;
    __asm__ __volatile__ (
        "vpbroadcastd (%2),%%ymm0\n"
        "2:\n"
        "vmovdqu (%0),%%ymm3\n"
        "vpsubd %%ymm0,%%ymm3,%%ymm3\n" // add offset
        "vmovdqu %%ymm3,(%0)\n"
        "addq $32,%0\n"
        "cmpq %1,%0\n"
        "jb 2b\n"
        "vzeroall\n"
        : "=r"(sit),"=r"(send),"=r"(offset_addr)
        :  "0"(sit), "1"(send), "2"(offset_addr)
        : "%ymm0","%ymm3","memory"
        );
    }
#endif
}
#endif

EXPORT void tGswTLweDecompH(IntPolynomial *result, const TLweSample *sample, const TGswParams *params) {
    const int32_t k = params->tlwe_params->k;
    const int32_t l = params->l;

    for (int32_t i = 0; i <= k; ++i) // b=a[k]
        tGswTorus32PolynomialDecompH(result + (i * l), &sample->a[i], params);
}

//Arithmetic operations on TLwe samples
/** result = (0,0) */
EXPORT void tLweClear(TLweSample *result, const TLweParams *params) {
    const int32_t k = params->k;

    for (int32_t i = 0; i < k; ++i) torusPolynomialClear(&result->a[i]);
    torusPolynomialClear(result->b);
    result->current_variance = 0.;
}




// Norme Euclidienne d'un IntPolynomial
EXPORT double intPolynomialNormSq2(const IntPolynomial *poly) {
    const int32_t N = poly->N;
    int32_t temp1 = 0;

    for (int32_t i = 0; i < N; ++i) {
        int32_t temp0 = poly->coefs[i] * poly->coefs[i];
        temp1 += temp0;
    }
    return temp1;
}


EXPORT void
tLweAddMulRTo(TLweSample *result, const IntPolynomial *p, const TLweSample *sample, const TLweParams *params) {
    const int32_t k = params->k;

    for (int32_t i = 0; i <= k; ++i)
        torusPolynomialAddMulR(result->a + i, p, sample->a + i);
    result->current_variance += intPolynomialNormSq2(p) * sample->current_variance;
}

EXPORT void tGswExternMulToTLwe(TLweSample *accum, const TGswSample *sample, const TGswParams *params) {
    const TLweParams *par = params->tlwe_params;
    const int32_t N = par->N;
    const int32_t kpl = params->kpl;
    //TODO: improve this new/delete
    IntPolynomial *dec = new_IntPolynomial_array(kpl, N);

    tGswTLweDecompH(dec, accum, params);
    tLweClear(accum, par);
    for (int32_t i = 0; i < kpl; i++) {
        tLweAddMulRTo(accum, &dec[i], &sample->all_sample[i], par);
    }

    //delete_IntPolynomial_array(kpl, dec);
}



//tfhe_MuxRotate

//result= (X^{a}-1)*source
EXPORT void torusPolynomialMulByXaiMinusOne(TorusPolynomial *result, int32_t a, const TorusPolynomial *source) {
    const int32_t N = source->N;
    Torus32 *out = result->coefsT;
    Torus32 *in = source->coefsT;

    //assert(a >= 0 && a < 2 * N);

    if (a < N) {
        for (int32_t i = 0; i < a; i++)//sur que i-a<0
            out[i] = -in[i - a + N] - in[i];
        for (int32_t i = a; i < N; i++)//sur que N>i-a>=0
            out[i] = in[i - a] - in[i];
    } else {
        const int32_t aa = a - N;
        for (int32_t i = 0; i < aa; i++)//sur que i-a<0
            out[i] = in[i - aa + N] - in[i];
        for (int32_t i = aa; i < N; i++)//sur que N>i-a>=0
            out[i] = -in[i - aa] - in[i];
    }
}

//mult externe de X^ai-1 par bki
EXPORT void tLweMulByXaiMinusOne(TLweSample *result, int32_t ai, const TLweSample *bk, const TLweParams *params) {
    const int32_t k = params->k;
    for (int32_t i = 0; i <= k; i++)
        torusPolynomialMulByXaiMinusOne(&result->a[i], ai, &bk->a[i]);
}

EXPORT void tLweAddTTo(TLweSample *result, const int32_t pos, const Torus32 x, const TLweParams *params) {
    result->a[pos].coefsT[0] += x;
}



EXPORT void tLweAddTo(TLweSample *result, const TLweSample *sample, const TLweParams *params) {
    const int32_t k = params->k;

    for (int32_t i = 0; i < k; ++i)
        torusPolynomialAddTo(&result->a[i], &sample->a[i]);
    torusPolynomialAddTo(result->b, sample->b);
    result->current_variance += sample->current_variance;
}


void tfhe_MuxRotate(TLweSample *result, const TLweSample *accum, const TGswSample *bki, const int32_t barai,
                    const TGswParams *bk_params) {
    // ACC = BKi*[(X^barai-1)*ACC]+ACC
    // temp = (X^barai-1)*ACC
    tLweMulByXaiMinusOne(result, barai, accum, bk_params->tlwe_params);
    // temp *= BKi
    tGswExternMulToTLwe(result, bki, bk_params);
    // ACC += temp
    tLweAddTo(result, accum, bk_params->tlwe_params);
}





//tfhe_blindRotate
/** result = sample */
EXPORT void tLweCopy(TLweSample *result, const TLweSample *sample, const TLweParams *params) {
    const int32_t k = params->k;
    const int32_t N = params->N;

    for (int32_t i = 0; i <= k; ++i)
        for (int32_t j = 0; j < N; ++j)
            result->a[i].coefsT[j] = sample->a[i].coefsT[j];

    result->current_variance = sample->current_variance;
}



EXPORT void
tfhe_blindRotate(TLweSample *accum, const TGswSample *bk, const int32_t *bara, const int32_t n, const TGswParams *bk_params) {

    //TGswSample* temp = new_TGswSample(bk_params);
    TLweSample *temp = new_TLweSample(bk_params->tlwe_params);
    TLweSample *temp2 = temp;
    TLweSample *temp3 = accum;

    for (int32_t i = 0; i < n; i++) {
        const int32_t barai = bara[i];
        if (barai == 0) continue; //indeed, this is an easy case!

        tfhe_MuxRotate(temp2, temp3, bk + i, barai, bk_params);
        //swap(temp2, temp3);
        temp = temp2;
        temp2 = temp3;
        temp3 = temp;

    }
    if (temp3 != accum) {
        tLweCopy(accum, temp3, bk_params->tlwe_params);
    }

    //delete_TLweSample(temp);
    //delete_TGswSample(temp);
}

//tfhe_blindRotateAndExtract





// TorusPolynomial = TorusPolynomial


//result= X^{a}*source
EXPORT void torusPolynomialMulByXai(TorusPolynomial *result, int32_t a, const TorusPolynomial *source) {
    const int32_t N = source->N;
    Torus32 *out = result->coefsT;
    Torus32 *in = source->coefsT;

    //assert(a >= 0 && a < 2 * N);
    //assert(result != source);

    if (a < N) {
        for (int32_t i = 0; i < a; i++)//sur que i-a<0
            out[i] = -in[i - a + N];
        for (int32_t i = a; i < N; i++)//sur que N>i-a>=0
            out[i] = in[i - a];
    } else {
        const int32_t aa = a - N;
        for (int32_t i = 0; i < aa; i++)//sur que i-a<0
            out[i] = in[i - aa + N];
        for (int32_t i = aa; i < N; i++)//sur que N>i-a>=0
            out[i] = -in[i - aa];
    }
}


EXPORT void torusPolynomialCopy(
        TorusPolynomial *result,
        const TorusPolynomial *sample) {
    //assert(result != sample);
    const int32_t N = result->N;
    const Torus32 *__restrict s = sample->coefsT;
    Torus32 *__restrict r = result->coefsT;

    for (int32_t i = 0; i < N; ++i) r[i] = s[i];
}

/** result = (0,mu) */
EXPORT void tLweNoiselessTrivial(TLweSample *result, const TorusPolynomial *mu, const TLweParams *params) {
    const int32_t k = params->k;

    for (int32_t i = 0; i < k; ++i) torusPolynomialClear(&result->a[i]);
    torusPolynomialCopy(result->b, mu);
    result->current_variance = 0.;
}


EXPORT void tLweExtractLweSampleIndex(LweSample* result, const TLweSample* x, const int32_t index, const LweParams* params,  const TLweParams* rparams) {
    const int32_t N = rparams->N;
    const int32_t k = rparams->k;
    //assert(params->n == k*N);

    for (int32_t i=0; i<k; i++) {
      for (int32_t j=0; j<=index; j++)
        result->a[i*N+j] = x->a[i].coefsT[index-j];
      for (int32_t j=index+1; j<N; j++)
        result->a[i*N+j] = -x->a[i].coefsT[N+index-j];
    }
    result->b = x->b->coefsT[index];
}


EXPORT void tLweExtractLweSample(LweSample* result, const TLweSample* x, const LweParams* params,  const TLweParams* rparams) {
    tLweExtractLweSampleIndex(result, x, 0, params, rparams);
}

// EXPORT TorusPolynomial* alloc_TorusPolynomial_array(int32_t nbelts) {
//     return (TorusPolynomial*) malloc(nbelts*sizeof(TorusPolynomial));
// }

// TorusPolynomial::TorusPolynomial(const int32_t N): N(N)
// {
//     this->coefsT = new Torus32[N]; 
// }

// EXPORT void init_TorusPolynomial_array(int32_t nbelts, TorusPolynomial* obj, const int32_t N) {
//     for (int32_t i=0; i<nbelts; i++) {
// 	new(obj+i) TorusPolynomial(N);
//     }
// }

// void init_TorusPolynomial_array(int32_t nbelts, struct TorusPolynomial* obj, const int32_t N) {
//     for (int32_t i = 0; i < nbelts; i++) {
//         TorusPolynomial_init(&obj[i], N);
//     }
// }



// EXPORT TorusPolynomial* new_TorusPolynomial_array(int32_t nbelts, const int32_t N) {
//     TorusPolynomial* obj = alloc_TorusPolynomial_array(nbelts);
//     init_TorusPolynomial_array(nbelts,obj,N);
//     return obj;
// }



EXPORT void tfhe_blindRotateAndExtract(LweSample *result,
                                       const TorusPolynomial *v,
                                       const TGswSample *bk,
                                       const int32_t barb,
                                       const int32_t *bara,
                                       const int32_t n,
                                       const TGswParams *bk_params) {

    const TLweParams *accum_params = bk_params->tlwe_params;
    const LweParams *extract_params = &accum_params->extracted_lweparams;
    const int32_t N = accum_params->N;
    const int32_t _2N = 2 * N;

    TorusPolynomial *testvectbis = new_TorusPolynomial(N);
    TLweSample *acc = new_TLweSample(accum_params);

    if (barb != 0) torusPolynomialMulByXai(testvectbis, _2N - barb, v);
    else torusPolynomialCopy(testvectbis, v);
    tLweNoiselessTrivial(acc, testvectbis, accum_params);
    tfhe_blindRotate(acc, bk, bara, n, bk_params);
    tLweExtractLweSample(result, acc, extract_params, accum_params);

    //delete_TLweSample(acc);
    //delete_TorusPolynomial(testvectbis);
}

//tfhe_bootstrap_woKS  



EXPORT void tfhe_bootstrap_woKS(LweSample *result,
                                const LweBootstrappingKey *bk,
                                Torus32 mu, const LweSample *x) {

    const TGswParams *bk_params = bk->bk_params;
    const TLweParams *accum_params = bk->accum_params;
    const LweParams *in_params = bk->in_out_params;
    const int32_t N = accum_params->N;
    const int32_t Nx2 = 2 * N;
    const int32_t n = in_params->n;

    TorusPolynomial *testvect = new_TorusPolynomial(N);
    //int32_t *bara = new int32_t[N];
    int32_t* bara = (int32_t*)malloc(N * sizeof(int32_t));

    int32_t barb = modSwitchFromTorus32(x->b, Nx2);
    for (int32_t i = 0; i < n; i++) {
        bara[i] = modSwitchFromTorus32(x->a[i], Nx2);
    }

    //the initial testvec = [mu,mu,mu,...,mu]
    for (int32_t i = 0; i < N; i++) testvect->coefsT[i] = mu;

    tfhe_blindRotateAndExtract(result, testvect, bk->bk, barb, bara, n, bk_params);

    //delete[] bara;
    //delete_TorusPolynomial(testvect);
}



//key switching
void lweKeySwitchTranslate_fromArray(LweSample* result, 
	const LweSample*** ks, const LweParams* params, 
	const Torus32* ai, 
	const int32_t n, const int32_t t, const int32_t basebit){
    const int32_t base=1<<basebit;       // base=2 in [CGGI16]
    const int32_t prec_offset=1<<(32-(1+basebit*t)); //precision
    const int32_t mask=base-1;

    for (int32_t i=0;i<n;i++){
	const uint32_t aibar=ai[i]+prec_offset;
	for (int32_t j=0;j<t;j++){
	    const uint32_t aij=(aibar>>(32-(j+1)*basebit)) & mask;
	    if(aij != 0) {lweSubTo(result,&ks[i][j][aij],params);}
	}
    }
}
//sample=(a',b')
EXPORT void lweKeySwitch(LweSample* result, const LweKeySwitchKey* ks, const LweSample* sample){
    const LweParams* params=ks->out_params;
    const int32_t n=ks->n;
    const int32_t basebit=ks->basebit;
    const int32_t t=ks->t;

    lweNoiselessTrivial(result,sample->b,params);
    lweKeySwitchTranslate_fromArray(result,
	    (const LweSample***) ks->ks, params,
	    sample->a, n, t, basebit);
}



struct LweSample * new_LweSample( const struct LweParams* params) {
    struct LweSample* sample;
    sample->a = (Torus32*)malloc(params->n * sizeof(Torus32));
    // 初始化其他成员
    sample->b = 0;
    sample->current_variance = 0.0;
    return sample;
}


EXPORT void tfhe_bootstrap(LweSample *result,
                           const LweBootstrappingKey *bk,
                           Torus32 mu, const LweSample *x) {

    struct LweSample *u = new_LweSample(&bk->accum_params->extracted_lweparams);
    // struct LweSample *u = new_gate_bootstrapping_ciphertext_array(1, bk->accum_params->extracted_lweparams);
    tfhe_bootstrap_woKS(u, bk, mu, x);
    // Key Switching
    lweKeySwitch(result, bk->ks, u);

    //delete_LweSample(u);
}

EXPORT Torus32 modSwitchToTorus32(int32_t mu, int32_t Msize){
    uint64_t interv = ((UINT64_C(1)<<63)/Msize)*2; // width of each intervall
    uint64_t phase64 = mu*interv;
    //floor to the nearest multiples of interv
    return phase64>>32;
}


// void LweSample_init(struct LweSample* sample, const struct LweParams* params) {
//     sample->a = (Torus32*)malloc(params->n * sizeof(Torus32));
//     // 初始化其他成员
//     sample->b = 0;
//     sample->current_variance = 0.0;
// }

void LweSample_init(struct LweSample* sample, const struct LweParams* params) {
    sample->a = (Torus32*)malloc(params->n * sizeof(Torus32));
    sample->b = 0;
    sample->current_variance = 0.0;
}

// 分配一个包含nbelts个LweSample的数组
struct LweSample* alloc_LweSample_array(int32_t nbelts) {
    return (struct LweSample*)malloc(nbelts * sizeof(struct LweSample));
}

// 初始化LweSample数组
void init_LweSample_array(int32_t nbelts, struct LweSample* obj, const struct LweParams* params) {
    for (int32_t i = 0; i < nbelts; i++) {
        LweSample_init(&obj[i], params);
    }
}

// 创建并初始化LweSample数组
struct LweSample* new_LweSample_array(int32_t nbelts, const struct LweParams* params) {
    struct LweSample* obj = alloc_LweSample_array(nbelts);
    init_LweSample_array(nbelts, obj, params);
    return obj;
}


// TFheGateBootstrappingParameterSet *read_new_tfheGateBootstrappingParameters(const Istream &F) {
//     int32_t ks_t, ks_basebit;
//     read_tfheGateBootstrappingProperParameters_section(F, ks_t, ks_basebit);
//     LweParams *in_out_params = read_new_lweParams(F);
//     TGswParams *bk_params = read_new_tGswParams(F);
//     TfheGarbageCollector::register_param(in_out_params);
//     TfheGarbageCollector::register_param(bk_params);
//     return new TFheGateBootstrappingParameterSet(ks_t, ks_basebit, in_out_params, bk_params);
// }


// TFheGateBootstrappingSecretKeySet *
// read_new_tfheGateBootstrappingSecretKeySet(const Istream &F, const TFheGateBootstrappingParameterSet *params = 0) {
//     if (params == 0) {
//         TFheGateBootstrappingParameterSet *tmp = read_new_tfheGateBootstrappingParameters(F);
//         TfheGarbageCollector::register_param(tmp);
//         params = tmp;
//     }
//     LweBootstrappingKey *bk = read_new_lweBootstrappingKey(F, params->in_out_params, params->tgsw_params);
//     LweKey *lwe_key = read_new_lweKey(F, params->in_out_params);
//     TGswKey *tgsw_key = read_new_tGswKey(F, params->tgsw_params);
//     LweBootstrappingKeyFFT *bkFFT = new_LweBootstrappingKeyFFT(bk);
//     return new TFheGateBootstrappingSecretKeySet(params, bk, bkFFT, lwe_key, tgsw_key);
// }


struct TGswParams*  TGswParams_init(struct TGswParams* params,int32_t l, int32_t Bgbit, const struct TLweParams* tlwe_params) {
    params->l = l;
    params->Bgbit = Bgbit;
    params->Bg = 1 << Bgbit;
    params->halfBg = params->Bg / 2;
    params->maskMod = params->Bg - 1;
    params->tlwe_params = tlwe_params;
    params->kpl = (tlwe_params->k + 1) * l;

    params->h = (Torus32*)malloc(l * sizeof(Torus32));
    for (int32_t i = 0; i < l; ++i) {
        int32_t kk = (32 - (i + 1) * Bgbit);
        params->h[i] = 1 << kk; // 1/(Bg^(i+1)) as a Torus32
    }

    // offset = Bg/2 * (2^(32-Bgbit) + 2^(32-2*Bgbit) + ... + 2^(32-l*Bgbit))
    uint32_t temp1 = 0;
    for (int32_t i = 0; i < l; ++i) {
        uint32_t temp0 = 1 << (32 - (i + 1) * Bgbit);
        temp1 += temp0;
    }
    params->offset = temp1 * params->halfBg;
    return params;
}


EXPORT void lweKeyGen(LweKey* result) {
  const int32_t n = result->params->n;
  //uniform_int_distribution<int32_t> distribution(0,1);

  for (int32_t i=0; i<n; i++) 
    // result->key[i]=distribution(generator);
    result->key[i]=i;

}



/** generate a tgsw key (in fact, a tlwe key) */
EXPORT void tGswKeyGen(TGswKey *result) {
    tLweKeyGen(&result->tlwe_key);
}
// TLwe
EXPORT void tLweKeyGen(TLweKey *result) {
    const int32_t N = result->params->N;
    const int32_t k = result->params->k;
    //uniform_int_distribution<int32_t> distribution(0, 1);

    for (int32_t i = 0; i < k; ++i)
        for (int32_t j = 0; j < N; ++j)
            //result->key[i].coefs[j] = distribution(generator);
            result->key[i].coefs[j] = i+j;
}


struct LweKey* new_LweKey( const struct LweParams* params) {
    struct LweKey* obj;
    obj->params = params;
    obj->key = (int32_t*)malloc(params->n * sizeof(int32_t));
    return obj;
}

struct TGswKey* new_TGswKey( const struct TGswParams* params) {
    struct TGswKey* obj;
    obj->params = params;
    obj->tlwe_params = params->tlwe_params;
    obj->key = (struct IntPolynomial*)malloc(params->kpl * sizeof(struct IntPolynomial));
    // 初始化其他成员
    return obj;
}
struct TGswSample * new_TGswSample( const struct TGswParams *params) {
    
    const int32_t k = params->tlwe_params->k;
    const int32_t l = params->l;

    struct TLweSample *all_sample = new_TLweSample_array((k + 1) * l, params->tlwe_params);

    struct TLweSample **bloc_sample = (struct TLweSample **)malloc((k + 1) * sizeof(struct TLweSample *));
    for (int32_t p = 0; p < k + 1; ++p) {
        bloc_sample[p] = all_sample + p * l;
    }
    struct TGswSample * obj = {params, bloc_sample, all_sample};
    return obj;
}
struct TGswSample *new_TGswSample_array(int32_t nbelts, const struct TGswParams *params) {
    struct TGswSample *obj = (struct TGswSample *)malloc(nbelts * sizeof(struct TGswSample));
    for (int32_t i = 0; i < nbelts; i++) {
        obj[i] = *new_TGswSample(params);
    }
    return obj;
}

void init_LweKeySwitchKey(struct LweKeySwitchKey *obj, int32_t n, int32_t t, int32_t basebit, const struct LweParams *out_params) {
    const int32_t base = 1 << basebit;
    struct LweSample *ks0_raw = new_LweSample_array(n * t * base, out_params);

    obj->n = n;
    obj->t = t;
    obj->basebit = basebit;
    obj->out_params = out_params;
    obj->ks0_raw = ks0_raw;
}

EXPORT LweKeySwitchKey* alloc_LweKeySwitchKey() {
    return (LweKeySwitchKey*) malloc(sizeof(LweKeySwitchKey));
}
EXPORT LweKeySwitchKey* new_LweKeySwitchKey(int32_t n, int32_t t, int32_t basebit, const LweParams* out_params) {
    LweKeySwitchKey* obj = alloc_LweKeySwitchKey();
    init_LweKeySwitchKey(obj, n,t,basebit,out_params);
    return obj;
}


EXPORT void init_LweBootstrappingKey(LweBootstrappingKey *obj, int32_t ks_t, int32_t ks_basebit, const LweParams *in_out_params,
                                     const TGswParams *bk_params) {
    const TLweParams *accum_params = bk_params->tlwe_params;
    const LweParams *extract_params = &accum_params->extracted_lweparams;
    const int32_t n = in_out_params->n;
    const int32_t N = extract_params->n;

    TGswSample *bk = new_TGswSample_array(n, bk_params);
    LweKeySwitchKey *ks = new_LweKeySwitchKey(N, ks_t, ks_basebit, in_out_params);
    //LweBootstrappingKey *obj = {in_out_params, bk_params, accum_params, extract_params, bk, ks};
    //new(obj) LweBootstrappingKey(in_out_params, bk_params, accum_params, extract_params, bk, ks);
    obj->in_out_params = in_out_params;
    obj->bk_params = bk_params;
    obj->accum_params = accum_params;
    obj->extract_params = extract_params;
    obj->bk = bk;
    obj->ks = ks;
}

EXPORT LweBootstrappingKey *alloc_LweBootstrappingKey() {
    return (LweBootstrappingKey *) malloc(sizeof(LweBootstrappingKey));
}

EXPORT LweBootstrappingKey *
new_LweBootstrappingKey(const int32_t ks_t, const int32_t ks_basebit, const LweParams *in_out_params,
                        const TGswParams *bk_params) {
    LweBootstrappingKey *obj = alloc_LweBootstrappingKey();
    init_LweBootstrappingKey(obj, ks_t, ks_basebit, in_out_params, bk_params);
    return obj;
}

EXPORT void tLweExtractKey(LweKey* result, const TLweKey* key) //sans doute un param supplémentaire
{
    const int32_t N = key->params->N;
    const int32_t k = key->params->k;
    //assert(result->params->n == k*N);
    for (int32_t i=0; i<k; i++) {
	for (int32_t j=0; j<N; j++)
	    result->key[i*N+j]=key->key[i].coefs[j];
    }
}


// EXPORT void lweNoiselessTrivial(LweSample* result, Torus32 mu, const LweParams* params){
//     const int32_t n = params->n;

//     for (int32_t i = 0; i < n; ++i) result->a[i] = 0;
//     result->b = mu;
//     result->current_variance = 0.;
// }

typedef int int32_t;
EXPORT Torus32 dtot32(double d) {
    return (Torus32)(int64_t)((d - (double)(int64_t)d) * (double)_two32);
    }
EXPORT void lweSymEncryptWithExternalNoise(LweSample* result, Torus32 message, double noise, double alpha, const LweKey* key){
    const int32_t n = key->params->n;

    result->b = message + dtot32(noise); 
    for (int32_t i = 0; i < n; ++i)
    {
        //result->a[i] = uniformTorus32_distrib(generator);
        result->a[i] = i;
        result->b += result->a[i]*key->key[i];
    }

    result->current_variance = alpha*alpha;
}


EXPORT void lweCreateKeySwitchKey(LweKeySwitchKey* result, const LweKey* in_key, const LweKey* out_key){
    const int32_t n = result->n;
    const int32_t t = result->t;
    const int32_t basebit = result->basebit;
    const int32_t base = 1<<basebit;
    const double alpha = out_key->params->alpha_min;
    const int32_t sizeks = n*t*(base-1);
    //const int32_t n_out = out_key->params->n;

    double err = 0;

    // chose a random vector of gaussian noises
    double* noise = (double*)malloc(sizeof(double) * sizeks);
    for (int32_t i = 0; i < sizeks; ++i){
        //normal_distribution<double> distribution(0.,alpha); 
        noise[i] = (double)i / (double)sizeks;;
        err += noise[i];
    }
    // recenter the noises
    err = err/sizeks;
    for (int32_t i = 0; i < sizeks; ++i) noise[i] -= err;


    // generate the ks
    int32_t index = 0; 
    for (int32_t i = 0; i < n; ++i) {
        for (int32_t j = 0; j < t; ++j) {

            // term h=0 as trivial encryption of 0 (it will not be used in the KeySwitching)
            lweNoiselessTrivial(&result->ks[i][j][0], 0, out_key->params);
            //lweSymEncrypt(&result->ks[i][j][0],0,alpha,out_key);

            for (int32_t h = 1; h < base; ++h) { // pas le terme en 0
                /*
                // noiseless encryption
                result->ks[i][j][h].b = (in_key->key[i]*h)*(1<<(32-(j+1)*basebit));
                for (int32_t p = 0; p < n_out; ++p) {
                    result->ks[i][j][h].a[p] = uniformTorus32_distrib(generator);
                    result->ks[i][j][h].b += result->ks[i][j][h].a[p] * out_key->key[p];
                }
                // add the noise 
                result->ks[i][j][h].b += dtot32(noise[index]);
                */
                Torus32 mess = (in_key->key[i]*h)*(1<<(32-(j+1)*basebit));
                lweSymEncryptWithExternalNoise(&result->ks[i][j][h], mess, noise[index], alpha, out_key);
                index += 1;
            }
        }
    }


    //delete[] noise; 
}

EXPORT void torusPolynomialUniform(TorusPolynomial *result) {
    const int32_t N = result->N;
    Torus32 *x = result->coefsT;

    for (int32_t i = 0; i < N; ++i)
        // x[i] = uniformTorus32_distrib(generator);
        x[i] = i;
}



// EXPORT void torusPolynomialAddMulRFFT(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2) {
//     const int32_t N = poly1->N;
//     LagrangeHalfCPolynomial* tmp = new_LagrangeHalfCPolynomial_array(3,N);
//     TorusPolynomial* tmpr = new_TorusPolynomial(N);
//     IntPolynomial_ifft(tmp+0,poly1);
//     TorusPolynomial_ifft(tmp+1,poly2);
//     LagrangeHalfCPolynomialMul(tmp+2,tmp+0,tmp+1);
//     TorusPolynomial_fft(tmpr, tmp+2);
//     torusPolynomialAddTo(result, tmpr);
//     delete_TorusPolynomial(tmpr);
//     delete_LagrangeHalfCPolynomial_array(3,tmp);
// }

EXPORT void tLweSymEncryptZero(TLweSample *result, double alpha, const TLweKey *key) {
    const int32_t N = key->params->N;
    const int32_t k = key->params->k;

    for (int32_t j = 0; j < N; ++j)
        // result->b->coefsT[j] = gaussian32(0, alpha);
        result->b->coefsT[j] = (double)j / (double)N;

    for (int32_t i = 0; i < k; ++i) {
        torusPolynomialUniform(&result->a[i]);
        torusPolynomialAddMulR(result->b, &key->key[i], &result->a[i]);
    }

    result->current_variance = alpha * alpha;
}

EXPORT void tGswEncryptZero(TGswSample *result, double alpha, const TGswKey *key) {
    const TLweKey *rlkey = &key->tlwe_key;
    const int32_t kpl = key->params->kpl;

    for (int32_t p = 0; p < kpl; ++p) {
        tLweSymEncryptZero(&result->all_sample[p], alpha, rlkey);
    }
}

EXPORT void tGswAddMuIntH(TGswSample *result, const int32_t message, const TGswParams *params) {
    const int32_t k = params->tlwe_params->k;
    const int32_t l = params->l;
    const Torus32 *h = params->h;

    // compute result += H
    for (int32_t bloc = 0; bloc <= k; ++bloc)
        for (int32_t i = 0; i < l; i++)
            result->bloc_sample[bloc][i].a[bloc].coefsT[0] += message * h[i];
}

EXPORT void tGswSymEncryptInt(TGswSample *result, const int32_t message, double alpha, const TGswKey *key) {
    tGswEncryptZero(result, alpha, key);
    tGswAddMuIntH(result, message, key->params);
}

EXPORT void tfhe_createLweBootstrappingKey(
        LweBootstrappingKey *bk,
        const LweKey *key_in,
        const TGswKey *rgsw_key) {
    // assert(bk->bk_params == rgsw_key->params);
    // assert(bk->in_out_params == key_in->params);

    const LweParams *in_out_params = bk->in_out_params;
    const TGswParams *bk_params = bk->bk_params;
    const TLweParams *accum_params = bk_params->tlwe_params;
    const LweParams *extract_params = &accum_params->extracted_lweparams;

    //LweKeySwitchKey* ks; ///< the keyswitch key (s'->s)
    const TLweKey *accum_key = &rgsw_key->tlwe_key;
    LweKey *extracted_key = new_LweKey(extract_params);
    tLweExtractKey(extracted_key, accum_key);
    lweCreateKeySwitchKey(bk->ks, extracted_key, key_in);
    //delete_LweKey(extracted_key);

    //TGswSample* bk; ///< the bootstrapping key (s->s")
    int32_t *kin = key_in->key;
    const double alpha = accum_params->alpha_min;
    const int32_t n = in_out_params->n;
    //const int32_t kpl = bk_params->kpl;
    //const int32_t k = accum_params->k;
    //const int32_t N = accum_params->N;
    //cout << "create the bootstrapping key bk ("  << "  " << n*kpl*(k+1)*N*4 << " bytes)" << endl;
    //cout << "  with noise_stdev: " << alpha << endl;
    for (int32_t i = 0; i < n; i++) {
        tGswSymEncryptInt(&bk->bk[i], kin[i], alpha, rgsw_key);
    }

}


EXPORT TFheGateBootstrappingSecretKeySet *
new_random_gate_bootstrapping_secret_keyset(const TFheGateBootstrappingParameterSet *params) {
    LweKey *lwe_key = new_LweKey(params->in_out_params);
    lweKeyGen(lwe_key);
    TGswKey *tgsw_key = new_TGswKey(params->tgsw_params);
    tGswKeyGen(tgsw_key);
    LweBootstrappingKey *bk = new_LweBootstrappingKey(params->ks_t, params->ks_basebit, params->in_out_params,
                                                      params->tgsw_params);
    tfhe_createLweBootstrappingKey(bk, lwe_key, tgsw_key);
    LweBootstrappingKeyFFT *bkFFT ;
    // = new_LweBootstrappingKeyFFT(bk);
    //return new TFheGateBootstrappingSecretKeySet(params, bk, bkFFT, lwe_key, tgsw_key);
    TFheGateBootstrappingSecretKeySet * key = {params, bk, bkFFT, lwe_key, tgsw_key};
    return key;
}




int main()
{   //we need note the  current ptr
    //struct TFheGateBootstrappingParameterSet *params;
    FILE * prams_file = fopen("size_prams","wb");
    //params = new_tfheGateBootstrappingParameterSet_fromFile(prams_file);
    int std_params = 1;

    // if (std_params) {
    //////
    static const int32_t N = 1024;//the prime value is 1024
    static const int32_t k = 1;
    static const int32_t n = 20;
    static const int32_t bk_l = 2;
    static const int32_t bk_Bgbit = 7;
    static const int32_t ks_basebit = 2;
    static const int32_t ks_length = 8;
    //static const double ks_stdev = pow(2.,-15); //standard deviation
    // static const double bk_stdev = pow(2.,-25);; //standard deviation
    static const double ks_stdev = 1.0 / (1 << 15); //standard deviation
    static const double bk_stdev = 1.0 / (1 << 25);; //standard deviation
    static const double max_stdev = 0.012467; //max standard deviation for a 1/4 msg space
    LweParams params_in_temp = {n, ks_stdev, max_stdev};
    LweParams *params_in = &params_in_temp;
    TLweParams params_accum_temp = {N, k, bk_stdev, max_stdev};
    TLweParams *params_accum = &params_accum_temp;

    //TGswParams params_bk_temp = {bk_l, bk_Bgbit, params_accum};
    TGswParams *params_bk = TGswParams_init(params_bk,bk_l, bk_Bgbit, params_accum);
    
    TFheGateBootstrappingParameterSet params_temp = {ks_length, ks_basebit, params_in, params_bk};
    TFheGateBootstrappingParameterSet *params = &params_temp;   
    //////
    // }

  //TFheGateBootstrappingSecretKeySet *key = new_random_gate_bootstrapping_secret_keyset(params);
    FILE * x7_file = fopen("size_bootstrap_secretkey","wb");
    //struct TFheGateBootstrappingSecretKeySet *key = new_tfheGateBootstrappingSecretKeySet_fromFile(x7_file);

    uint32_t seed[] = {214, 1592, 657};
    //tfhe_random_generator_setSeed(seed, 3);
    TFheGateBootstrappingSecretKeySet *key = new_random_gate_bootstrapping_secret_keyset(params);

    const struct TFheGateBootstrappingCloudKeySet *bk = &(key->cloud); 

    struct LweSample *temp_result = new_LweSample_array(1,params->in_out_params);

    const Torus32 MU = modSwitchToTorus32(1, 8);
    FILE* temp_result_file = fopen("temp_result","wb");
    //import_lweSample_fromFile(temp_result_file,temp_result,params->in_out_params);
    //here we must need know what operation we want to do
    struct LweSample *result = new_LweSample_array(1,params->in_out_params);
    tfhe_bootstrap(result,bk->bk,MU,temp_result);

    //bool c_de = bootsSymDecrypt(c_en,key);
    //std::cout<<"output and data c_de :" << c_de << std::endl;
    //bootsOR(c_en,a_en,b_en,bk);



    // time_t t1 = clock();
    //decrypt
    //bool c_de = bootsSymDecrypt(c_en,key);
    //delete_gate_bootstrapping_secret_keyset(key);
    //delete_gate_bootstrapping_parameters(params);
    return 0;
}