#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>



#include "lagrangehalfc.h"
/*
#include "ni.h"
#include "ni.c"
#include "sdram.h"
*/
typedef int32_t Torus32; //avant uint32_t





























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


struct IntPolynomial {
   const int32_t N;
   int32_t* coefs;

#ifdef __cplusplus   
   IntPolynomial(const int32_t N);
   ~IntPolynomial();
   IntPolynomial(const IntPolynomial&) = delete; //forbidden 
   IntPolynomial* operator=(const IntPolynomial&) = delete; //forbidden
#endif
};


/** This structure represents an torus polynomial modulo X^N+1 */
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

struct TLweSample {
    TorusPolynomial *a; ///< array of length k+1: mask + right term
    TorusPolynomial *b; ///< alias of a[k] to get the right term
    double current_variance; ///< avg variance of the sample
    const int32_t k;
#ifdef __cplusplus

    TLweSample(const TLweParams *params);

    ~TLweSample();

    TLweSample(const TLweSample &) = delete;

    void operator=(const TLweSample &) = delete;

#endif
};






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

EXPORT void IntPolynomial_ifft(LagrangeHalfCPolynomial* result, const IntPolynomial* p) {
    fp1024_fftw.execute_reverse_int(((LagrangeHalfCPolynomial_IMPL*)result)->coefsC, p->coefs);
}
EXPORT void TorusPolynomial_ifft(LagrangeHalfCPolynomial* result, const TorusPolynomial* p) {
    fp1024_fftw.execute_reverse_torus32(((LagrangeHalfCPolynomial_IMPL*)result)->coefsC, p->coefsT);
}
EXPORT void TorusPolynomial_fft(TorusPolynomial* result, const LagrangeHalfCPolynomial* p) {
    fp1024_fftw.execute_direct_Torus32(result->coefsT, ((LagrangeHalfCPolynomial_IMPL*)p)->coefsC);
}


EXPORT void torusPolynomialAddMulRFFT(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2) {
    const int32_t N = poly1->N;
    LagrangeHalfCPolynomial* tmp = new_LagrangeHalfCPolynomial_array(3,N);
    TorusPolynomial* tmpr = new_TorusPolynomial(N);
    IntPolynomial_ifft(tmp+0,poly1);
    TorusPolynomial_ifft(tmp+1,poly2);
    LagrangeHalfCPolynomialMul(tmp+2,tmp+0,tmp+1);
    TorusPolynomial_fft(tmpr, tmp+2);
    torusPolynomialAddTo(result, tmpr);
    delete_TorusPolynomial(tmpr);
    delete_LagrangeHalfCPolynomial_array(3,tmp);

//tGswExternMulToTLwe

EXPORT IntPolynomial* alloc_IntPolynomial_array(int32_t nbelts) {
    return (IntPolynomial*) malloc(nbelts*sizeof(IntPolynomial));
}

IntPolynomial::IntPolynomial(const int32_t N): N(N)
{
    this->coefs = new int32_t[N]; 
}

EXPORT void init_IntPolynomial_array(int32_t nbelts, IntPolynomial* obj, const int32_t N) {
    for (int32_t i=0; i<nbelts; i++) {
	new(obj+i) IntPolynomial(N);
    }
}


EXPORT IntPolynomial* new_IntPolynomial_array(int32_t nbelts, const int32_t N) {
    IntPolynomial* obj = alloc_IntPolynomial_array(nbelts);
    init_IntPolynomial_array(nbelts,obj,N);
    return obj;
}



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


#define torusPolynomialAddMulR torusPolynomialAddMulRFFT

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

    delete_IntPolynomial_array(kpl, dec);
}



//tfhe_MuxRotate

//result= (X^{a}-1)*source
EXPORT void torusPolynomialMulByXaiMinusOne(TorusPolynomial *result, int32_t a, const TorusPolynomial *source) {
    const int32_t N = source->N;
    Torus32 *out = result->coefsT;
    Torus32 *in = source->coefsT;

    assert(a >= 0 && a < 2 * N);

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
        swap(temp2, temp3);

    }
    if (temp3 != accum) {
        tLweCopy(accum, temp3, bk_params->tlwe_params);
    }

    delete_TLweSample(temp);
    //delete_TGswSample(temp);
}

//tfhe_blindRotateAndExtract





// TorusPolynomial = TorusPolynomial


//result= X^{a}*source
EXPORT void torusPolynomialMulByXai(TorusPolynomial *result, int32_t a, const TorusPolynomial *source) {
    const int32_t N = source->N;
    Torus32 *out = result->coefsT;
    Torus32 *in = source->coefsT;

    assert(a >= 0 && a < 2 * N);
    assert(result != source);

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
    assert(result != sample);
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
    assert(params->n == k*N);

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

TorusPolynomial::TorusPolynomial(const int32_t N): N(N)
{
    this->coefsT = new Torus32[N]; 
}

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



EXPORT TorusPolynomial* new_TorusPolynomial_array(int32_t nbelts, const int32_t N) {
    TorusPolynomial* obj = alloc_TorusPolynomial_array(nbelts);
    init_TorusPolynomial_array(nbelts,obj,N);
    return obj;
}



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

    delete_TLweSample(acc);
    delete_TorusPolynomial(testvectbis);
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
    int32_t *bara = new int32_t[N];

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



EXPORT void tfhe_bootstrap(LweSample *result,
                           const LweBootstrappingKey *bk,
                           Torus32 mu, const LweSample *x) {

    LweSample *u = new_LweSample(&bk->accum_params->extracted_lweparams);

    tfhe_bootstrap_woKS(u, bk, mu, x);
    // Key Switching
    lweKeySwitch(result, bk->ks, u);

    //delete_LweSample(u);
}






int main()
{   //we need note the  current ptr
    TFheGateBootstrappingParameterSet *params;
    FILE * prams_file = fopen("size_prams","wb");
    params = new_tfheGateBootstrappingParameterSet_fromFile(prams_file);

  //TFheGateBootstrappingSecretKeySet *key = new_random_gate_bootstrapping_secret_keyset(params);
    FILE * x7_file = fopen("size_bootstrap_secretkey","wb");
    TFheGateBootstrappingSecretKeySet *key = new_tfheGateBootstrappingSecretKeySet_fromFile(x7_file);
    const TFheGateBootstrappingCloudKeySet *bk = &(key->cloud); 

    LweSample *temp_result = new_gate_bootstrapping_ciphertext_array(1,params);

    static const Torus32 MU = modSwitchToTorus32(1, 8);
    LweSample *c_en = new_gate_bootstrapping_ciphertext_array(1,params);
    FILE* temp_result_file = fopen("temp_result","wb");
    temp_result = import_lweSample_fromFile(temp_result_file,params->in_out_params,temp_result);
    //here we must need know what operation we want to do
    tfhe_bootstrap(result,bk->bk,MU,temp_result);

    //bool c_de = bootsSymDecrypt(c_en,key);
    std::cout<<"output and data c_de :" << c_de << std::endl;
    //bootsOR(c_en,a_en,b_en,bk);



    // time_t t1 = clock();
    //decrypt
    bool c_de = bootsSymDecrypt(c_en,key);
    //delete_gate_bootstrapping_secret_keyset(key);
    //delete_gate_bootstrapping_parameters(params);
}