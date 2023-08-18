#include <stdio.h>
#include <stdlib.h>
//#include <math.h>
//#include <fftw3.h>

#define N_value 1024
#define n_value 630
#define MEMORY_SIZE 50*1024


#define INT64_C(c) (c ## LL)
#define UINT64_C(c) (c ## ULL)

typedef int32_t Torus32; //avant uint32_t
typedef unsigned int		uint32_t;
typedef unsigned long long uint64_t;
#define  star_disk_address 0x00000000


char memory[MEMORY_SIZE];
uint32_t used_memory = 0;
uint64_t rest_memory = MEMORY_SIZE;
//#include "lagrangehalfc_impl.h"
/*
#include "ni.h"
#include "ni.c"
#include "sdram.h"
*/
// 简化的内存块结构
typedef struct Block {
    size_t size;
    struct Block* next;
} Block;

typedef struct Variable {
    char name[32];
    size_t size;
    void* value;
} Variable;
const char *filename = "sdram.txt";

// 内存池起始地址
static void* memory_start = NULL;
// 空闲内存链表
static Block* free_list = NULL;


char sdram[512*1024*1024];
static void * curr_disk_address = sdram;
static void * sd_start = NULL;
static Block* sd_free_list = NULL;


void initMemoryPool(void* start, size_t size) {
    memory_start = start;
    free_list = (Block*)start;
    free_list->size = size - sizeof(Block);
    free_list->next = NULL;
}





void init_sdram(void* start, size_t size) {
    sd_start = start;
    sd_free_list = (Block*)start;
    sd_free_list->size = size - sizeof(Block);
    sd_free_list->next = NULL;
}


// void initMemoryPool_512MB(void* start, size_t size) {
//     memory_start = start;
//     free_list = (Block*)start;
//     free_list->size = size - sizeof(Block);
//     free_list->next = NULL;
// }


char* myStrncpy(char* dest, const char* src, size_t n) {
    size_t i;
    for (i = 0; i < n && src[i] != '\0'; ++i) {
        dest[i] = src[i];
    }
    for (; i < n; ++i) {
        dest[i] = '\0';
    }
    return dest;
}


int malloc_ram_size =0;
void* mymalloc(size_t size) {
    if (size == 0) {
        return NULL;
    }
    used_memory += size ;
    used_memory += sizeof(struct Block);
    // printf("the size of block is %d\n",sizeof(struct Block));/////// size of block is 16B
    printf("the used memory is %d\n",used_memory);
    // rest_memory -= (size +sizeof(Block));
    // printf("the rest memory is %d",rest_memory);
    // printf("    the rest of memory is %ld  ~~  and ~~  ",rest_memory);
    printf("    malloc size is %d  ~~  and ~~  ", size);
    malloc_ram_size += size;
    printf("    malloc_ram_size is %d\n", malloc_ram_size);
    Block* prev = NULL;
    Block* current = free_list;

    while (current != NULL) {
        if (current->size >= size) {
            if (current->size > size + sizeof(Block)) {
                Block* next = (Block*)((char*)current + size + sizeof(Block));
                next->size = current->size - size - sizeof(Block);
                next->next = current->next;
                current->size = size;
                current->next = next;
            }

            if (prev == NULL) {
                free_list = current->next;
            } else {
                prev->next = current->next;
            }
            //printf("the address of current is %p\n", current);
            return (char*)current + sizeof(Block);
        }
        
        prev = current;
        current = current->next;
    }
    printf("we haven't enough memory\n");
    abort();

    return NULL; // 没有足够的空闲内存
}



void mergeFreeBlocks() {
    Block* prev = NULL;
    Block* current = free_list;
    
    while (current != NULL && current->next != NULL) {
        Block* nextBlock = current->next;
        
        if ((char*)current + current->size + sizeof(Block) == (char*)nextBlock) {
            current->size += nextBlock->size + sizeof(Block);
            current->next = nextBlock->next;
        } else {
            prev = current;
        }
        
        current = current->next;
    }
}

void myFree(void* ptr) {
    if (ptr != NULL) {
        
        Block* block = (Block*)((char*)ptr - sizeof(Block));
        rest_memory += (block->size + sizeof(Block));
        used_memory -= block->size;
        used_memory -= sizeof(struct Block);
        printf("used memory is recovered to %d",used_memory);
        printf("free size is %ld memory and now the rest of memory may be is  %ld\n",block->size,rest_memory);
        block->next = free_list;
        free_list = block;
        mergeFreeBlocks();
    }
}




// #define INT64_C(c)  c ## LL

//#include "tfhe_io.h"
//typedef long long int int64_t;
uint32_t startAddr = 0x1000;



void * write_data_to_disk( void *data, size_t dataSize,uint32_t initialed_value);


#define torusPolynomialMulR torusPolynomialMultFFT
// #define torusPolynomialAddMulR torusPolynomialAddMulRFFT
#define torusPolynomialSubMulR torusPolynomialSubMulRFFT

static const long _two32 = INT64_C(1) << 32; // 2^32



//get the number of multiplications for a bootstrapping
static int mult_count = 0;










struct TorusPolynomial {
    int32_t N;
   Torus32* coefsT;
};



struct LweParams {
	 int32_t n;
	 double alpha_min;//le plus petit bruit tq sur
	 double alpha_max;//le plus gd bruit qui permet le déchiffrement
//since all members are declared constant, a constructor is 
//required in the structure.
};

struct LweSample {
	Torus32* a; //-- the n coefs of the mask
    Torus32 b;  //
   	double current_variance; //-- average noise of the sample
};

struct TLweParams {
     int32_t N; ///< a power of 2: degree of the polynomials
     int32_t k; ///< number of polynomials in the mask
     double alpha_min; ///< minimal noise s.t. the sample is secure
     double alpha_max; ///< maximal noise s.t. we can decrypt
    struct LweParams extracted_lweparams; ///< lwe params if one extracts
};

struct TLweSample {
    struct TorusPolynomial *a; ///< array of length k+1: mask + right term
    struct TorusPolynomial *b; ///< alias of a[k] to get the right term
    double current_variance; ///< avg variance of the sample
    int32_t k;
};

struct IntPolynomial {
    int32_t N;
   int32_t* coefs;
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
};

struct TGswSample {
    struct TLweSample *all_sample; ///< TLweSample* all_sample; (k+1)l TLwe Sample
    struct TLweSample **bloc_sample;///< accès optionnel aux différents blocs de taille l.
    // double current_variance;
     int32_t k;
     int32_t l;
};


struct LweKeySwitchKey {
    int32_t n; ///< length of the input key: s'
    int32_t t; ///< decomposition length
    int32_t basebit; ///< log_2(base)
    int32_t base; ///< decomposition base: a power of 2 
     struct  LweParams* out_params; ///< params of the output key s 
    struct LweSample* ks0_raw; //tableau qui contient tout les Lwe samples de taille nlbase
    struct LweSample** ks1_raw;// de taille nl  pointe vers un tableau ks0_raw dont les cases sont espaceés de base positions
    struct LweSample*** ks; ///< the keyswitch elements: a n.l.base matrix
    // de taille n pointe vers ks1 un tableau dont les cases sont espaceés de ell positions
};

struct LweBootstrappingKey{
     struct LweParams* in_out_params; ///< paramètre de l'input et de l'output. key: s
      struct TGswParams* bk_params; ///< params of the Gsw elems in bk. key: s"
     struct TLweParams* accum_params; ///< params of the accum variable key: s"
     struct LweParams* extract_params; ///< params after extraction: key: s' 
    struct TGswSample* bk; ///< the bootstrapping key (s->s")
    struct LweKeySwitchKey* ks; ///< the keyswitch key (s'->s)
};


struct TFheGateBootstrappingParameterSet {
     int32_t ks_t;
     int32_t ks_basebit;
     struct LweParams * in_out_params;
     struct TGswParams * tgsw_params;
};

struct LweKey {
    struct LweParams* params;
   int32_t* key;
};
struct TLweKey {
     struct TLweParams *params; ///< the parameters of the key
    struct IntPolynomial *key; ///< the key (i.e k binary polynomials)
};

struct TGswKey {
     struct TGswParams *params; ///< the parameters
     struct TLweParams *tlwe_params; ///< the tlwe params of each rows
    struct IntPolynomial *key; ///< the key (array of k polynomials)
    struct TLweKey tlwe_key;
};
// struct TLweSampleFFT {
//     struct LagrangeHalfCPolynomial *a; ///< array of length k+1: mask + right term
//     struct LagrangeHalfCPolynomial *b; ///< alias of a[k] to get the right term
//     double current_variance; ///< avg variance of the sample
//      int32_t k; //required during the destructor call...
// #ifdef __cplusplus

//     TLweSampleFFT(const TLweParams *params, LagrangeHalfCPolynomial *a, double current_variance);

//     ~TLweSampleFFT();

//     TLweSampleFFT(const TLweSampleFFT &) = delete;

//     void operator=(const TLweSampleFFT &) = delete;

// #endif
// };
// struct TGswSampleFFT {
//     struct TLweSampleFFT *all_samples; ///< TLweSample* all_sample; (k+1)l TLwe Sample
//     struct TLweSampleFFT **sample; ///< accès optionnel aux différents blocs de taille l.
//     //double current_variance;
//      int32_t k;
//      int32_t l;

// #ifdef __cplusplus

//     TGswSampleFFT(const TGswParams *params, TLweSampleFFT *all_samples);

//     ~TGswSampleFFT();

//     TGswSampleFFT(const TGswSampleFFT &) = delete;

//     void operator=(const TGswSampleFFT &) = delete;

// #endif
// };
// struct LweBootstrappingKeyFFT {
//     const struct LweParams* in_out_params; ///< paramètre de l'input et de l'output. key: s
//     const struct TGswParams* bk_params; ///< params of the Gsw elems in bk. key: s"
//     const struct TLweParams* accum_params; ///< params of the accum variable key: s"
//     const struct LweParams* extract_params; ///< params after extraction: key: s' 
//     const struct TGswSampleFFT* bkFFT; ///< the bootstrapping key (s->s")
//     const struct LweKeySwitchKey* ks; ///< the keyswitch key (s'->s)


// #ifdef __cplusplus
//    LweBootstrappingKeyFFT(const LweParams* in_out_params, 
//     const TGswParams* bk_params,
//     const TLweParams* accum_params,
//     const LweParams* extract_params, 
//     const TGswSampleFFT* bkFFT,
//     const LweKeySwitchKey* ks);
//     ~LweBootstrappingKeyFFT();
//     LweBootstrappingKeyFFT(const LweBootstrappingKeyFFT&) = delete;
//     void operator=(const LweBootstrappingKeyFFT&) = delete;
  
// #endif


// };

struct TFheGateBootstrappingCloudKeySet {
     struct TFheGateBootstrappingParameterSet * params;
     struct LweBootstrappingKey * bk;
    //const LweBootstrappingKeyFFT *const bkFFT;
};
struct TFheGateBootstrappingCloudKeySet * new_TFheGateBootstrappingCloudKeySet(const struct TFheGateBootstrappingParameterSet *params,struct  LweBootstrappingKey *bk){
    struct TFheGateBootstrappingCloudKeySet *result = (struct TFheGateBootstrappingCloudKeySet *)  mymalloc(sizeof(struct TFheGateBootstrappingCloudKeySet));
    result->params = params;
    result->bk = bk;
    return result;
}


struct TFheGateBootstrappingSecretKeySet {
     struct TFheGateBootstrappingParameterSet *params;
     struct LweKey *lwe_key;
     struct TGswKey *tgsw_key;
     struct TFheGateBootstrappingCloudKeySet cloud;
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
// void TorusPolynomial_init(TorusPolynomial* poly,  int32_t N) {
//     poly->N = N;
//     poly->coefsT = (Torus32*) mymalloc(N * sizeof(Torus32));
//     write_data_to_disk(poly->coefsT, sizeof(Torus32)*N);
//     myFree(poly->coefsT);
//     // write_data_to_disk(poly, sizeof(TorusPolynomial));
// }

// // TorusPolynomial的清理函数（类似于析构函数）
// void TorusPolynomial_destroy(struct TorusPolynomial* poly) {
//     free(poly->coefsT);
// }








//#include "tfhe.h"
//#include "tfhe_io.h"


//  wle functionfnew_TorusPolynomial_array

  void lweNoiselessTrivial(LweSample* result, Torus32 mu, const LweParams* params){
    // printf("2.1\n");
    const int32_t n = params->n;
    // printf("2.5\n");
    for (int32_t i = 0; i < n; ++i) result->a[i] = 0;
    result->b = mu;
    result->current_variance = 0.;
}


/** result = result - sample */
  void lweSubTo(LweSample* result, const LweSample* sample, const LweParams* params){
    const int32_t n = params->n;
    const Torus32* __restrict sa = sample->a;
    Torus32* __restrict ra = result->a;
    result->b -= sample->b;
    result->current_variance += sample->current_variance; 
}



  void lweAddTo(LweSample* result, const LweSample* sample, const LweParams* params){
    const int32_t n = params->n;

    for (int32_t i = 0; i < n; ++i) result->a[i] += sample->a[i];
    result->b += sample->b;
    result->current_variance += sample->current_variance; 
}


//typedef unsigned long uint64_t1;



int32_t modSwitchFromTorus32(Torus32 phase1, int32_t Msize){
    uint64_t interv = ((UINT64_C(1)<<63)/Msize)*2; // width of each intervall
    uint64_t half_interval = interv/2; // begin of the first intervall
    uint64_t phase64 = (uint64_t)phase1;
    phase64 = (phase64 << 32) + half_interval;
   // << 32) + half_interval;
    //floor to the nearest multiples of interv
    return phase64/interv;
}

// return the disk address of the TorusPolynomial
  TorusPolynomial* new_TorusPolynomial(const int32_t N) {
    TorusPolynomial* tmp = (TorusPolynomial*)  mymalloc(sizeof(TorusPolynomial));


    //tmp = (TorusPolynomial*)malloc(sizeof(TorusPolynomial));
    //tmp = (TorusPolynomial*)write_data_to_disk(tmp_ram, sizeof(TorusPolynomial),1);

    tmp->N = N;
    //tmp->coefsT = (Torus32*)mymalloc(N*sizeof(Torus32));
    tmp->coefsT = write_data_to_disk(tmp->coefsT, sizeof(Torus32)*N,0);
    //myFree(tmp->coefsT);
    myFree(tmp);
    tmp = write_data_to_disk(tmp, sizeof(TorusPolynomial),1);
    // tmp->coefsT = (Torus32*)malloc(N*sizeof(Torus32));
    return tmp;
}

  TorusPolynomial* new_TorusPolynomial_array(int32_t nbelts, const int32_t N) {
    TorusPolynomial* obj = (TorusPolynomial*)mymalloc(nbelts*sizeof(TorusPolynomial));
    myFree(obj);
    //    if(obj == NULL) 
    //     printf("error !!!!failed  mymalloc a TorusPolynomial array\n "); 
    //obj = (TorusPolynomial*)malloc(nbelts*sizeof(TorusPolynomial));
    //obj = (TorusPolynomial*)write_data_to_disk(obj, sizeof(TorusPolynomial)*nbelts,0);

    for (int32_t i = 0; i < nbelts; i++) {
        TorusPolynomial* temp =new_TorusPolynomial(N);
        //write_data_to_disk(temp, sizeof(TorusPolynomial));
        obj[i] = *temp;
        //myFree(temp);
        //obj[i] = (TorusPolynomial*)write_data_to_disk(&obj[i], sizeof(TorusPolynomial),1);
        //myFree(obj);
        // myFree(obj[i].coefsT);
        //myFree(obj[i].coefsT);
    }
    // init_TorusPolynomial_array(nbelts,obj,N);
    // write_data_to_disk( obj, sizeof(TorusPolynomial));
    obj = (TorusPolynomial*)write_data_to_disk(obj, sizeof(TorusPolynomial)*nbelts,1);
    // myFree(obj);
    return obj;
}


//   void init_LagrangeHalfCPolynomial_array(int32_t nbelts, LagrangeHalfCPolynomial* obj, const int32_t N) {
//     for (int32_t i=0; i<nbelts; i++) {
// 	//  LagrangeHalfCPolynomial_IMPL_init(&(obj[i]),N);
//     }
//     }
//   struct LagrangeHalfCPolynomial* new_LagrangeHalfCPolynomial_array(int32_t nbelts, const int32_t N) {
//     struct LagrangeHalfCPolynomial* obj =(struct LagrangeHalfCPolynomial*)  mymalloc(nbelts*sizeof(struct LagrangeHalfCPolynomial));
//     init_LagrangeHalfCPolynomial_array(nbelts,obj,N);
//     return obj;
// }


struct TLweSample *new_TLweSample(const TLweParams* params) {
    const int32_t k = params->k;
    TLweSample* result = (TLweSample*)mymalloc(sizeof(TLweSample));
    result->a = new_TorusPolynomial_array(k+1, params->N);
    result->b = result->a + k;
    result->current_variance = 0.;
    myFree(result);

    result = (TLweSample*)write_data_to_disk(result, sizeof(TLweSample),1);
    return result;
}
struct TLweSample *new_TLweSample_array(int32_t nbelts, const struct TLweParams *params) {
    struct TLweSample *obj = (struct TLweSample *) mymalloc(nbelts * sizeof(struct TLweSample));
    for (int32_t i = 0; i < nbelts; i++) {
        struct TLweSample* temp = new_TLweSample(params);
        obj[i] = *temp;
        //myFree(temp);
    }
    myFree(obj);
    obj = (struct TLweSample *)write_data_to_disk(obj, sizeof(struct TLweSample)*nbelts,1);

    return obj;
}


// struct TLweSampleFFT *new_TLweSampleFFT(const TLweParams* params) {
//     const int32_t k = params->k;
//     TLweSampleFFT* result = (TLweSampleFFT*)  mymalloc(sizeof(TLweSampleFFT));
//     // result->a = new_TorusPolynomialFFT_array(k+1, params->N);
//     result->a = new_LagrangeHalfCPolynomial_array(k+1, params->N);
//     result->b = result->a + k;
//     result->current_variance = 0.;
//     return result;
// }
// struct TLweSampleFFT *new_TLweSampleFFT_array(int32_t nbelts, const struct TLweParams *params) {
//     struct TLweSampleFFT *obj = (struct TLweSampleFFT *) mymalloc(nbelts * sizeof(struct TLweSampleFFT));
//     for (int32_t i = 0; i < nbelts; i++) {
//         obj[i] = *new_TLweSampleFFT(params);
//     }
//     return obj;
// }

// TorusPolynomial = 0
  void torusPolynomialClear(TorusPolynomial *result) {
    const int32_t N = result->N;

    for (int32_t i = 0; i < N; ++i) result->coefsT[i] = 0;
}

//   void tfhe_blindRotate(TLweSample* accum, const TGswSample* bk, const int32_t* bara, const int32_t n, const TGswParams* bk_params);
//   void tfhe_blindRotateAndExtract(LweSample* result, const TorusPolynomial* v, const TGswSample* bk, const int32_t barb, const int32_t* bara, const int32_t n, const TGswParams* bk_params);
//   void tfhe_bootstrap_woKS(LweSample* result, const LweBootstrappingKey* bk, Torus32 mu, const LweSample* x);
//   void tfhe_bootstrap(LweSample* result, const LweBootstrappingKey* bk, Torus32 mu, const LweSample* x);
//   void tfhe_createLweBootstrappingKey(LweBootstrappingKey* bk, const LweKey* key_in, const TGswKey* rgsw_key);




  //void die_dramatically(const char* message);

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



//   void IntPolynomial_ifft(LagrangeHalfCPolynomial* result, const IntPolynomial* p) {
//     fp1024_fftw.execute_reverse_int = execute_reverse_int;
//     fp1024_fftw.execute_reverse_int(((struct LagrangeHalfCPolynomial_IMPL*)result)->coefsC, p->coefs);
// }
//   void TorusPolynomial_ifft(LagrangeHalfCPolynomial* result, const TorusPolynomial* p) {
//     fp1024_fftw.execute_reverse_torus32 = execute_reverse_torus32;  
//     fp1024_fftw.execute_reverse_torus32(((struct LagrangeHalfCPolynomial_IMPL*)result)->coefsC, p->coefsT);
// }
//   void TorusPolynomial_fft(TorusPolynomial* result, const LagrangeHalfCPolynomial* p) {
//     fp1024_fftw.execute_direct_Torus32 = execute_direct_Torus32;
//     fp1024_fftw.execute_direct_Torus32(result->coefsT, ((struct LagrangeHalfCPolynomial_IMPL*)p)->coefsC);
// }



// void LagrangeHalfCPolynomial_IMPL_init(struct LagrangeHalfCPolynomial_IMPL* obj, const int32_t N) {
//     ////assert(N == 1024);
//     //obj->N = N;

//     // 分配内存并初始化 coefsC 数组
//     obj->coefsC = (struct cplx*) mymalloc((N / 2) * sizeof(struct cplx));
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




// TorusPolynomial += TorusPolynomial
  void torusPolynomialAddTo(TorusPolynomial *result, const TorusPolynomial *poly2) {
    const int32_t N = poly2->N;
    Torus32 *r = result->coefsT;
    const Torus32 *b = poly2->coefsT;

    for (int32_t i = 0; i < N; ++i)
        r[i] += b[i];
}


//   void torusPolynomialAddMulRFFT(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2) {
//     return;
//     const int32_t N = poly1->N;
//     struct LagrangeHalfCPolynomial* tmp = new_LagrangeHalfCPolynomial_array(3,N);
//     struct TorusPolynomial* tmpr = new_TorusPolynomial(N);
//     IntPolynomial_ifft(tmp+0,poly1);
//     TorusPolynomial_ifft(tmp+1,poly2);
//     LagrangeHalfCPolynomialMul(tmp+2,tmp+0,tmp+1);
//     TorusPolynomial_fft(tmpr, tmp+2);
//     torusPolynomialAddTo(result, tmpr);
//     delete_TorusPolynomial(tmpr);
//     delete_LagrangeHalfCPolynomial_array(3,tmp);
// }
//tGswExternMulToTLwe

  IntPolynomial* alloc_IntPolynomial_array(int32_t nbelts) {
    return (IntPolynomial*)  mymalloc(nbelts*sizeof(IntPolynomial));
}

// IntPolynomial::IntPolynomial(const int32_t N): N(N)
// {
//     this->coefs = new int32_t[N]; 
// }

//   void init_IntPolynomial_array(int32_t nbelts, IntPolynomial* obj, const int32_t N) {
//     for (int32_t i=0; i<nbelts; i++) {
// 	new(obj+i) IntPolynomial(N);
//     }
// }

void init_IntPolynomial_array(int32_t nbelts, struct IntPolynomial* obj, const int32_t N) {
    for (int32_t i = 0; i < nbelts; i++) {
        // 为每个元素分配内存
        obj[i].coefs = (int32_t*) mymalloc(N * sizeof(int32_t));

        // 初始化 N
        obj[i].N = N;

        // 在这里可以根据需要对 coefs 数组进行初始化
        // 例如，可以将 coefs 数组的元素全部初始化为0
        for (int32_t j = 0; j < N; j++) {
            obj[i].coefs[j] = 0;
        }
    }
}

  IntPolynomial* new_IntPolynomial_array(int32_t nbelts, const int32_t N) {
    IntPolynomial* obj = (IntPolynomial*)  mymalloc(nbelts*sizeof(IntPolynomial));

    init_IntPolynomial_array(nbelts,obj,N);
    return obj;
}



#define INCLUDE_TGSW_TORUS32POLYNOMIAL_DECOMP_H 


#if defined INCLUDE_ALL || defined INCLUDE_TGSW_TORUS32POLYNOMIAL_DECOMP_H
#undef INCLUDE_TGSW_TORUS32POLYNOMIAL_DECOMP_H
  void
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

  void tGswTLweDecompH(IntPolynomial *result, const TLweSample *sample, const TGswParams *params) {
    const int32_t k = params->tlwe_params->k;
    const int32_t l = params->l;

    for (int32_t i = 0; i <= k; ++i) // b=a[k]
        tGswTorus32PolynomialDecompH(result + (i * l), &sample->a[i], params);
}

//Arithmetic operations on TLwe samples
/** result = (0,0) */
  void tLweClear(TLweSample *result, const TLweParams *params) {
    const int32_t k = params->k;

    for (int32_t i = 0; i < k; ++i) torusPolynomialClear(&result->a[i]);
    torusPolynomialClear(result->b);
    result->current_variance = 0.;
}




// Norme Euclidienne d'un IntPolynomial
  double intPolynomialNormSq2(const IntPolynomial *poly) {
    const int32_t N = poly->N;
    int32_t temp1 = 0;

    for (int32_t i = 0; i < N; ++i) {
        int32_t temp0 = poly->coefs[i] * poly->coefs[i];
        temp1 += temp0;
    }
    return temp1;
}


void * torusPolynomialAddMulR(TorusPolynomial *result, const IntPolynomial *poly1,
                                                const TorusPolynomial *poly2) {
    extern int mult_count;
    const int32_t N = poly1->N;
    const int32_t *a = poly1->coefs;
    const Torus32 *b = poly2->coefsT;
    Torus32 *res = result->coefsT;

    for (int32_t i = 0; i < N; ++i) {
        res[i] += a[i] * b[i];
        mult_count++;
    }
    //return result;
}



  void
tLweAddMulRTo(TLweSample *result, const IntPolynomial *p, const TLweSample *sample, const TLweParams *params) {
    const int32_t k = params->k;

    for (int32_t i = 0; i <= k; ++i)
        torusPolynomialAddMulR(result->a + i, p, sample->a + i);
    result->current_variance += intPolynomialNormSq2(p) * sample->current_variance;
}

  void tGswExternMulToTLwe(TLweSample *accum, const TGswSample *sample, const TGswParams *params) {
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
  void torusPolynomialMulByXaiMinusOne(TorusPolynomial *result, int32_t a, const TorusPolynomial *source) {
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
  void tLweMulByXaiMinusOne(TLweSample *result, int32_t ai, const TLweSample *bk, const TLweParams *params) {
    const int32_t k = params->k;
    for (int32_t i = 0; i <= k; i++)
        torusPolynomialMulByXaiMinusOne(&result->a[i], ai, &bk->a[i]);
}

  void tLweAddTTo(TLweSample *result, const int32_t pos, const Torus32 x, const TLweParams *params) {
    result->a[pos].coefsT[0] += x;
}



  void tLweAddTo(TLweSample *result, const TLweSample *sample, const TLweParams *params) {
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
  void tLweCopy(TLweSample *result, const TLweSample *sample, const TLweParams *params) {
    const int32_t k = params->k;
    const int32_t N = params->N;

    for (int32_t i = 0; i <= k; ++i)
        for (int32_t j = 0; j < N; ++j)
            result->a[i].coefsT[j] = sample->a[i].coefsT[j];

    result->current_variance = sample->current_variance;
}



  void
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
  void torusPolynomialMulByXai(TorusPolynomial *result, int32_t a, const TorusPolynomial *source) {
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


  void torusPolynomialCopy(
        TorusPolynomial *result,
        const TorusPolynomial *sample) {
    //assert(result != sample);
    const int32_t N = result->N;
    const Torus32 *__restrict s = sample->coefsT;
    Torus32 *__restrict r = result->coefsT;

    for (int32_t i = 0; i < N; ++i) r[i] = s[i];
}

/** result = (0,mu) */
  void tLweNoiselessTrivial(TLweSample *result, const TorusPolynomial *mu, const TLweParams *params) {
    const int32_t k = params->k;

    for (int32_t i = 0; i < k; ++i) torusPolynomialClear(&result->a[i]);
    torusPolynomialCopy(result->b, mu);
    result->current_variance = 0.;
}


  void tLweExtractLweSampleIndex(LweSample* result, const TLweSample* x, const int32_t index, const LweParams* params,  const TLweParams* rparams) {
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


  void tLweExtractLweSample(LweSample* result, const TLweSample* x, const LweParams* params,  const TLweParams* rparams) {
    tLweExtractLweSampleIndex(result, x, 0, params, rparams);
}

//   TorusPolynomial* alloc_TorusPolynomial_array(int32_t nbelts) {
//     return (TorusPolynomial*)  mymalloc(nbelts*sizeof(TorusPolynomial));
// }

// TorusPolynomial::TorusPolynomial(const int32_t N): N(N)
// {
//     this->coefsT = new Torus32[N]; 
// }

//   void init_TorusPolynomial_array(int32_t nbelts, TorusPolynomial* obj, const int32_t N) {
//     for (int32_t i=0; i<nbelts; i++) {
// 	new(obj+i) TorusPolynomial(N);
//     }
// }

// void init_TorusPolynomial_array(int32_t nbelts, struct TorusPolynomial* obj, const int32_t N) {
//     for (int32_t i = 0; i < nbelts; i++) {
//         TorusPolynomial_init(&obj[i], N);
//     }
// }



//   TorusPolynomial* new_TorusPolynomial_array(int32_t nbelts, const int32_t N) {
//     TorusPolynomial* obj = alloc_TorusPolynomial_array(nbelts);
//     init_TorusPolynomial_array(nbelts,obj,N);
//     return obj;
// }



  void tfhe_blindRotateAndExtract(LweSample *result,
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



  void tfhe_bootstrap_woKS(LweSample *result,
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
    int32_t* bara = (int32_t*) mymalloc(N * sizeof(int32_t));
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
  void lweKeySwitch(LweSample* result, const LweKeySwitchKey* ks, const LweSample* sample){
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
    struct LweSample* sample = (struct LweSample*)  mymalloc(sizeof(struct LweSample));
    sample->a = (Torus32*) mymalloc(params->n * sizeof(Torus32));
    // 初始化其他成员
    sample->b = 0;
    sample->current_variance = 0.0;
    return sample;
}


  void tfhe_bootstrap(LweSample *result,
                           const LweBootstrappingKey *bk,
                           Torus32 mu, const LweSample *x) {

    struct LweSample *u = new_LweSample(&bk->accum_params->extracted_lweparams);
    
    // struct LweSample *u = new_gate_bootstrapping_ciphertext_array(1, bk->accum_params->extracted_lweparams);
    // printf("1");
    tfhe_bootstrap_woKS(u, bk, mu, x);
    // Key Switching
    // printf("2");
    lweKeySwitch(result, bk->ks, u);

    //delete_LweSample(u);
}

  Torus32 modSwitchToTorus32(int32_t mu, int32_t Msize){
    uint64_t interv = ((UINT64_C(1)<<63)/Msize)*2; // width of each intervall
    uint64_t phase64 = mu*interv;
    //floor to the nearest multiples of interv
    return phase64>>32;
}


// void LweSample_init(struct LweSample* sample, const struct LweParams* params) {
//     sample->a = (Torus32*) mymalloc(params->n * sizeof(Torus32));
//     // 初始化其他成员
//     sample->b = 0;
//     sample->current_variance = 0.0;
// }

void LweSample_init(struct LweSample* sample, const struct LweParams* params) {
    sample->a = (Torus32*) mymalloc(params->n * sizeof(Torus32));//4*630B
    myFree(sample->a);
    sample->a = (struct Torus32*) write_data_to_disk( sample->a,params->n * sizeof(Torus32) ,0);
    sample->b = 0;
    sample->current_variance = 0.0;
}

// 分配一个包含nbelts个LweSample的数组
struct LweSample* alloc_LweSample_array(int32_t nbelts) {
    return (struct LweSample*) mymalloc(nbelts * sizeof(struct LweSample));//16KB
}

// 初始化LweSample数组
void init_LweSample_array(int32_t nbelts, struct LweSample* obj, const struct LweParams* params) {
    for (int32_t i = 0; i < nbelts; i++) {
        LweSample_init(&obj[i], params);
    }
}

// 创建并初始化LweSample数组
struct LweSample* new_LweSample_array(int32_t nbelts, const struct LweParams* params) {
// #ifdef SUNFLOWER
    //printf("hi~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
    struct LweSample* obj = (struct LweSample*) malloc(nbelts * sizeof(struct LweSample));//16KB
    init_LweSample_array(nbelts, obj, params);
    for(int i=0;i<nbelts;i++){   
        //here have some bug needed to be fixed
        struct LweSample* temp_sample = (struct LweSample*)  mymalloc(sizeof(struct LweSample));
        LweSample_init(temp_sample, params);
        //myFree(temp_sample->a);

       // write_data_to_disk( temp_sample,sizeof(struct LweSample) );
        myFree(temp_sample);
    }
// #else
// printf("mmmmmmmmmmmmmmmmmmmmmmmmmmm\n");
//     struct LweSample* obj = curr_disk_address;
//     // init_LweSample_array(nbelts, obj, params);
//         for(int i=0;i<nbelts;i++){   
//         //here have some bug needed to be fixed
//         struct LweSample* temp_sample = (struct LweSample*)  mymalloc(sizeof(struct LweSample));
//         LweSample_init(temp_sample, params);
//         //write_data_to_disk( temp_sample,sizeof(struct LweSample) );
//         myFree(temp_sample);
//     }
// #endif
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

void TGswParams_init(struct TGswParams* params,int32_t l, int32_t Bgbit, const struct TLweParams* tlwe_params) {
    //struct TGswParams* params = (struct TGswParams*) mymalloc(sizeof(struct TGswParams));
    params->l = l;
    params->Bgbit = Bgbit;
    params->Bg = 1 << Bgbit;
    params->halfBg = params->Bg / 2;
    params->maskMod = params->Bg - 1;
    params->tlwe_params = tlwe_params;
    params->kpl = (tlwe_params->k + 1) * l;

    params->h = (Torus32*) mymalloc(l * sizeof(Torus32));
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
    //return params;
}


  void lweKeyGen(LweKey* result) {

int32_t alpha = result->params->alpha_min;
    int32_t n = result->params->n;
    
    // int32_t n = 1024;
  //uniform_int_distribution<int32_t> distribution(0,1);
  for (int32_t i=0; i<n; i++) 
    // result->key[i]=distribution(generator);
    result->key[i]=i;

}

struct TLweKey * new_tlwe_key(const struct TLweParams *params) {
    struct TLweKey *obj = (struct TLweKey *)  mymalloc(sizeof(struct TLweKey));
    obj->params = params;
    // obj->key = new_TorusPolynomial_array(params->k, params->N);
    obj->key = new_IntPolynomial_array(params->k, params->N);
    return obj;
}



// TLwe
  void tLweKeyGen(TLweKey *result) {
     int32_t N = result->params->N;
     int32_t k = result->params->k;
    //uniform_int_distribution<int32_t> distribution(0, 1);

    for (int32_t i = 0; i < k; ++i)
        for (int32_t j = 0; j < N; ++j)
            //result->key[i].coefs[j] = distribution(generator);
            result->key[i].coefs[j] = i+j;
}


  void tGswKeyGen(TGswKey *result) {
    tLweKeyGen(&result->tlwe_key);
}

/** generate a tgsw key (in fact, a tlwe key) */



struct LweKey* new_LweKey( const struct LweParams* params) {
    struct LweKey* obj = (struct LweKey*) mymalloc(sizeof(struct LweKey));
    //struct LweKey* obj ;

    //int32_t *key = (int32_t*) mymalloc(params->n * sizeof(int32_t));
    obj->params = params;
    obj->key = (int32_t*) mymalloc(params->n * sizeof(int32_t));
    myFree(obj->key);
    obj->key = (int32_t*) write_data_to_disk(obj->key,params->n * sizeof(int32_t),0);

    obj = (struct LweKey*) write_data_to_disk(obj,sizeof(struct LweKey),1);
    //struct LweKey obj = {params, key};
    return obj;
}

struct TGswKey* new_TGswKey( const struct TGswParams* params) {
    // struct TGswKey* obj;
    // obj->params = params;
    // obj->tlwe_params = params->tlwe_params;
    struct TGswKey* obj = (struct TGswKey*) mymalloc(sizeof(struct TGswKey));
    struct IntPolynomial *key = (struct IntPolynomial*) mymalloc(params->kpl * sizeof(struct IntPolynomial));
    // struct TGswKey obj = {params, params->tlwe_params, key,params->tlwe_params};
    obj->params = params;
    obj->tlwe_params = params->tlwe_params;
    obj->key = key;
    obj->tlwe_key = *new_tlwe_key(params->tlwe_params);//
    //obj->tlwe_key = obj->tlwe_params;
    // 初始化其他成员
    return obj;
}
struct TGswSample * new_TGswSample( const struct TGswParams *params) {
    struct TGswSample * obj = (struct TGswSample *) mymalloc(sizeof(struct TGswSample));
    const int32_t k = params->tlwe_params->k;
    const int32_t l = params->l;
    struct TLweSample *all_sample = new_TLweSample_array((k + 1) * l, params->tlwe_params);
    if(all_sample == NULL)
        printf("err!!!!  failed  mymalloc all_sample~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    // for (int i=0;i<k+1;i++){
    //     struct TLweSample * temp = all_sample + p * l; 
    // }
    struct TLweSample **bloc_sample = (struct TLweSample **) mymalloc((k + 1) * sizeof(struct TLweSample *));
    myFree(bloc_sample);
    bloc_sample =(struct TLweSample *) write_data_to_disk(bloc_sample,(k + 1) *sizeof(struct TLweSample *),0);

    printf("the value of bloc_sample is %x \n",bloc_sample);
    // Block * block = (Block *)((char *)bloc_sample - sizeof(Block));
    //block->size = (k + 1) * sizeof(struct TLweSample *);
    // printf("the value of block is %x \n",block);/
    // printf("the size of blocsample is %d \n",block->size);
    // myFree(bloc_sample);
    // struct TLweSample **bloc_sample = (struct TLweSample **) malloc((k + 1) * sizeof(struct TLweSample *));
    // bloc_sample = (struct TLweSample **) mymalloc((k + 1) * sizeof(struct TLweSample *));

    printf("the size of bloc_sample is %d \n",sizeof(bloc_sample));
    // printf("the malloc bsaple size is %d",)

    printf("the size of b_sample is %d \n",(k + 1) * sizeof(struct TLweSample *));
    for (int32_t p = 0; p < k + 1; ++p) {
        bloc_sample[p] = all_sample + p * l;
    }
    // struct TGswSample * obj = {params, bloc_sample, all_sample};
    // obj->params = params;
    obj->bloc_sample = bloc_sample;
    obj->all_sample = all_sample;
    obj->k = params->tlwe_params->k;
    obj->l = params->l;
    myFree(obj);
    obj = (struct TGswSample *) write_data_to_disk(obj,sizeof(struct TGswSample),1);

    //myFree(bloc_sample);
    // for (int i = 0; i < (k + 1) * l; i++) {
    //     myFree(&all_sample[i]);
    // }
    printf("hiiiiiiiiiiiii");
    //myFree(all_sample);
    //myFree(obj);
    return obj;
}



struct TGswSample *new_TGswSample_array(int32_t nbelts, const struct TGswParams *params) {
    // printf("the size of array is %d\n",nbelts);
    // #ifndef SUNFLOWER
    struct TGswSample *obj = (struct TGswSample *) mymalloc(nbelts * sizeof(struct TGswSample));//16KB
    myFree(obj);
    obj = (struct TGswSample *) write_data_to_disk(obj,nbelts * sizeof(struct TGswSample),0);
    // obj = (struct TGswSample *) malloc(nbelts * sizeof(struct TGswSample));//16KB
    
    // #else
    // struct TGswSample *obj = (struct TGswSample *) mymalloc(nbelts * sizeof(struct TGswSample));//16KB
    // #endif

    // printf("the size of array is %d\n",nbelts);
    if(obj == NULL)
        printf("err!!!!  failed  mymalloc tgswSample~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    else
        printf("success  mymalloc tgswSample array\n");
    for (int  i = 0; i < nbelts; i++) {
        printf("    i is %d\n",i);    
        struct TGswSample * temp = new_TGswSample(params);
        //myFree(temp);
        obj[i] = *temp;
        //write_data_to_disk(&obj[i],sizeof(struct TGswSample));
    }
    //myFree(obj);
    //obj = (struct TGswSample *) write_data_to_disk(obj,nbelts * sizeof(struct TGswSample),1);
    return obj;
}
// struct TGswSampleFFT *new_TGswSampleFFT_array(int32_t nbelts, const struct TGswParams *params) {
//     struct TGswSampleFFT *obj = (struct TGswSampleFFT *) mymalloc(nbelts * sizeof(struct TGswSampleFFT));
//     for (int32_t i = 0; i < nbelts; i++) {
//         obj[i] = *new_TGswSampleFFT(params);
//     }
//     return obj;
// }



void init_LweKeySwitchKey(struct LweKeySwitchKey *obj, int32_t n, int32_t t, int32_t basebit, const struct LweParams *out_params) {
    const int32_t base = 1 << basebit;
    // struct LweSample *ks0_raw = new_LweSample_array(n * t * base, out_params);
    //
    //struct LweSample *ks0_raw = (struct LweSample *) mymalloc(n * t * base * sizeof(struct LweSample));// here we must replace it with real address
    //myFree(ks0_raw);
    struct LweSample* ks0_raw = (struct LweSample *) write_data_to_disk(ks0_raw,n * t * base * sizeof(struct LweSample),0);
    // for (int32_t p = 0; p < n * t * base; ++p) {
    //     struct LweSample * temp = (struct LweSample *)  mymalloc(sizeof(struct LweSample));
    //     myFree(temp);
    //     temp = (struct LweSample *) write_data_to_disk(temp,sizeof(struct LweSample),0);
    //     struct LweSample temp_data = *temp;
    //     ks0_raw[p] = *temp;
    //     // struct LweSample * data_sd_addr = write_data_to_disk(&ks0_raw[p],sizeof(struct LweSample),1);
    //     // ks0_raw[p] = *data_sd_addr;
    // }
    // // struct LweSample *ks0_raw = new_LweSample_array(n * t * base, out_params);

    printf("2222222222222222222222222222222222222\n");
    obj->n = n;
    obj->t = t;
    obj->basebit = basebit;
    obj->out_params = out_params;
    obj->ks0_raw = ks0_raw;
    //obj->ks1_raw = (struct LweSample **)  malloc(n * t * sizeof(struct LweSample *));//4KB   //we need chere the size
    obj->ks1_raw = (struct LweSample **)  write_data_to_disk(obj->ks1_raw,n * t * sizeof(struct LweSample *),0);
    // obj->ks = (LweSample***) malloc(n * sizeof(LweSample**));//4KB
    obj->ks = (struct LweSample***) write_data_to_disk(obj->ks,n * sizeof(struct LweSample**),0);
    for (int32_t p = 0; p < n*t; ++p)
	    obj->ks1_raw[p] = obj->ks0_raw + base*p;  // direct modify the real addr
	for (int32_t p = 0; p < n; ++p)
	    obj->ks[p] = obj->ks1_raw + t*p;
}
  LweKeySwitchKey* alloc_LweKeySwitchKey() {
    return (LweKeySwitchKey*)mymalloc(sizeof(LweKeySwitchKey));
}
  LweKeySwitchKey* new_LweKeySwitchKey(int32_t n, int32_t t, int32_t basebit, const LweParams* out_params) {
    LweKeySwitchKey* obj = alloc_LweKeySwitchKey();
    myFree(obj);
    // obj = (LweKeySwitchKey*)  malloc(sizeof(LweKeySwitchKey));
    obj = (struct LweKeySwitchKey *)write_data_to_disk(obj,sizeof(LweKeySwitchKey),0);
    // printf("success  mymalloc lwekey\n");
    init_LweKeySwitchKey(obj, n,t,basebit,out_params);

    // if(obj->out_params == NULL)
    //     // printf("out_params is null\n");
    //     else
            // printf("out_params is not null\n");
    return obj;
}


  void init_LweBootstrappingKey(LweBootstrappingKey *obj, int32_t ks_t, int32_t ks_basebit, const LweParams *in_out_params,
                                     const TGswParams *bk_params) {
     TLweParams *accum_params = bk_params->tlwe_params;
     LweParams *extract_params = &accum_params->extracted_lweparams;
     int32_t n = in_out_params->n;
     int32_t N = extract_params->n;
    // printf("the value of n is %d\n",n);//
    //here we need to record the sdram address of the bk and ks
    TGswSample *bk = new_TGswSample_array(n, bk_params);
    printf("111");
    LweKeySwitchKey *ks = new_LweKeySwitchKey(N, ks_t, ks_basebit, in_out_params);
    // printf("222");
    //LweBootstrappingKey *obj = {in_out_params, bk_params, accum_params, extract_params, bk, ks};
    //new(obj) LweBootstrappingKey(in_out_params, bk_params, accum_params, extract_params, bk, ks);
    obj->in_out_params = in_out_params;
    obj->bk_params = bk_params;
    obj->accum_params = accum_params;
    obj->extract_params = extract_params;
    obj->bk = bk;
    obj->ks = ks;
}

  LweBootstrappingKey *alloc_LweBootstrappingKey() {
    return (LweBootstrappingKey *)  mymalloc(sizeof(LweBootstrappingKey));
}

  LweBootstrappingKey *
new_LweBootstrappingKey(const int32_t ks_t, const int32_t ks_basebit, const LweParams *in_out_params,
                        const TGswParams *bk_params) {
    LweBootstrappingKey *obj = (LweBootstrappingKey *)  mymalloc(sizeof(LweBootstrappingKey));
   // myFree(obj);
    //obj = (LweBootstrappingKey *) write_data_to_disk(obj,sizeof(LweBootstrappingKey),0);
    init_LweBootstrappingKey(obj, ks_t, ks_basebit, in_out_params, bk_params);
    printf("finish new_LweBootstrappingKey\n");
    
    // if(obj->accum_params == NULL)
    //     printf("accum_params is null\n");
    //     else
    //         printf("accum_params is not null\n");
    // if(obj->extract_params == NULL)
    //     printf("extract_params is null\n");
    //     else
    //         printf("extract_params is not null\n");
    return obj;
}

  void tLweExtractKey(LweKey* result, const TLweKey* key) //sans doute un param supplémentaire
{
    const int32_t N = key->params->N;
    const int32_t k = key->params->k;
    //assert(result->params->n == k*N);
    for (int32_t i=0; i<k; i++) {
	for (int32_t j=0; j<N; j++)
	    result->key[i*N+j]=key->key[i].coefs[j];
    }
}


//   void lweNoiselessTrivial(LweSample* result, Torus32 mu, const LweParams* params){
//     const int32_t n = params->n;

//     for (int32_t i = 0; i < n; ++i) result->a[i] = 0;
//     result->b = mu;
//     result->current_variance = 0.;
// }

typedef int int32_t;
  Torus32 dtot32(double d) {
    return (Torus32)(int64_t)((d - (double)(int64_t)d) * (double)_two32);
    }
  void lweSymEncryptWithExternalNoise(LweSample* result, Torus32 message, double noise, double alpha, const LweKey* key){
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


  void lweCreateKeySwitchKey(LweKeySwitchKey* result, const LweKey* in_key, const LweKey* out_key){
    const int32_t n = result->n;
    const int32_t t = result->t;
    const int32_t basebit = result->basebit;
    const int32_t base = 1<<basebit;
    const double alpha = out_key->params->alpha_min;
    const int32_t sizeks = n*t*(base-1);
    //const int32_t n_out = out_key->params->n;
    // printf("2.1\n");
    double err = 0;
    // chose a random vector of gaussian noises   // note here the malloc space is too huge
    // double* noise = (double*) malloc(sizeof(double) * sizeks);
    double * noise = (double *) write_data_to_disk(noise,sizeof(double) * sizeks,0);


    for (int32_t i = 0; i < sizeks; ++i){
        //normal_distribution<double> distribution(0.,alpha); 
        noise[i] = (double)i / (double)sizeks;
        err += noise[i];
    }
    // recenter the noises
    // printf("2.2\n");
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

  void torusPolynomialUniform(TorusPolynomial *result) {
    const int32_t N = result->N;
    Torus32 *x = result->coefsT;

    for (int32_t i = 0; i < N; ++i)
        // x[i] = uniformTorus32_distrib(generator);
        x[i] = i;
}



//   void torusPolynomialAddMulRFFT(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2) {
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

  void tLweSymEncryptZero(TLweSample *result, double alpha, const TLweKey *key) {
    // if(key->params == NULL)
    //     printf("key->params is NULL\n");
    //     else
    //         printf("key->params is not NULL\n");
    // if(key->key == NULL)
    //     printf("key->key is NULL\n");
    //     else
    //         printf("key->key is not NULL\n");

    // if(result->a == NULL)
    //     printf("result->a is NULL\n");
    //     else
    //         printf("result->a is not NULL\n");
    // if(result->b == NULL)
    //     printf("result->b is NULL\n");
    //     else
    //         printf("result->b is not NULL\n");
    const int32_t N = key->params->N;
    const int32_t k = key->params->k;

    for (int32_t j = 0; j < N; ++j){
        // result->b->coefsT[j] = gaussian32(0, alpha);
        //printf("tlwesymencrypt\n");}
        result->b->coefsT[j] = j;}
    for (int32_t i = 0; i < k; ++i) {
        torusPolynomialUniform(&result->a[i]);
        torusPolynomialAddMulR(result->b, &key->key[i], &result->a[i]);
    }

    result->current_variance = alpha * alpha;
}

  void tGswEncryptZero(TGswSample *result, double alpha, const TGswKey *key) {
    const TLweKey *rlkey = &key->tlwe_key;
    const int32_t kpl = key->params->kpl;

    for (int32_t p = 0; p < kpl; ++p) {
        tLweSymEncryptZero(&result->all_sample[p], alpha, rlkey);
    }
}

  void tGswAddMuIntH(TGswSample *result, const int32_t message, const TGswParams *params) {
    const int32_t k = params->tlwe_params->k;
    const int32_t l = params->l;
    const Torus32 *h = params->h;

    // compute result += H
    for (int32_t bloc = 0; bloc <= k; ++bloc)
        for (int32_t i = 0; i < l; i++)
            result->bloc_sample[bloc][i].a[bloc].coefsT[0] += message * h[i];
}

  void tGswSymEncryptInt(TGswSample *result, const int32_t message, double alpha, const TGswKey *key) {
    tGswEncryptZero(result, alpha, key);
    tGswAddMuIntH(result, message, key->params);
}

  void tfhe_createLweBootstrappingKey(
        LweBootstrappingKey *bk,
        const LweKey *key_in,
        const TGswKey *rgsw_key) {
    // assert(bk->bk_params == rgsw_key->params);
    // assert(bk->in_out_params == key_in->params);
    const LweParams *in_out_params = (const LweParams *)mymalloc(sizeof(LweParams));
    const TGswParams *bk_params = (const TGswParams *)mymalloc(sizeof(TGswParams));
    const TLweParams *accum_params = (const TLweParams *)mymalloc(sizeof(TLweParams));
    const LweParams *extract_params = (const LweParams *)mymalloc(sizeof(LweParams));
    const TLweKey *accum_key = (const TLweKey *)mymalloc(sizeof(TLweKey));
    in_out_params = bk->in_out_params;
    bk_params = bk->bk_params;
    accum_params = bk_params->tlwe_params;
    extract_params = &accum_params->extracted_lweparams;

    //LweKeySwitchKey* ks; ///< the keyswitch key (s'->s)
    accum_key = &rgsw_key->tlwe_key;

    LweKey *extracted_key = new_LweKey(extract_params);//4KB
    // printf("1\n");
    tLweExtractKey(extracted_key, accum_key);
    // printf("2\n");
    // if(bk->ks == NULL)
    //     printf("bk->ks is null\n");
    //     else printf("bk->ks is not null\n");
    lweCreateKeySwitchKey(bk->ks, extracted_key, key_in);
    // printf("3\n");
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
    // printf("4\n");


}
// LweBootstrappingKeyFFT* new_LweBootstrappingKeyFFT(const struct TFheGateBootstrappingParameterSet *params){
//     LweBootstrappingKeyFFT* ptr = (LweBootstrappingKeyFFT*)  mymalloc(sizeof(LweBootstrappingKeyFFT));
//     ptr->accum_params = params->tgsw_params->tlwe_params;
//     ptr->extract_params = &params->tgsw_params->tlwe_params->extracted_lweparams;
//     ptr->bk_params = params->tgsw_params;
//     ptr->in_out_params = params->in_out_params;
//     ptr->ks = new_LweKeySwitchKey(ptr->extract_params->n, params->ks_t ,params->ks_basebit, ptr->in_out_params);
//     ptr->bkFFT = new_TGswSampleFFT_array(ptr->in_out_params->n, ptr->bk_params);


//     return ptr;


// }

  TFheGateBootstrappingSecretKeySet *
new_random_gate_bootstrapping_secret_keyset(const TFheGateBootstrappingParameterSet *params) {
    TFheGateBootstrappingSecretKeySet *obj = (TFheGateBootstrappingSecretKeySet *)  mymalloc(
            sizeof(TFheGateBootstrappingSecretKeySet));

    // printf("the value of params->in_out_params->n is : %d\n",params->in_out_params->n);
    printf("generate the new_LweBootstrappingKey");
    LweBootstrappingKey *bk = new_LweBootstrappingKey(params->ks_t, params->ks_basebit, params->in_out_params,
                                                      params->tgsw_params);
    printf("finish the new_LweBootstrappingKey");
    // if(bk == NULL)

    //     printf("bk is null\n");
    //     else printf(" mymalloc the bootstrap key successfully\n");
    printf("the address of params->in_out_params is : %p\n",params->in_out_params);
    LweKey *lwe_key = new_LweKey(params->in_out_params);
    lweKeyGen(lwe_key);
    printf("the address of lwe_key is : %p\n",lwe_key);

    TGswKey *tgsw_key = new_TGswKey(params->tgsw_params);
    printf("the address of tgsw_key is : %p\n",tgsw_key);
    tGswKeyGen(tgsw_key);
    tfhe_createLweBootstrappingKey(bk, lwe_key, tgsw_key);
    myFree(lwe_key);
    myFree(tgsw_key);

    printf("init the bk successfully\n");
    // LweBootstrappingKeyFFT *bkFFT  = new_LweBootstrappingKeyFFT(params);
    // = new_LweBootstrappingKeyFFT(bk);
    //return new TFheGateBootstrappingSecretKeySet(params, bk, bkFFT, lwe_key, tgsw_key);
    // printf("4\n");
    obj->params = params;
    obj->lwe_key = lwe_key; 
    obj->tgsw_key = tgsw_key;   
    obj->cloud = *new_TFheGateBootstrappingCloudKeySet(params,bk);

    //TFheGateBootstrappingSecretKeySet  key = {params, bk, bkFFT, lwe_key, tgsw_key};
    return obj;
}

void init_parms(struct TFheGateBootstrappingParameterSet * params){

    static const int32_t N = N_value;//the prime value is 1024
    static const int32_t k = 1;
    static const int32_t n = n_value;
    static const int32_t bk_l = 2;
    static const int32_t bk_Bgbit = 7;
    static const int32_t ks_basebit = 2;
    static const int32_t ks_length = 8;
    //static const double ks_stdev = pow(2.,-15); //standard deviation
    // static const double bk_stdev = pow(2.,-25);; //standard deviation
    static const double ks_stdev = 1.0 / (1 << 15); //standard deviation
    static const double bk_stdev = 1.0 / (1 << 25);; //standard deviation
    static const double max_stdev = 0.012467; //max standard deviation for a 1/4 msg space
    //LweParams params_in_temp = {n, ks_stdev, max_stdev};
    LweParams *params_in = (LweParams *) mymalloc(sizeof(LweParams));
    printf("the address of params_in is : %p\n",params_in);
    params_in->n = n;
    params_in->alpha_min = ks_stdev;
    params_in->alpha_max = max_stdev;
    //LweParams extracted_lweparams = {N * k, bk_stdev, max_stdev};
    // TLweParams params_accum_temp = {N, k, bk_stdev, max_stdev,extracted_lweparams};
    // TLweParams *params_accum = &params_accum_temp;
    TLweParams *params_accum = (TLweParams *) mymalloc(sizeof(TLweParams));
    params_accum->N = N;
    params_accum->k = k;
    params_accum->alpha_min = bk_stdev;
    params_accum->alpha_max = max_stdev;
    params_accum->extracted_lweparams.alpha_max = max_stdev;
    params_accum->extracted_lweparams.alpha_min = bk_stdev;
    params_accum->extracted_lweparams.n = N*k;

    //TGswParams params_bk_temp = {bk_l, bk_Bgbit, params_accum};
    TGswParams *params_bk = (TGswParams *) mymalloc(sizeof(TGswParams));
    TGswParams_init(params_bk,bk_l, bk_Bgbit, params_accum);
    // printf("the params_bk->tlwe_params->k is : %d\n",params_bk->tlwe_params->k);
    //TFheGateBootstrappingParameterSet params_temp = {ks_length, ks_basebit, params_in, params_bk};
    params->in_out_params = params_in;
    params->tgsw_params = params_bk;
    params->ks_t = ks_length;
    params->ks_basebit = ks_basebit;
    // myyFree(params_accum);
    // myyFree(params_bk);
    // myyFree(params_in);

}

struct Record {
    uint32_t address;
    uint32_t dataSize;
    char data[];
};



//use this function to simulate the operation transfer data from dataram to sdram
void* write_data_to_disk(void *data, size_t dataSize,uint32_t initialed_value) {
    // extern uint32_t curr_disk_address;
    if(initialed_value == 1){
        printf("the address of curr_disk_address is : %p\n",curr_disk_address);
        FILE *file = fopen(filename, "ab");
        if (file == NULL) {
            printf("Error opening file.\n");
            return;
        }
    // fwrite(&address, sizeof(uint32_t), 1, file);
        //fwrite(&dataSize, sizeof(size_t), 1, file);
        fwrite(data, dataSize, 1, file);
        fwrite(&curr_disk_address, sizeof(uint32_t), 1, file);
        fwrite(&dataSize, sizeof(size_t), 1, file);
        char newline = '\n';
        fwrite(&newline, sizeof(char), 1, file);
        fclose(file);
    }
    curr_disk_address += dataSize+sizeof(uint32_t)+sizeof(size_t)+sizeof(char);

    return curr_disk_address;
}

void* read_data_from_disk(size_t* dataSize) {
    FILE* file = fopen(filename, "rb");

    if (file == NULL) {
        printf("Error opening file.\n");
        return NULL;
    }

    fseek(file, curr_disk_address, SEEK_SET);
    fread(&curr_disk_address, sizeof(uint32_t), 1, file);

    fread(dataSize, sizeof(size_t), 1, file);
    void* data = malloc(*dataSize);

    fread(data, *dataSize, 1, file);

    fclose(file);
    return data;
}





int main()
{   //we need note the  current ptr
    //struct TFheGateBootstrappingParameterSet *params;
    printf("executing the main function~~\n");
    // FILE * prams_file = fopen("size_prams","wb");
    //params = new_tfheGateBootstrappingParameterSet_fromFile(prams_file);
    // int std_params = 1;
    // char memory[MEMORY_SIZE];
    printf("the start of memory is : %p   ,and the end of memory is :%p  \n",memory,memory + sizeof(memory));
    // char sdram[4*1024*1024];


    initMemoryPool(memory, sizeof(memory));
    init_sdram(sdram, sizeof(sdram));
    TFheGateBootstrappingParameterSet *params = (TFheGateBootstrappingParameterSet *) mymalloc(sizeof(TFheGateBootstrappingParameterSet));
    init_parms(params);
    TFheGateBootstrappingParameterSet * params_SD = (struct TFheGateBootstrappingParameterSet *)write_data_to_disk(params, sizeof(TFheGateBootstrappingParameterSet),1);
    //printf("the size of params is : %d\n",sizeof(TFheGateBootstrappingParameterSet));
    //printf("the address of params is :%x\n",params);
    int N = params->tgsw_params->tlwe_params->N;
    int n = params->in_out_params->n;
    double alpha = params->in_out_params->alpha_min;
    printf("the test of alpha is : %f\n",alpha);
    printf("the N is : %d , and the n is : %d\n ",N,n);
    printf("the product of N * n = : %d\n",N * n);
    double test = params->tgsw_params->tlwe_params->alpha_max;
    printf("the value of test is : %f\n",test);
    //////
    // }
   //printf("hello world\n");
  //TFheGateBootstrappingSecretKeySet *key = new_random_gate_bootstrapping_secret_keyset(params);
    //FILE * x7_file = fopen("size_bootstrap_secretkey","wb");
    //struct TFheGateBootstrappingSecretKeySet *key = new_tfheGateBootstrappingSecretKeySet_fromFile(x7_file);

    uint32_t seed[] = {214, 1592, 657};
    //tfhe_random_generator_setSeed(seed, 3);
    printf("the value of in_out_params->n is : %d\n",params->in_out_params->n);
    printf("generate the thfeGateBootstrappingSecretKeySet\n");
    TFheGateBootstrappingSecretKeySet *key = new_random_gate_bootstrapping_secret_keyset(params);
    printf("executring the bootstrap function,successfully generate key \n");
    const struct TFheGateBootstrappingCloudKeySet *bk = &(key->cloud); 

    struct LweSample *temp_result = new_LweSample_array(1,params->in_out_params);

    const Torus32 MU = modSwitchToTorus32(1, 8);
    //FILE* temp_result_file = fopen("temp_result","wb");
    //import_lweSample_fromFile(temp_result_file,temp_result,params->in_out_params);
    //here we must need know what operation we want to do
    struct LweSample *result = new_LweSample_array(1,params->in_out_params);
    printf("executing the bootstrap function\n");
    tfhe_bootstrap(result,bk->bk,MU,temp_result);
    printf("bootsrapping is done !! \n");

    printf("the mul conut is %d\n",mult_count);
    printf("the rate of mul is %f\n\n",mult_count / (N * n * 1.0));
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