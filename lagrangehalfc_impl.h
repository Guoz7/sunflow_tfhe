#ifndef LAGRANGEHALFC_IMPL_H
#define LAGRANGEHALFC_IMPL_H
typedef int32_t Torus32; //avant uint32_t
typedef unsigned int		uint32_t;
// #include <cassert>
//#include <cmath>
//#include <ccomplex>
//typedef double _Complex cplx;
// typedef std::complex<double> cplx;






#define INT64_C(c)  c ## LL

struct cplx {
    double real; // 实部
    double imag; // 虚部
};


//#include <fftw3.h>
//#include "tfhe.h"
//#include "polynomials.h"



// 定义FFT_Processor_fftw结构体
struct FFT_Processor_fftw {
    const int32_t _2N;
    const int32_t N;    
    const int32_t Ns2;
    double* rev_in;
    struct cplx* rev_out;
    struct cplx* in;
    double* out;
    void (*plan_fftw)();
    struct cplx* omegaxminus1;
    void (*execute_reverse_int)(struct cplx* res, const int32_t* a);
    void (*execute_reverse_torus32)(struct cplx* res, const Torus32* a);
    void (*execute_direct_Torus32)(Torus32* res, const struct cplx* a);
};

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

//     ~FFT_Processor_fftw();
// };

//extern  struct FFT_Processor_fftw fp1024_fftw;

/**
 * structure that represents a real polynomial P mod X^N+1
 * as the N/2 complex numbers:
 * P(w), P(w^3), ..., P(w^(N-1))
 * where w is exp(i.pi/N)
 */
struct LagrangeHalfCPolynomial_IMPL
{
   struct cplx* coefsC;
   struct FFT_Processor_fftw* proc;

//    LagrangeHalfCPolynomial_IMPL(int32_t N);
//    ~LagrangeHalfCPolynomial_IMPL();
};

#endif // LAGRANGEHALFC_IMPL_H



void plan_fftw(struct FFT_Processor_fftw* processor) {
    // Ensure FFTW plan thread safety
    static int initialized = 0;
    //static fftw_plan rev_p;
    //static fftw_plan p;
    // static fftw_mutex_t mutex;

    // fftw_execute_with_flops(1); // ensure fftw is initialized (dummy fftw call)

    // if (!initialized) {
    //     fftw_mutex_init(&mutex);
    //     initialized = 1;
    // }

    // fftw_mutex_lock(&mutex);
    // if (!rev_p) {
    //     processor->rev_p = fftw_plan_dft_r2c_1d(processor->_2N, processor->rev_in, processor->rev_out, FFTW_ESTIMATE);
    // }
    // if (!p) {
    //     processor->p = fftw_plan_dft_c2r_1d(processor->_2N, processor->in, processor->out, FFTW_ESTIMATE);
    // }
    // fftw_mutex_unlock(&mutex);
}


// Function to perform reverse FFT for int32_t array
void execute_reverse_int(struct cplx* res, const int32_t* a, const struct FFT_Processor_fftw* processor) {
    struct cplx* rev_out_cplx = processor->rev_out;
    double* rev_in = processor->rev_in;
    int32_t N = processor->N;

    for (int32_t i = 0; i < N; i++) {
        rev_in[i] = a[i] / 2.0;
    }

    for (int32_t i = 0; i < N; i++) {
        rev_in[N + i] = -rev_in[i];
    }

    processor->plan_fftw();

    for (int32_t i = 0; i < processor->Ns2; i++) {
        res[i] = rev_out_cplx[2 * i + 1];
    }

    for (int32_t i = 0; i <= processor->Ns2; i++) {
        //assert(fabs(rev_out_cplx[2 * i].real) < 1e-20);
        //assert(fabs(rev_out_cplx[2 * i].imag) < 1e-20);
    }
}

// Function to perform reverse FFT for Torus32 array
void execute_reverse_torus32(struct cplx* res, const Torus32* a, const struct FFT_Processor_fftw* processor) {
    static const double _2pm33 = 1.0 / (double)(INT64_C(1) << 33);
    int32_t* aa = (int32_t*)a;
    struct cplx* rev_out_cplx = processor->rev_out;
    double* rev_in = processor->rev_in;
    int32_t N = processor->N;

    for (int32_t i = 0; i < N; i++) {
        rev_in[i] = aa[i] * _2pm33;
    }

    for (int32_t i = 0; i < N; i++) {
        rev_in[N + i] = -rev_in[i];
    }

    processor->plan_fftw();

    for (int32_t i = 0; i < processor->Ns2; i++) {
        res[i] = rev_out_cplx[2 * i + 1];
    }

    for (int32_t i = 0; i <= processor->Ns2; i++) {
        //assert(fabs(rev_out_cplx[2 * i].real) < 1e-20);
        //assert(fabs(rev_out_cplx[2 * i].imag) < 1e-20);
    }
}

// Function to perform direct FFT for Torus32 array
void execute_direct_Torus32(Torus32* res, const struct cplx* a, const struct FFT_Processor_fftw* processor) {
    static const double _2p32 = (double)(INT64_C(1) << 32);
    const double _1sN = (double)(1) / (double)(processor->N);
    struct cplx* in_cplx = processor->in;
    double* out = processor->out;
    int32_t N = processor->N;
    int32_t Ns2 = processor->Ns2;

    for (int32_t i = 0; i <= Ns2; i++) {
        in_cplx[2 * i].real = 0.0;
        in_cplx[2 * i].imag = 0.0;
    }

    for (int32_t i = 0; i < Ns2; i++) {
        in_cplx[2 * i + 1] = a[i];
    }

    processor->plan_fftw();

    for (int32_t i = 0; i < N; i++) {
        res[i] = (Torus32)((int64_t)(out[i] * _1sN * _2p32));
    }

    for (int32_t i = 0; i < N; i++) {
        //assert(fabs(out[N + i] + out[i]) < 1e-20);
    }
}