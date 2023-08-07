#ifndef LAGRANGEHALFC_IMPL_H
#define LAGRANGEHALFC_IMPL_H
typedef int32_t Torus32; //avant uint32_t
//#include <cassert>
//#include <cmath>
//#include <ccomplex>
//typedef double _Complex cplx;
// typedef std::complex<double> cplx;


struct cplx {
    double real; // 实部
    double imag; // 虚部
};


#include <fftw3.h>
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
    void (*execute_reverse_int)(struct cplx* res, const int32_t* a);
    void (*execute_reverse_torus32)(struct cplx* res, const Torus32* a);
    void (*execute_direct_Torus32)(Torus32* res, const struct cplx* a);
//     ~FFT_Processor_fftw();
// };

extern  struct FFT_Processor_fftw fp1024_fftw;

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



void execute_reverse_int(struct cplx* res, const int* a) {
    cplx* rev_out_cplx = (cplx*) rev_out; //fftw_complex and cplx are layout-compatible
    for (int32_t i=0; i<N; i++) rev_in[i]=a[i]/2.;
    for (int32_t i=0; i<N; i++) rev_in[N+i]=-rev_in[i];
    fftw_execute(rev_p);
    for (int32_t i=0; i<Ns2; i++) res[i]=rev_out_cplx[2*i+1];
    for (int32_t i=0; i<=Ns2; i++) assert(abs(rev_out_cplx[2*i])<1e-20);
}
void execute_reverse_torus32(struct cplx* res, const Torus32* a) {
    static const double _2pm33 = 1./double(INT64_C(1)<<33);
    int32_t* aa = (int32_t*) a;
    cplx* rev_out_cplx = (cplx*) rev_out; //fftw_complex and cplx are layout-compatible
    for (int32_t i=0; i<N; i++) rev_in[i]=aa[i]*_2pm33;
    for (int32_t i=0; i<N; i++) rev_in[N+i]=-rev_in[i];
    fftw_execute(rev_p);
    for (int32_t i=0; i<Ns2; i++) res[i]=rev_out_cplx[2*i+1];
    for (int32_t i=0; i<=Ns2; i++) assert(abs(rev_out_cplx[2*i])<1e-20);
}
void execute_direct_Torus32(Torus32* res, const struct cplx* a) {
    static const double _2p32 = double(INT64_C(1)<<32);
    static const double _1sN = double(1)/double(N);
    cplx* in_cplx = (cplx*) in; //fftw_complex and cplx are layout-compatible
    for (int32_t i=0; i<=Ns2; i++) in_cplx[2*i]=0;
    for (int32_t i=0; i<Ns2; i++) in_cplx[2*i+1]=a[i];
    fftw_execute(p);
    for (int32_t i=0; i<N; i++) res[i]=Torus32(int64_t(out[i]*_1sN*_2p32));
    //pas besoin du fmod... Torus32(int64_t(fmod(rev_out[i]*_1sN,1.)*_2p32));
    for (int32_t i=0; i<N; i++) assert(fabs(out[N+i]+out[i])<1e-20);
}