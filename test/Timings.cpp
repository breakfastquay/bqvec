/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

#include "bqvec/VectorOpsComplex.h"
#include "bqvec/Allocators.h"

#include <iostream>
#include <cstdlib>

#include <time.h>

using namespace std;
using namespace breakfastquay;

#ifdef _WIN32
#define drand48() (-1+2*((float)rand())/RAND_MAX)
#endif

bool
testMultiply()
{
    cerr << "testVectorOps: testing v_multiply complex" << endl;

    const int iterations = 50000;
    const int N = 1024;

    bq_complex_t *target = allocate<bq_complex_t>(N);
    bq_complex_t *src1 = allocate<bq_complex_t>(N);
    bq_complex_t *src2 = allocate<bq_complex_t>(N);

    for (int i = 0; i < N; ++i) {
	src1[i].re = drand48();
	src1[i].im = drand48();
	src2[i].re = drand48();
	src2[i].im = drand48();
    }

    float divisor = float(CLOCKS_PER_SEC) / 1000.f;
    clock_t start = clock();

    double first, last, total = 0;
    for (int j = 0; j < iterations; ++j) {
        for (int i = 0; i < N; ++i) {
            bq_complex_t result;
            c_multiply(result, src1[i], src2[i]);
            if (i == 0) first = result.re;
            if (i == N-1) last = result.im;
            total += result.re;
            total += result.im;
            total += j;
        }
    }

    clock_t end = clock();

    cerr << "repeated c_multiply: first = " << first << ", last = " << last
         << ", total = " << total << endl;
    cerr << "time for repeated c_multiply: " << float(end - start)/divisor
         << endl;

    start = clock();

    first = last = total = 0;

    for (int j = 0; j < iterations; ++j) {
        v_multiply_to(target, src1, src2, N);
        for (int i = 0; i < N; ++i) {
            if (i == 0) first = target[i].re;
            if (i == N-1) last = target[i].im;
            total += target[i].re;
            total += target[i].im;
            total += j;
        }
    }

    end = clock();

    cerr << "v_multiply_to: first = " << first << ", last = " << last
         << ", total = " << total << endl;
    cerr << "time for v_multiply_to: " << float(end - start)/divisor << endl;

    deallocate(target);
    deallocate(src1);
    deallocate(src2);
    
    cerr << endl;
    return true;
}

bool
testPolarToCart()
{
    cerr << "testVectorOps: testing v_polar_to_cartesian" << endl;

    const int iterations = 50000;
    const int N = 1024;
    
    bq_complex_t *target = allocate<bq_complex_t>(N);
    bq_complex_element_t *mag = allocate<bq_complex_element_t>(N);
    bq_complex_element_t *phase = allocate<bq_complex_element_t>(N);

    for (int i = 0; i < N; ++i) {
	mag[i] = drand48();
	phase[i] = (drand48() * M_PI * 2) - M_PI;
    }

    float divisor = float(CLOCKS_PER_SEC) / 1000.f;
    clock_t start = clock();

    double first, last, total = 0;
    for (int j = 0; j < iterations; ++j) {
        for (int i = 0; i < N; ++i) {
            double real = mag[i] * cos(phase[i]);
            double imag = mag[i] * sin(phase[i]);
            if (i == 0) first = real;
            if (i == N-1) last = imag;
            total += real;
            total += imag;
            total += j;
        }
    }

    clock_t end = clock();

    cerr << "naive method: first = " << first << ", last = " << last
         << ", total = " << total << endl;
    cerr << "time for naive method: " << float(end - start)/divisor
         << endl;
    
    start = clock();

    first = last = total = 0;
    
    for (int j = 0; j < iterations; ++j) {
        v_polar_to_cartesian(target, mag, phase, N);
        for (int i = 0; i < N; ++i) {
            if (i == 0) first = target[i].re;
            if (i == N-1) last = target[i].im;
            total += target[i].re;
            total += target[i].im;
            total += j;
        }
    }

    end = clock();

    cerr << "v_polar_to_cartesian: first = " << first << ", last = " << last
         << ", total = " << total << endl;
    cerr << "time for v_polar_to_cartesian: " << float(end - start)/divisor
         << endl;
    
    deallocate(target);
    deallocate(mag);
    deallocate(phase);
    
    cerr << endl;
    return true;
}

bool
testPolarToCartInterleaved()
{
    cerr << "testVectorOps: testing v_polar_interleaved_to_cartesian" << endl;

    const int iterations = 50000;
    const int N = 1024;

    bq_complex_t *target = allocate<bq_complex_t>(N);
    bq_complex_element_t *source = allocate<bq_complex_element_t>(N*2);

    for (int i = 0; i < N; ++i) {
	source[i*2] = drand48();
	source[i*2+1] = (drand48() * M_PI * 2) - M_PI;
    }

    float divisor = float(CLOCKS_PER_SEC) / 1000.f;
    clock_t start = clock();

    double first, last, total = 0;
    for (int j = 0; j < iterations; ++j) {
        for (int i = 0; i < N; ++i) {
            double real = source[i*2] * cos(source[i*2+1]);
            double imag = source[i*2] * sin(source[i*2+1]);
            if (i == 0) first = real;
            if (i == N-1) last = imag;
            total += real;
            total += imag;
            total += j;
        }
    }

    clock_t end = clock();

    cerr << "naive method: first = " << first << ", last = " << last
         << ", total = " << total << endl;
    cerr << "time for naive method: " << float(end - start)/divisor
         << endl;
    
    start = clock();

    first = last = total = 0;
    
    for (int j = 0; j < iterations; ++j) {
        v_polar_interleaved_to_cartesian(target, source, N);
        for (int i = 0; i < N; ++i) {
            if (i == 0) first = target[i].re;
            if (i == N-1) last = target[i].im;
            total += target[i].re;
            total += target[i].im;
            total += j;
        }
    }

    end = clock();

    cerr << "v_polar_interleaved_to_cartesian: first = " << first << ", last = " << last
         << ", total = " << total << endl;
    cerr << "time for v_polar_interleaved_to_cartesian: " << float(end - start)/divisor
         << endl;
    
    deallocate(target);
    deallocate(source);
    
    cerr << endl;
    return true;
}

bool
testCartToPolarComplexTypes()
{
    cerr << "testVectorOps: testing v_cartesian_to_polar [complex types]" << endl;

    const int iterations = 50000;
    const int N = 1024;

    bq_complex_t *source = allocate<bq_complex_t>(N);
    bq_complex_element_t *mag = allocate<bq_complex_element_t>(N);
    bq_complex_element_t *phase = allocate<bq_complex_element_t>(N);

    for (int i = 0; i < N; ++i) {
        source[i].re = (drand48() * 2.0) - 1.0;
        source[i].im = (drand48() * 2.0) - 1.0;
    }

    float divisor = float(CLOCKS_PER_SEC) / 1000.f;
    clock_t start = clock();

    double first, last, total = 0;
    for (int j = 0; j < iterations; ++j) {
        for (int i = 0; i < N; ++i) {
            double mag = sqrt(source[i].re * source[i].re + source[i].im * source[i].im);
            double phase = atan2(source[i].im, source[i].re);
            if (i == 0) first = mag;
            if (i == N-1) last = phase;
            total += mag;
            total += j * phase;
        }
    }
    
    clock_t end = clock();

    cerr << "naive method: first = " << first << ", last = " << last
         << ", total = " << total << endl;
    cerr << "time for naive method: " << float(end - start)/divisor << endl;

    start = clock();

    first = last = total = 0;
    for (int j = 0; j < iterations; ++j) {
        v_cartesian_to_polar(mag, phase, source, N);
        for (int i = 0; i < N; ++i) {
            if (i == 0) first = mag[i];
            if (i == N-1) last = phase[i];
            total += mag[i];
            total += j * phase[i];
        }
    }

    end = clock();

    cerr << "v_cartesian_to_polar: first = " << first << ", last = " << last
         << ", total = " << total << endl;
    cerr << "time for v_cartesian_to_polar: " << float(end - start)/divisor << endl;

    deallocate(source);
    deallocate(mag);
    deallocate(phase);
    
    cerr << endl;
    return true;
}

bool
testCartToPolar()
{
    cerr << "testVectorOps: testing v_cartesian_to_polar" << endl;

    const int iterations = 50000;
    const int N = 1024;

    bq_complex_element_t *real = allocate<bq_complex_element_t>(N);
    bq_complex_element_t *imag = allocate<bq_complex_element_t>(N);
    bq_complex_element_t *mag = allocate<bq_complex_element_t>(N);
    bq_complex_element_t *phase = allocate<bq_complex_element_t>(N);

    for (int i = 0; i < N; ++i) {
        real[i] = (drand48() * 2.0) - 1.0;
        imag[i] = (drand48() * 2.0) - 1.0;
    }

    float divisor = float(CLOCKS_PER_SEC) / 1000.f;
    clock_t start = clock();

    double first, last, total = 0;
    for (int j = 0; j < iterations; ++j) {
        for (int i = 0; i < N; ++i) {
            double mag = sqrt(real[i] * real[i] + imag[i] * imag[i]);
            double phase = atan2(imag[i], real[i]);
            if (i == 0) first = mag;
            if (i == N-1) last = phase;
            total += mag;
            total += j * phase;
        }
    }
    
    clock_t end = clock();

    cerr << "naive method: first = " << first << ", last = " << last
         << ", total = " << total << endl;
    cerr << "time for naive method: " << float(end - start)/divisor << endl;

    start = clock();

    first = last = total = 0;
    for (int j = 0; j < iterations; ++j) {
        v_cartesian_to_polar(mag, phase, real, imag, N);
        for (int i = 0; i < N; ++i) {
            if (i == 0) first = mag[i];
            if (i == N-1) last = phase[i];
            total += mag[i];
            total += j * phase[i];
        }
    }

    end = clock();

    cerr << "v_cartesian_to_polar: first = " << first << ", last = " << last
         << ", total = " << total << endl;
    cerr << "time for v_cartesian_to_polar: " << float(end - start)/divisor << endl;

    deallocate(real);
    deallocate(imag);
    deallocate(mag);
    deallocate(phase);
    
    cerr << endl;
    return true;
}

int main(int, char **)
{
    cerr << endl;
    if (!testMultiply()) return 1;
    if (!testPolarToCart()) return 1;
    if (!testPolarToCartInterleaved()) return 1;
    if (!testCartToPolarComplexTypes()) return 1;
    if (!testCartToPolar()) return 1;
    return 0;
}

