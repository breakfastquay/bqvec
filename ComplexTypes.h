/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */
/* Copyright Chris Cannam - All Rights Reserved */

#ifndef BQ_COMPLEX_TYPES_H
#define BQ_COMPLEX_TYPES_H

namespace breakfastquay {

#ifndef NO_COMPLEX_TYPES
// Convertible with other complex types that store re+im consecutively
typedef double bq_complex_element_t;
typedef struct {
    bq_complex_element_t re;
    bq_complex_element_t im;
} bq_complex_t;
#endif

}

#endif
