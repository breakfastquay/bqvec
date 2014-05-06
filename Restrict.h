/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */
/* Copyright Chris Cannam - All Rights Reserved */

#ifndef BQ_RESTRICT_H
#define BQ_RESTRICT_H

#ifdef __MSVC__
#define R__ __restrict
#endif

#ifdef __GNUC__
#define R__ __restrict__
#endif

#ifndef R__
#define R__
#endif

#endif
