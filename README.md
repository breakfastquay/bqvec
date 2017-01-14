
bqvec
=====

A small library for vector management and arithmetic in C++ using raw
C pointer arrays, designed for simple audio buffer-shuffling. Uses
vector arithmetic helpers in places, plus loops written with an eye to
compiler auto-vectorisation. Also includes aligned malloc wrappers and
a lock-free ring buffer.

This code originated as part of the Rubber Band Library written by the
same authors (see https://bitbucket.org/breakfastquay/rubberband/).
It has been pulled out into a separate library and relicensed under a
more permissive licence.

Generally expected to be vendored in to local project builds rather
than being installed as a system library.

C++ standard required: C++98 (does not use C++11 or newer features)

 * To compile on Linux: make test
 * To compile on macOS: make -f build/Makefile.osx test

[![Build Status](https://travis-ci.org/breakfastquay/bqvec.svg?branch=master)](https://travis-ci.org/breakfastquay/bqvec)

Copyright 2007-2017 Particular Programs Ltd, see COPYING for
(BSD/MIT-style) licence terms.

