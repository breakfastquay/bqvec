
# Add to VECTOR_DEFINES the relevant options for your desired
# third-party library support.
#
# Available options are
#
#  -DHAVE_IPP    Intel's Integrated Performance Primitives are available
#  -DHAVE_VDSP   Apple's Accelerate framework is available
#
# The above are optional (they affect performance, not function) and
# you may define more than one of them.
#
# The following two options trade off speed against precision for single-
# precision paths in cases where IPP and VDSP are not available:
#
#  -DUSE_POMMIER_MATHFUN    Use Julien Pommier's SSE/NEON implementation
#                of sincos in 32-bit polar-to-cartesian conversion
#  -DUSE_APPROXIMATE_ATAN2  Use a quick atan2 approximation in 32-bit
#                cartesian-to-polar conversion
#
# Add any relevant -I flags for include paths as well.
#
# Note that you must supply the same flags when including bqvec
# headers later as you are using now when compiling the library. (You
# may find it simplest to just add the bqvec source files to your
# application's build system and not build a bqvec library at all.)

VECTOR_DEFINES	:= # -DHAVE_IPP

# Add any related includes and libraries here
#
THIRD_PARTY_INCLUDES	:= # -I/opt/intel/ipp/include
THIRD_PARTY_LIBS	:= # -L/opt/intel/ipp/lib/intel64_lin -Wl,-Bstatic -lipps -lippvm -lippcore -Wl,-Bdynamic


# Add to ALLOCATOR_DEFINES options relating to aligned malloc.
# These are not usually necessary.
#
# Available options are
#
#  -DHAVE_POSIX_MEMALIGN       The posix_memalign call is available in sys/mman.h
#  -DLACK_POSIX_MEMALIGN       The posix_memalign call is not available
#
#  -DMALLOC_IS_ALIGNED         The malloc call already returns aligned memory
#  -DMALLOC_IS_NOT_ALIGNED     The malloc call does not return aligned memory
#
#  -DUSE_OWN_ALIGNED_MALLOC    If no aligned malloc is available, roll your own
#  -DAVOID_OWN_ALIGNED_MALLOC  If no aligned malloc is available, refuse to build
#
#  -DLACK_BAD_ALLOC            The C++ library lacks the std::bad_alloc exception
#
# Here "aligned" is assumed to mean "aligned enough for whatever
# vector stuff the space will be used for" which likely means at least
# 16-byte alignment.
#
# If no options are provided, we will use IPP functions if HAVE_IPP is
# defined, or else use _aligned_malloc when building with Visual C++
# on Windows, roll our own when building with some other compiler on
# Windows, use system malloc when building on OS/X, and use
# posix_memalign elsewhere.
#
# Note that you must supply the same flags when including bqvec
# headers later as you are using now when compiling the library. (You
# may find it simplest to just add the bqvec source files to your
# application's build system and not build a bqvec library at all.)

ALLOCATOR_DEFINES := 


SRC_DIR	:= src
TEST_DIR := test
HEADER_DIR := bqvec

SOURCES	:= $(wildcard $(SRC_DIR)/*.cpp)
HEADERS	:= $(wildcard $(HEADER_DIR)/*.h) $(wildcard $(SRC_DIR)/*.h)

OBJECTS	:= $(SOURCES:.cpp=.o)
OBJECTS	:= $(OBJECTS:.c=.o)

TIMINGS_SOURCES	:= $(TEST_DIR)/Timings.cpp
TIMINGS_OBJECTS	:= $(TIMINGS_SOURCES:.cpp=.o)

TEST_SOURCES	:= $(wildcard $(TEST_DIR)/*.cpp)
TEST_OBJECTS	:= $(TEST_SOURCES:.cpp=.o)

CXXFLAGS := $(VECTOR_DEFINES) $(ALLOCATOR_DEFINES) -I. $(THIRD_PARTY_INCLUDES) -I$(HEADER_DIR) -O3 -ffast-math -Wall -Werror -fpic -std=c++98

LIBRARY	:= libbqvec.a

all:	$(LIBRARY) timings

test:	$(LIBRARY) timings test-allocators test-vectorops test-vectorops-complex
	./test-allocators && ./test-vectorops && ./test-vectorops-complex

valgrind:	$(LIBRARY) timings test-allocators test-vectorops test-vectorops-complex
	valgrind ./test-allocators && valgrind ./test-vectorops && valgrind ./test-vectorops-complex

$(LIBRARY):	$(OBJECTS)
	$(AR) rc $@ $^

timings: $(TIMINGS_OBJECTS) $(LIBRARY)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(THIRD_PARTY_LIBS)

test-allocators:	test/TestAllocators.o $(LIBRARY)
	$(CXX) $(CXXFLAGS) -o $@ $^ -lboost_unit_test_framework -L. -lbqvec $(THIRD_PARTY_LIBS)

test-vectorops:	test/TestVectorOps.o $(LIBRARY)
	$(CXX) $(CXXFLAGS) -o $@ $^ -lboost_unit_test_framework -L. -lbqvec $(THIRD_PARTY_LIBS)

test-vectorops-complex:	test/TestVectorOpsComplex.o $(LIBRARY)
	$(CXX) $(CXXFLAGS) -o $@ $^ -lboost_unit_test_framework -L. -lbqvec $(THIRD_PARTY_LIBS)

clean:		
	rm -f $(OBJECTS) $(TEST_OBJECTS) $(TIMINGS_OBJECTS)

distclean:	clean
	rm -f $(LIBRARY) test-allocators test-vectorops test-vectorops-complex

depend:
	makedepend -Y -fMakefile $(SOURCES) $(TIMINGS_SOURCES) $(TEST_SOURCES) $(HEADERS)


# DO NOT DELETE

test/Timings.o: bqvec/VectorOpsComplex.h bqvec/VectorOps.h bqvec/Restrict.h
test/Timings.o: bqvec/ComplexTypes.h
test/Timings.o: bqvec/VectorOpsComplex.h bqvec/VectorOps.h bqvec/Restrict.h
test/Timings.o: bqvec/ComplexTypes.h
test/TestVectorOpsComplex.o: bqvec/VectorOpsComplex.h bqvec/VectorOps.h
test/TestVectorOpsComplex.o: bqvec/Restrict.h bqvec/ComplexTypes.h
test/TestVectorOpsComplex.o: bqvec/VectorOps.h
test/TestVectorOps.o: bqvec/VectorOps.h
test/TestAllocators.o: bqvec/VectorOps.h bqvec/Allocators.h
bqvec/RingBuffer.o: bqvec/Barrier.h bqvec/Allocators.h bqvec/Restrict.h
bqvec/RingBuffer.o: bqvec/VectorOps.h
bqvec/VectorOpsComplex.o: bqvec/VectorOps.h bqvec/Restrict.h
bqvec/VectorOpsComplex.o: bqvec/ComplexTypes.h
bqvec/VectorOps.o: bqvec/Restrict.h
