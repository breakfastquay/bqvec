
# Add to VECTOR_DEFINES the relevant options for your desired
# third-party library support.
#
# Available options are
#
#  -DHAVE_IPP    Intel's Integrated Performance Primitives are available
#  -DHAVE_VDSP   Apple's Accelerate framework is available
#
# These are optional (they affect performance, not function) and you
# may define more than one of them.
# 
# Add any relevant -I flags for include paths as well.
#
# Note that you must supply the same flags when including bqvec
# headers later as you are using now when compiling the library. (You
# may find it simplest to just add the bqvec source files to your
# application's build system and not build a bqvec library at all.)

VECTOR_DEFINES	:= 

# Add any related includes and libraries here
#
THIRD_PARTY_INCLUDES	:=                   # e.g. -I/opt/intel/ipp/include
THIRD_PARTY_LIBS	:=                   # e.g. -L/opt/intel/ipp/lib/intel64_lin -Wl,-Bstatic -lipps -lippvm -lippcore -Wl,-Bdynamic


# Add to ALLOCATOR_DEFINES options relating to aligned malloc.
#
# Available options are
#
#  -DHAVE_POSIX_MEMALIGN     The posix_memalign call is available in sys/mman.h
#  -DLACK_POSIX_MEMALIGN     The posix_memalign call is not available
#
#  -DMALLOC_IS_ALIGNED       The malloc call already returns aligned memory
#  -DMALLOC_IS_NOT_ALIGNED   The malloc call does not return aligned memory
#
#  -DUSE_OWN_ALIGNED_MALLOC  No aligned malloc is available, roll your own
#
#  -DLACK_BAD_ALLOC          The C++ library lacks the std::bad_alloc exception
#
# Here "aligned" is assumed to mean "aligned enough for whatever
# vector stuff the space will be used for" which most likely means
# 16-byte alignment.
#
# The default is to use _aligned_malloc when building with Visual C++,
# system malloc when building on OS/X, and posix_memalign otherwise.
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

test:	$(LIBRARY) timings test-vectorops test-vectorops-complex
	./test-vectorops && ./test-vectorops-complex

$(LIBRARY):	$(OBJECTS)
	$(AR) rc $@ $^

timings: $(TIMINGS_OBJECTS) $(LIBRARY)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(THIRD_PARTY_LIBS)

test-vectorops:	test/TestVectorOps.o
	$(CXX) $(CXXFLAGS) -o $@ $^ -lboost_unit_test_framework $(THIRD_PARTY_LIBS)

test-vectorops-complex:	test/TestVectorOpsComplex.o
	$(CXX) $(CXXFLAGS) -o $@ $^ -lboost_unit_test_framework $(THIRD_PARTY_LIBS)

clean:		
	rm -f $(OBJECTS) $(TEST_OBJECTS) $(TIMINGS_OBJECTS)

distclean:	clean
	rm -f $(LIBRARY) test-vectorops test-vectorops-complex

depend:
	makedepend -Y -fMakefile $(SOURCES) $(TIMINGS_SOURCES) $(TEST_SOURCES) $(HEADERS)


# DO NOT DELETE

test/Timings.o: bqvec/VectorOpsComplex.h bqvec/VectorOps.h bqvec/Restrict.h
test/Timings.o: bqvec/ComplexTypes.h
test/Timings.o: bqvec/VectorOpsComplex.h bqvec/VectorOps.h bqvec/Restrict.h
test/Timings.o: bqvec/ComplexTypes.h
test/TestVectorOps.o: bqvec/VectorOps.h bqvec/Restrict.h
bqvec/RingBuffer.o: bqvec/Barrier.h bqvec/Allocators.h bqvec/Restrict.h
bqvec/RingBuffer.o: bqvec/VectorOps.h
bqvec/VectorOpsComplex.o: bqvec/VectorOps.h bqvec/Restrict.h
bqvec/VectorOpsComplex.o: bqvec/ComplexTypes.h
bqvec/VectorOps.o: bqvec/Restrict.h
