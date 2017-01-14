
VECTOR_DEFINES		:= -DHAVE_IPP

THIRD_PARTY_INCLUDES	:= -I/opt/intel/ipp/include
THIRD_PARTY_LIBS	:= -L/opt/intel/ipp/lib/intel64_lin -Wl,-Bstatic -lipps -lippvm -lippcore -Wl,-Bdynamic

ALLOCATOR_DEFINES 	:= 

include build/Makefile.inc

