cmake_minimum_required( VERSION 3.10 )

project( Unittest LANGUAGES CXX )

set( ITENSOR_DIR "/root/itensor" )

# Include 3rd party ITensor utilities
include_directories( ${ITENSOR_DIR}/ "/root/itensor.utility/" )

# Include Itensor library
add_library( itensor STATIC IMPORTED )
set_property( TARGET
    itensor PROPERTY IMPORTED_LOCATION ${ITENSOR_DIR}/lib/libitensor.a
)

# Define sources and output executable
add_executable( test.exe
    test_one_particle_basis.cpp
)

# Include Catch2 framework for unit test
find_package(Catch2 3 REQUIRED)

# Link executable with libraries
set( ITENSOR_LIBFLAGS "-lpthread -L/usr/lib -lblas -llapack -fopenmp" )
target_link_libraries( test.exe
    PRIVATE itensor "${ITENSOR_LIBFLAGS}" Catch2::Catch2WithMain
)
