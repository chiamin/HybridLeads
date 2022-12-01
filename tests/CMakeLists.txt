cmake_minimum_required(VERSION 3.10)

project(Unittest LANGUAGES CXX)

set(ITENSOR_DIR "/root/itensor")

# Include 3rd party ITensor utilities
include_directories(${ITENSOR_DIR}/ "/root/itensor.utility/")

# Include Itensor library
add_library(itensor STATIC IMPORTED)
set_property(TARGET
    itensor PROPERTY IMPORTED_LOCATION ${ITENSOR_DIR}/lib/libitensor.a
)

# Include armadillo
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

# Define sources and output executable
add_executable(test.exe
    test_one_particle_basis.cpp
)

# Compiling flags
target_compile_features(test.exe PRIVATE cxx_std_17)
target_compile_options(test.exe PRIVATE -m64 -fconcepts -fPIC)
target_link_libraries(test.exe PRIVATE pthread blas lapack)

# Include Catch2 framework for unit test
find_package(Catch2 3 REQUIRED)

# Link executable with libraries
set(ITENSOR_LIBFLAGS "-lpthread -L/usr/lib -lblas -llapack -fopenmp")
target_link_libraries(test.exe
    PRIVATE itensor "${ITENSOR_LIBFLAGS}"
    ${ARMADILLO_LIBRARIES}
    Catch2::Catch2WithMain
)
