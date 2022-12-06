cmake_minimum_required(VERSION 3.10)

project(Unittest LANGUAGES CXX)

set(ITENSOR_DIR "/root/itensor")

# Include 3rd party headers from itensor.utility
include_directories("/root/itensor.utility/")

# Include Itensor library
add_library(itensor STATIC IMPORTED)
set_property(TARGET
    itensor PROPERTY IMPORTED_LOCATION ${ITENSOR_DIR}/lib/libitensor.a
)
include_directories(${ITENSOR_DIR}/)

# Include armadillo
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

# Include Catch2 framework for unit test
find_package(Catch2 3 REQUIRED)

# Name all your sources filename here without extension
set(SOURCES_FILES
    test_one_particle_basis
    test_hamiltonian
)
foreach(name ${SOURCES_FILES})
    # Define sources and output executable
    add_executable(${name}.exe
        ${name}.cpp
    )

    # Compiling flags
    target_compile_features(${name}.exe PRIVATE cxx_std_17)
    target_compile_options(${name}.exe PRIVATE -m64 -fconcepts -fPIC)
    target_link_libraries(${name}.exe PRIVATE pthread blas lapack)
    set(ITENSOR_LIBFLAGS "-lpthread -L/usr/lib -lblas -llapack -fopenmp")

    # Link executable with libraries
    target_link_libraries(${name}.exe
        PRIVATE itensor "${ITENSOR_LIBFLAGS}"
        ${ARMADILLO_LIBRARIES}
        Catch2::Catch2WithMain
    )
endforeach()
