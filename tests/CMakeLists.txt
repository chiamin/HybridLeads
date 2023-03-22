cmake_minimum_required(VERSION 3.10)

project(Unittest LANGUAGES CXX)
enable_testing()
set(ITENSOR_DIR "/root/itensor")

# Include module headers
include_directories("/home/hybridleads/")

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

# Include glog for logging
find_package (glog 0.6.0 REQUIRED)

# Include Catch2 framework for unit test
find_package(Catch2 3 REQUIRED)

# Include trompeloeil for mocking in unit test
find_package(trompeloeil REQUIRED)

# Name all your sources filename here without extension
set(SOURCES_FILES
    kbasis/test_one_particle_basis.cpp
    kbasis/test_hamiltonian.cpp
    # itdvp/test_itdvp.cpp
    hybridbasis/test_mpo_model.cpp
)
option(USE_COVERAGE "Use coverage reports" OFF)
foreach(SOURCE_PATH ${SOURCES_FILES})
    # Define sources and output executable
    get_filename_component(SOURCE ${SOURCE_PATH} NAME)
    string(REGEX REPLACE "[.]cpp$" ".exe" OUTPUT_EXE ${SOURCE})
    add_executable(${OUTPUT_EXE} ${SOURCE_PATH})

    # Compiling flags
    target_compile_features(${OUTPUT_EXE} PRIVATE cxx_std_17)
    target_compile_options(${OUTPUT_EXE} PRIVATE -m64 -fconcepts -fPIC)
    target_link_libraries(${OUTPUT_EXE} PRIVATE pthread blas lapack)
    set(ITENSOR_LIBFLAGS "-lpthread -L/usr/lib -lblas -llapack -fopenmp")
    if(USE_COVERAGE)
        SET(GCC_COVERAGE_COMPILE_FLAGS "-fprofile-arcs -ftest-coverage")
        SET(GCC_COVERAGE_LINK_FLAGS "--coverage")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")
        SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}")
    endif(USE_COVERAGE)

    # Link executable with libraries
    target_link_libraries(${OUTPUT_EXE}
        PRIVATE itensor "${ITENSOR_LIBFLAGS}"
        ${ARMADILLO_LIBRARIES}
        glog::glog
        Catch2::Catch2WithMain
        trompeloeil::trompeloeil
    )

    add_test(NAME ${OUTPUT_EXE} COMMAND ${OUTPUT_EXE})
endforeach()
unset(USE_COVERAGE CACHE)
