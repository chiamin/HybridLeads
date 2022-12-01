# HybridLeads
[![docker-build](https://github.com/tanlin2013/HybridLeads/actions/workflows/build.yml/badge.svg)](https://github.com/tanlin2013/HybridLeads/actions/workflows/build.yml)
[![tests](https://github.com/tanlin2013/HybridLeads/actions/workflows/tests.yml/badge.svg)](https://github.com/tanlin2013/HybridLeads/actions/workflows/tests.yml)

Installation
------------

1. Container

    * Pull the image:
    ```
    docker pull ghcr.io/chiamin/hybridleads:latest
    ```
    * Build the image by yourself:
    ```
    docker build --no-cache --force-rm -t hybridleads .
    ```
    * Run the container:
    ```
    docker run --rm -it -v $(pwd)/tests:/home/tests ghcr.io/chiamin/hybridleads:latest
    ```
    **Note**: replace the image name `ghcr.io/chiamin/hybridleads:latest` by `hybridleads` if you're building the image by yourself.
    * Package version:
      * `gcc 12.2.0`
      * `c++17`
      * `cmake 3.10^`
      * `ITensor v3.1.11`
      * `armadillo 11.4.x`
      * `Catch2 v3.2.0`
    * Environment variables:
      The compiling flags for ITensor,
      * `CCCOM="g++ -m64 -std=c++17 -fconcepts -fPIC"`
      * `PLATFORM="lapack"`
      * `BLAS_LAPACK_LIBFLAGS="-lpthread -L/usr/lib -lblas -llapack"`

2. Build from scratch

    There're few prior requirements,
    * lapack, blas
    * [ITensor](https://itensor.org/)
    * [itensor.utility](https://github.com/chiamin/itensor.utility)
    * [armadillo](https://arma.sourceforge.net/)


Run the tests
-------------
[Catch2](https://github.com/catchorg/Catch2) is also the framework that *ITensor* adopted for unit test. It requires to compile the sources with cmake.

To compile the tests, one can do (within ```tests/``` folder)

```
cmake -B build
make -C build
```

The resulting executable is ```test.exe``` (in ```tests/build/```).


Run other main executables
--------------------------

Compiling flags are pre-set as enviroment variables in the Docker container, one can do
```
make -e
```
to use those flags.

Then, for instance, one can run the executable ```itdvp/itdvp.exe``` by
```
./itdvp.exe input
```
with the parameters been assigned in ```itdvp/input```.
