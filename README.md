# HybridLeads

Installation
------------

1. Container

    * Pull the image:
    ```
    docker pull ghcr.io/chiamin/hybridleads:main
    ```
    * Build the image by yourself:
    ```
    docker build --no-cache --force-rm -t hybridleads .
    ```
    * Run the container:
    ```
    docker run --rm -it -v $(pwd)/tests:/home/tests hybridleads
    ```
    * Package version: 
      * `gcc 12.2.0`
      * `c++17`
      * `ITensor v3`

2. Build from scratch
    
    There're few prior requirements,
    * lapack, blas
    * [ITensor](https://itensor.org/)
    * [itensor.utility](https://github.com/chiamin/itensor.utility)
