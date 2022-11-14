FROM gcc:12.2
LABEL maintainer="ChiaMin chiaminchung@gmail.com, TaoLin tanlin2013@gmail.com"

ARG WORKDIR=/home
ARG PKGDIR=/root

# Set up compiling flags for ITensor
ENV CCCOM="g++ -m64 -std=c++17 -fconcepts -fPIC" \
    PLATFORM="lapack" \
    BLAS_LAPACK_LIBFLAGS="-lpthread -L/usr/lib -lblas -llapack"

# Install cmake, lapack, blas
RUN apt update && \
    apt-get install -y --no-install-recommends \
    cmake \
    liblapack-dev \
    liblapacke-dev \
    libopenblas-dev

# Install ITensor
RUN git clone --branch v3 https://github.com/ITensor/ITensor.git $PKGDIR/itensor && \
    cd $PKGDIR/itensor && \
    cp options.mk.sample options.mk && \
    make -e

# Download 3rd party ITensor utilities
RUN git clone https://github.com/chiamin/itensor.utility.git $PKGDIR/itensor.utility

# Install Catch2 framework for unit test
RUN git clone https://github.com/catchorg/Catch2.git $PKGDIR/catch2 && \
    cd $PKGDIR/catch2 && \
    cmake -Bbuild -H. -DBUILD_TESTING=OFF && \
    cmake --build build/ --target install

RUN apt-get -y clean && \
    rm -rf /var/lib/apt/lists/*

WORKDIR $WORKDIR
COPY . $WORKDIR

ENTRYPOINT /bin/bash
