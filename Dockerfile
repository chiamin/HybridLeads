FROM gcc:12.2
LABEL maintainer="ChiaMin chiaminchung@gmail.com, TaoLin tanlin2013@gmail.com"

ARG WORKDIR=/home
ARG PKGDIR=/root
WORKDIR $WORKDIR
COPY . $WORKDIR

# Set up compiling flags for ITensor, and the shared lib path
ENV CCCOM="g++ -m64 -std=c++17 -fconcepts -fPIC" \
    PLATFORM="lapack" \
    BLAS_LAPACK_LIBFLAGS="-lpthread -L/usr/lib -lblas -llapack" \
    LD_LIBRARY_PATH="/usr/local/lib"

# Install cmake, lapack, blas
RUN apt update && \
    apt-get install -y --no-install-recommends \
    cmake \
    # gdb \
    liblapack-dev \
    liblapacke-dev \
    libopenblas-dev

# Copy external dependencies from git submodules into $PKGDIR
COPY ext $PKGDIR

# Install ITensor
RUN cd $PKGDIR/itensor && \
    cp options.mk.sample options.mk && \
    make -e

# Install armadillo
RUN cd $PKGDIR/armadillo && \
    cmake . && \
    make install

# Install glog
RUN cd $PKGDIR/glog && \
    cmake -B build -G "Unix Makefiles" && \
    cmake --build build --target install

# Install Catch2 framework for unit test
RUN cd $PKGDIR/catch2 && \
    cmake -Bbuild -H. -DBUILD_TESTING=OFF && \
    cmake --build build --target install

# Install trompeloeil
RUN cd $PKGDIR/trompeloeil && \
    cmake -B build -G "Unix Makefiles" && \
    cmake --build build --target install

RUN apt-get -y clean && \
    rm -rf /var/lib/apt/lists/*

ENTRYPOINT /bin/bash
