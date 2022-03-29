from fedora:33

RUN dnf -y install \
    gcc-fortran \
    gcc-c++ \
    gcc \
    cmake \
    make \
    netcdf-fortran-devel \
    git \
    valgrind \
  && dnf clean all

COPY . /aerosol/

# command line argument for option build of aerosol code
ARG BUILD=ALL

# build the aerosol demo code
RUN mkdir build \
    && cd build \
    && if [ "$BUILD" = "TESTS_ONLY" ] ; then \
         echo "Building tests only" \
         && cmake -D TESTS_ONLY:BOOL=TRUE ../aerosol \
       ; else \
         echo "Building all" \
         && cmake ../aerosol \
       ; fi \
    && make