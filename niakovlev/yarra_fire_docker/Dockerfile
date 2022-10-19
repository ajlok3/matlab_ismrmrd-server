FROM niakovlev/firematlab
USER root
RUN apt update
RUN apt install -y libfftw3-dev liblapacke-dev libpng-dev libopenblas-dev

COPY bart-0.6.00 /bart-0.6.00
COPY libgpuNUFFT_ATM_f.so /usr/lib
COPY libcufft.so.9.1 /usr/lib/x86_64-linux-gnu
RUN cp /bart-0.6.00/bart /
USER 1001

ENTRYPOINT ["/usr/bin/mlrtapp/fire_matlab_ismrmrd_server"]

