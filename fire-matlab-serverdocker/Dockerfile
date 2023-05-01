FROM matlabruntime/r2020b/release/update2/c06e401208000000

RUN apt update
RUN apt install -y libfftw3-dev liblapacke-dev libpng-dev libopenblas-dev

COPY ./applicationFilesForMATLABCompiler /usr/bin/mlrtapp
COPY ./files/bart-0.6.00/ /bart-0.6.00/
COPY ./files/libgpuNUFFT_ATM_f.so /usr/lib
COPY ./files/libcufft.so.9.1 /usr/lib/x86_64-linux-gnu

RUN cp /bart-0.6.00/bart /

ENTRYPOINT ["/usr/bin/mlrtapp/fire_matlab_ismrmrd_server"]
