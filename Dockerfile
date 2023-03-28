FROM ubuntu:22.04

MAINTAINER Dmitrii Smirnov "Dmitrii.Smirnov@skoltech.edu"

SHELL ["/bin/sh", "-c"]


RUN apt-get update && apt-get -y --no-install-recommends install \
    build-essential \
    ssh \
    sudo \
    clang \
    cmake \
    mpich \
    libboost-all-dev \
    wget

# Install miniconda
ENV CONDA_DIR /opt/conda
RUN wget --no-check-certificate --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && chmod a+x ~/miniconda.sh && ~/miniconda.sh -b -p /opt/conda

# Put conda in path so we can use conda activate
ENV PATH=$CONDA_DIR/bin:$PATH

RUN conda -V


RUN wget --no-check-certificate --quiet http://arcuda.skoltech.ru/~d.smirnov/optimalTAD.tar.gz && tar xzf optimalTAD.tar.gz && cd optimalTAD/armatus && mpicxx -std=c++11 -w -L /usr/include/boost/lib/ -I include/ -I /usr/include/boost/include/ -O3 -o binaries/armatus src/*.cpp -lboost_iostreams -lboost_program_options -lboost_system && wget --no-check-certificate --quiet http://arcuda.skoltech.ru/~d.smirnov/chr2L_control.txt.gz && cp binaries/armatus ../optimalTAD/ && cd .. && pip install .

#RUN groupadd -r sample && useradd -r -g sample sample
#USER sample

ENV user lg

#RUN useradd -m -d /home/${user} ${user} && \
#    chown -R ${user} /home/${user} && \
#    adduser ${user} sudo && \
#    echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers

#USER ${user}

RUN cd optimalTAD && optimalTAD -v

WORKDIR optimalTAD/

#RUN optimalTAD run

RUN mpirun --allow-run-as-root -np 1 ./optimalTAD/armatus -r 20000 -i armatus/chr2L_control.txt.gz -g 1.2 -o output/test_gamma -m -s 0.2

CMD ["optimalTAD", "run"]





