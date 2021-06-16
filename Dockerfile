FROM ubuntu:groovy

RUN apt-get update && apt-get install -y wget build-essential cmake git python3-dev
RUN apt-get install -y libgmp-dev zlib1g-dev nlohmann-json3-dev libssl-dev
RUN apt-get install -y --no-install-recommends libntl-dev

WORKDIR /home
RUN wget https://boostorg.jfrog.io/artifactory/main/release/1.73.0/source/boost_1_73_0.tar.bz2
RUN tar --bzip2 -xf boost_1_73_0.tar.bz2 && rm boost_1_73_0.tar.bz2
RUN cd /home/boost_1_73_0 && ./bootstrap.sh --with-python=/usr/bin/python3
RUN cd /home/boost_1_73_0 && ./b2 install

WORKDIR /home
RUN git clone https://github.com/rogersce/cnpy.git \
 && cd cnpy \
 && mkdir build \
 && cd build \
 && cmake .. \
 && make \
 && make install

WORKDIR /home
RUN wget http://hms.isi.jhu.edu/acsc/libpaillier/libpaillier-0.8.tar.gz \
 && tar -xzvf libpaillier-0.8.tar.gz \
 && rm libpaillier-0.8.tar.gz \
 && cd libpaillier-0.8 \
 && ./configure --with-gmp-include=/usr/include/x86_64-linux-gnu --with-gmp-lib=/usr/lib/x86_64-linux-gnu \
 && make install

WORKDIR /home
RUN git clone https://github.com/encryptogroup/ABY.git \
 && cd ABY \
 && mkdir build \
 && cd build \
 && cmake .. \
 && make \
 && make install \
 && sed '9d' -i /usr/local/lib/cmake/ABY/ABYConfig.cmake

WORKDIR /home
RUN wget https://github.com/relic-toolkit/relic/archive/refs/tags/relic-toolkit-0.5.0.tar.gz \
 && tar -xzvf relic-toolkit-0.5.0.tar.gz \
 && rm relic-toolkit-0.5.0.tar.gz \
 && cd relic-relic-toolkit-0.5.0 \
 && mkdir build \
 && cd build \
 && cmake .. \
 && make \
 && make install

RUN apt-get install -y python3-pip
RUN pip install numpy scikit-learn kneed

CMD ["bash"]
