FROM ubuntu:18.04

MAINTAINER Bernat del Olmo bdelolmo@gencardio.com

ENV DEBIAN_FRONTEND noninteractive

# install required packages
RUN apt-get update && apt-get install -y \
    build-essential \
    g++ \
    gfortran \
    git \
    make \
    make \
    zlib1g \
    libncurses5-dev \
    libncursesw5-dev \
    libboost-program-options-dev \
    libboost-system-dev \
    libboost-filesystem-dev \
    libboost-all-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    python3-pip \
    python-pip \
    cpanminus \
    samtools \
    tabix \
    bcftools \
    gnupg2 \
    gnupg1

# Utilities
RUN apt-get update && apt-get install -y gnupg ca-certificates pandoc

RUN apt-get install add-apt-key

# Add R 3.5 repository
RUN echo 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/' >> /etc/apt/sources.list.d/cran35.sources.list
RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9

RUN apt-get update --allow-unauthenticated
RUN apt-get -y install r-base r-base-core r-recommended

RUN apt-get -y install wget

RUN apt-get -y install python3.5

RUN apt-get -y install python-numpy

RUN cpanm Parallel::ForkManager \
    Statistics::Descriptive \
    Sort::Key::Natural \
    Chart::Gnuplot \
    JSON::MaybeXS \
    DBD::SQLite \
    Excel::Writer::XLSX

#RUN pip install --upgrade pip
#RUN pip install numpy
RUN apt-get -y install macs
#RUN pip install MACS2
RUN R -e "install.packages('BiocManager',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN Rscript -e 'BiocManager::install("DNAcopy")'
RUN Rscript -e 'BiocManager::install("aroma.light")'
RUN R -e "install.packages('argparse',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('gridExtra',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('ggplot2',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('RColorBrewer',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('reshape2',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('egg',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('PSCBS',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('tidyr',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('data.table',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('dplyr',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('GenomicRanges',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN Rscript -e 'BiocManager::install("GenomicRanges")'

RUN cd /opt \
    && git clone https://github.com/arq5x/bedtools2.git \
    && cd bedtools2/ \
    && make \
    && make install


#install GRAPES latest version
RUN cd /opt \
    && git clone https://github.com/bdolmo/GRAPES.git \
    && cd /opt/GRAPES/ \
    && perl INSTALL.PL \
    && echo "finished" && echo "done" && echo "redone333"

# Workdir
WORKDIR /root/

# Add GRAPES to PATH
ENV PATH="/opt/GRAPES:${PATH}"

# by default GRAPES is executed
CMD ["GRAPES"]
