############################################################
# Dockerfile to build R image
############################################################

FROM rocker/r-base:latest

RUN R -e 'install.packages("BiocManager")' \
&& R -e 'BiocManager::install("DESeq2")'
#&& R --vanilla
#&& ls 
#&& ls /usr \
#&& ls /usr/local \
#&& ls /usr/local/bin
ENTRYPOINT R --vanilla










# Set the base image to Ubuntu
#FROM ubuntu:22.04

#RUN ln -snf /usr/share/zoneinfo/$CONTAINER_TIMEZONE /etc/localtime && echo $CONTAINER_TIMEZONE > /etc/timezone \
#&& apt-get update --fix-missing \
#&& apt-get install -y libgomp1 libcurl4-gnutls-dev libpcre2-dev g++ gfortran wget gcc mono-mcs make libbz2-dev zlib1g zlib1g-dev liblzma5 liblzma-dev libncurses5 libncurses5-dev bzip2 \
#&& cd /usr/local/ \
#&& wget https://cran.r-project.org/src/base/R-4/R-4.2.1.tar.gz \ 
#&& tar -xzvf R-4.2.1.tar.gz \
#&& rm -rf R-4.2.1.tar.gz \
#&& cd R-4.2.1 \
#&& ./configure --with-readline=no --with-x=no \
#&& make \
#&& make install \
#&& ls \
#&& cd /usr/local/bin/ \
#&& ls \
#&& cd /usr/local/ \
#&& chmod +x /usr/local/bin/R \
#&& rm -rf R-4.2.1 \
#&& apt-get remove -y wget gcc make libbz2-dev zlib1g-dev liblzma-dev libncurses-dev bzip2 \
#&& apt-get autoremove -y \
#&& apt-get clean \
#&& rm -rf /var/lib/apt/lists/*
#ENTRYPOINT ["/usr/local/bin/R"]











# Install wget and gnupg commands
#RUN  apt-get update \
 #&& apt-get install -y wget gnupg2 gnupg gnupg1 \
 #&& rm -rf /var/lib/apt/lists/* \


# Download (wget command) the public GPG key and onvert it in a format APT understands (gpg --dearmor command)
#&& wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | gpg --dearmor -o /usr/share/keyrings/r-project.gpg \

# Add R source list to sources.list.d directory = where apt looks for new sources
#&& echo "deb [signed-by=/usr/share/keyrings/r-project.gpg] https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" | tee -a /etc/apt/sources.list.d/r-project.list \


# Update the package list so APT can read the new R source
#&& apt update \

# Install R
#&& apt install -y --no-install-recommends r-base \

# Install package DESeq2
#&& R -e 'install.packages("BiocManager")' \
#&& R -e 'BiocManager::install("DESeq2")'

