FROM --platform=linux/amd64 ubuntu:22.04
RUN apt-get update --fix-missing \
&& apt-get install -y wget \
&& cd /usr/local/ \
&& wget https://github.com/ShiLab-Bioinformatics/subread/releases/download/2.0.2/subread-2.0.2-Linux-x86_64.tar.gz \
&& tar -xzvf subread-2.0.2-Linux-x86_64.tar.gz \
&& cp subread-2.0.2-Linux-x86_64/bin/featureCounts /usr/local/bin/featureCounts \
&& chmod +x /usr/local/bin/featureCounts \
&& rm subread-2.0.2-Linux-x86_64.tar.gz \
&& rm -rf subread-2.0.2-Linux-x86_64
ENTRYPOINT ["/usr/local/bin/featureCounts"]
