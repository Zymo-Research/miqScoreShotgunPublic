FROM python:3.6.7

MAINTAINER Michael M. Weinstein, Zymo Research
LABEL version="0.0.1"

WORKDIR /

RUN apt-get update && \
    apt install -y autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libssl-dev libncurses5-dev && \
    cd tmp && \
    wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
    tar -xjvf samtools-1.9.tar.bz2 && \
    rm samtools-1.9.tar.bz2 && \
    cd samtools-1.9 && \
    ./configure && \
    make && \
    make prefix=/opt install && \
    cd /tmp && \
    rm -rf samtools-1.9

ENV PATH "$PATH:/opt/bin"


#set up BWA
RUN cd /tmp &&\
    wget https://github.com/lh3/bwa/archive/v0.7.16.tar.gz && \
    tar -xvf v0.7.16.tar.gz && \
    rm v0.7.16.tar.gz && \
    cd bwa-0.7.16 &&\
    make && \
    cd /tmp && \
    cp -r bwa-0.7.16 /opt/ && \
    rm -rf bwa-0.7.16

ENV PATH "$PATH:/opt/bwa-0.7.16"


#set up minimap2
RUN cd /tmp &&\
    wget https://github.com/lh3/minimap2/archive/v2.17.tar.gz && \
    tar -xvf v2.17.tar.gz && \
    rm v2.17.tar.gz && \
    cd minimap2-2.17 && \
    make && \
    cd /tmp && \
    cp -r minimap2-2.17 /opt/ && \
    rm -rf minimap2-2.17

ENV PATH "$PATH:/opt/minimap2-2.17"


#set up scripts
RUN cd /opt && \
    mkdir referenceBuild && \
    mkdir referenceBuild/reference && \
    mkdir miqScoreShotgun

COPY ./reference /opt/referenceBuild/reference

COPY ./requirements.txt /opt/referenceBuild

#doing expensive and unlikely to change build processes here to speed up testing builds
RUN cd /opt/referenceBuild/ && \
    pip3 install -r requirements.txt && \
    cd reference && \
    echo "Indexing standard genome" && \
    bwa index zrCommunityStandard.fa

#doing cheaper and likely to change build steps now
COPY . /opt/miqScoreShotgun

RUN cd /opt/miqScoreShotgun && \
    rm -rf reference && \
    mv /opt/referenceBuild/reference . && \
    rm -rf /opt/referenceBuild

CMD ["python3", "/opt/miqScoreShotgun/analyzeStandardReads.py"]
