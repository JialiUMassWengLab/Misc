FROM cogniteev/oracle-java:java8
#FROM ubuntu:12.04

## Install dependencies
RUN apt-get update && \
    apt-get install -y \
    byacc \
    bzip2 \
    g++ \
    gcc \
    git \
    make \
    python2.7-dev \
    python-matplotlib \
    python-numpy \
    python-pip \
    python-software-properties \
    r-base \
#    samtools \
    unzip \
    wget \
    zlib1g-dev 

#install HTSeq
RUN wget https://pypi.python.org/packages/72/0f/566afae6c149762af301a19686cd5fd1876deb2b48d09546dbd5caebbb78/HTSeq-0.6.1.tar.gz#md5=b7f4f38a9f4278b9b7f948d1efbc1f05 && \
    tar -zxvf HTSeq-0.6.1.tar.gz && \
    pip install pysam && \
    cd HTSeq-0.6.1/ && \
    python setup.py build && \
    python setup.py install
    
#install trim-galore
RUN wget http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/trim_galore_v0.4.2.zip && \
    unzip trim_galore_v0.4.2.zip

#install picard
#RUN apt-get update -y && \
#    add-apt-repository ppa:webupd8team/java && \
#    apt-get update && \
#    apt-get install oracle-java8-installer -y && \
RUN git clone https://github.com/broadinstitute/picard.git && \
    cd picard/ && \
    ./gradlew shadowJar

#install PEAR
RUN wget http://sco.h-its.org/exelixis/web/software/pear/files/pear-0.9.10-bin-64.tar.gz && \
    tar -zxvf pear-0.9.10-bin-64.tar.gz

#install bioawk
RUN git clone git://github.com/lh3/bioawk.git && \
    cd bioawk && \
    make

#install STAR
RUN git clone https://github.com/alexdobin/STAR.git
#RUN git clone https://github.com/alexdobin/STAR.git && \
#    cd STAR/source/ && \
#    make

#install samtools 1.3.1
RUN wget http://downloads.sourceforge.net/project/samtools/samtools/1.3.1/samtools-1.3.1.tar.bz2 && \
    bunzip2 samtools-1.3.1.tar.bz2 && \
    tar -xvf samtools-1.3.1.tar && \
    cd samtools-1.3.1 && \
    ./configure && \
    make && \
    make install

#install RSeQC
RUN pip install RSeQC
#RUN wget https://sourceforge.net/projects/rseqc/files/latest/download && \
#    mv download RSeQC-2.6.4.tar.gz && \
#    tar -zxvf RSeQC-2.6.4.tar.gz && \
#    cd RSeQC-2.6.4 && \
#    pip install -U pip setuptools && \
#    python setup.py install
    
#install FastQC
RUN wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip && \
    unzip fastqc_v0.11.5.zip && \
    cd FastQC/ && \
    chmod 755 fastqc

#install cutadapt
RUN pip install --user --upgrade cutadapt

#install septk
RUN git clone https://github.com/lh3/seqtk.git && \
    cd seqtk/ && \
    make


##Add scripts
COPY plotXYHist.r /data
COPY Dockerfile /data


##Delete intermediate files
RUN rm -R LICENSE \
       	  *.txt \
	  *.pdf \
	  *.zip \
	  *.gz \
	  *.tar


##Move executables to /bin
RUN mkdir /rbin/ && \
    cp bioawk/bioawk /rbin && \
    cp ~/.local/bin/cutadapt /rbin && \
    cp FastQC/fastqc /rbin && \
    cp STAR/bin/Linux_x86_64_static/STAR /rbin && \
    cp trim_galore /rbin && \
    cp samtools-1.3.1/samtools /rbin && \
    cp pear-0.9.10-bin-64/pear-0.9.10-bin-64 /rbin && \
    cp picard/build/libs/picard.jar /rbin && \
    cp seqtk/seqtk /rbin


##External software version info
RUN echo "    bioawk : Aug 27,2013-commit-5e8b41d\n \
    cutadapt : version 1.11\n \
    FastQC : version 0.11.5\n \
    HTSeq : version 0.6.1\n \
    STAR : version 2.5.2b\n \
    trim_galore : version 0.4.2\n \
    samtools : version 1.3.1\n \
    pear : version 0.9.10\n \
    picard : version 2.7.1\n \
    RSeQC : version 2.6.4\n \
    seqtk : version 1.2-r95-dirty\n" \
    | cat > versionInfo
	 

CMD ["sh", "-c", "export PATH=$PATH:/data/bioawk:/data/STAR/bin/Linux_x86_64:/data:/data/pear-0.9.10-bin-64"]
