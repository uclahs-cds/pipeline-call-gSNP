FROM perl:5.26.2 as builder

MAINTAINER Alfonso Lam <alfonsolam@mednet.ucla.edu>
ENV TERM=vt100

# Updating and Installing Linux Packages needed for this tool:
RUN sed -i 's/# \(.*multiverse$\)/\1/g' /etc/apt/sources.list && \
  apt-get update && apt-get -y upgrade && apt-get install -y \
  ant autoconf  automake  autotools-dev  build-essential  byobu  curl  gcc  git \
  htop  less  libbz2-dev  libcurl4-gnutls-dev  libgd-dev liblzma-dev  libncurses5-dev \
  libpcre++-dev  libssl-dev  libtool  llvm-dev  make  man  openjdk-8-jdk  pkg-config  \
  r-base software-properties-common unzip  uuid-dev  vim  wget  zlib1g-dev && \
  apt-get clean && \
  rm -rf /var/lib/apt/lists/* && rm -rf /var/cache/oracle-jdk8-installer

# Set PERL enviroment to make cpan to always answer yes when dependencies need to be installed:
ENV PERL_MM_USE_DEFAULT=1

# Installing Third Parties Perl Modules
# =====================================
RUN cpanm --notest Class::Load
RUN cpanm --notest Class::Load::XS
RUN cpanm --notest Moose
RUN cpanm --notest List::MoreUtils
RUN cpanm --notest YAML
RUN cpanm --notest Digest::MD5::File
RUN cpanm --notest DateTime
RUN cpanm --notest HPCI
RUN cpanm --notest IO::String
RUN cpanm --notest Carp::Always
RUN cpanm --notest MooseX::Types::Path::Class
RUN cpanm --notest namespace::autoclean
RUN cpanm --notest Path::Class

# Installing some Perl Modules that require some changes since they are too old:
# ==============================================================================

WORKDIR /src/
RUN wget https://sourceforge.net/projects/samtools/files/samtools/0.1.17/samtools-0.1.17.tar.bz2 && \
   tar -xjf samtools-0.1.17.tar.bz2 && rm -rf samtools-0.1.17.tar.bz2 && cd samtools-0.1.17 && \
   sed -i 's/-D_CURSES_LIB=1/-D_CURSES_LIB=1 -fPIC/g' Makefile &&  make

WORKDIR /src/
RUN  wget https://cpan.metacpan.org/authors/id/L/LD/LDS/Bio-SamTools-1.43.tar.gz && \
     tar -zxvf Bio-SamTools-1.43.tar.gz  && rm -rf Bio-SamTools-1.43.tar.gz  &&  cd Bio-SamTools-1.43/ && \
     sed -i 's/Please enter the location of the bam.h and compiled libbam.a files: /\/src\/samtools-0.1.17\//g'  Build.PL && \
     sed -i 's/prompt($prompt)/$prompt/g' Build.PL && \
     perl Build.PL && ./Build installdeps && ./Build manifest && ./Build install

RUN wget https://cpan.metacpan.org/authors/id/M/MI/MIROD/XML-DOM-XPath-0.14.tar.gz && \
    tar -zxvf XML-DOM-XPath-0.14.tar.gz && cd XML-DOM-XPath-0.14/ && \
    perl Makefile.PL &&  make  && make install

RUN wget https://cpan.metacpan.org/authors/id/M/MA/MANWAR/Filesys-DiskUsage-0.13.tar.gz && \
    tar -zxvf Filesys-DiskUsage-0.13.tar.gz && rm -rf Filesys-DiskUsage-0.13.tar.gz && cd Filesys-DiskUsage-0.13/ && \
    perl Makefile.PL &&  make  && make install

RUN wget https://cpan.metacpan.org/authors/id/C/CD/CDRAUG/BioPerl-1.7.6.tar.gz && \
    tar -zxvf BioPerl-1.7.6.tar.gz && rm -rf BioPerl-1.7.6.tar.gz && cd BioPerl-1.7.6/ && \
    perl Makefile.PL &&  make  && make install

ARG SGE_ROOT=/sched/sge/sge-2011.11
ADD SGE_tools /SGE_tools
WORKDIR /SGE_tools
RUN perl Build.PL && ./Build installdeps && ./Build install

ADD HPCD-SGE-0.57 /HPCD-SGE-0.57
WORKDIR /HPCD-SGE-0.57
RUN perl Build.PL && ./Build installdeps && ./Build install

WORKDIR /src/
# Set paths for GATK and Picard
# =============================
ENV GATKROOT=/src
ENV PICARDROOT=/src

# Installing GATK 3.7
# ===================
#RUN wget "https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.7-0-gcfedb67"
#RUN tar -jxf "auth?package=GATK-archive&version=3.7-0-gcfedb67"

ADD GATK-3.7-0-gcfedb67 /src/GATK-3.7-0-gcfedb67
RUN cp /src/GATK-3.7-0-gcfedb67/GenomeAnalysisTK.jar /src/GenomeAnalysisTK.jar

# Picard 2.7.1
# ============
RUN wget https://github.com/broadinstitute/picard/releases/download/2.7.1/picard.jar

WORKDIR /src/
RUN wget https://github.com/samtools/htslib/archive/1.9.tar.gz && tar -xzvf 1.9.tar.gz && rm 1.9.tar.gz
WORKDIR /src/htslib-1.9
RUN autoheader && autoconf && ./configure && make && make install

WORKDIR /src/
RUN wget https://github.com/samtools/samtools/archive/1.9.tar.gz && tar -xzvf 1.9.tar.gz && rm 1.9.tar.gz
WORKDIR /src/samtools-1.9
RUN  autoheader && autoconf -Wno-syntax && ./configure && make && make install

# Install some required BL modules
# ================================

ADD BoutrosLab-Utilities-General /src/BoutrosLab-Utilities-General
WORKDIR /src/BoutrosLab-Utilities-General
RUN perl Build.PL && ./Build && ./Build installdeps && ./Build install

ADD BoutrosLab-Base /src/BoutrosLab-Base
WORKDIR /src/BoutrosLab-Base
RUN perl Build.PL && ./Build && ./Build installdeps && ./Build install

ADD NGS-Tools-Utils /src/NGS-Tools-Utils
WORKDIR /src/NGS-Tools-Utils
RUN perl Build.PL &&  ./Build &&  ./Build installdeps &&  ./Build install

ADD BLSupport-Roles /src/BLSupport-Roles
WORKDIR /src/BLSupport-Roles
RUN perl Build.PL && ./Build && ./Build installdeps && ./Build install

# Install the NGS-Tools Modules:
# ==============================

ADD NGS-Tools-Picard /src/NGS-Tools-Picard
WORKDIR /src/NGS-Tools-Picard
RUN perl Build.PL && ./Build && ./Build installdeps && ./Build install

ADD NGS-Tools-ContEst /src/NGS-Tools-ContEst
WORKDIR /src/NGS-Tools-ContEst
RUN perl Build.PL && ./Build && ./Build installdeps && ./Build install

ADD NGS-Tools-GATK /src/NGS-Tools-GATK
RUN sed -i 's/\x27-an MQRankSum\x27/\#\x27-an MQRankSum\x27/g' /src/NGS-Tools-GATK/lib/NGS/Tools/GATK.pm
RUN sed -i 's/\x27-an ReadPosRankSum\x27/\#\x27-an ReadPosRankSum\x27/g' /src/NGS-Tools-GATK/lib/NGS/Tools/GATK.pm
WORKDIR /src/NGS-Tools-GATK
RUN perl Build.PL && ./Build && ./Build installdeps && ./Build install

RUN cp /src/samtools-1.9/samtools  /usr/bin/samtools 

WORKDIR /src/

# Defining the Command Prompt:

CMD ["/bin/bash -c"]
