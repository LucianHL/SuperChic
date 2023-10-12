FROM fedora:38
ENV LC_ALL=C
RUN  yum -y  install  dnf-plugins-core \
                      bc make binutils git wget diffutils file sed gawk grep which autoconf automake libtool python3-devel python3-lhapdf lhapdf lhapdf-devel\
                      gcc-gfortran gcc-c++ bzip2   openssl-devel openssl && yum -y clean all &&  lhapdf update &&  lhapdf install MSHT20qed_nnlo && \
git clone https://github.com/scarrazza/apfel.git  && cd  apfel && ./configure --prefix=/usr && make && make install && cd  ../ && rm -rf apfel 
