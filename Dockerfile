FROM fedora:38
ENV LC_ALL=C
ENV PYTHON=/usr/bin/python3
RUN  yum -y  install  dnf-plugins-core \
                      bc make binutils git wget cmake  diffutils file sed gawk grep which autoconf automake libtool rpm-build python3 python3-devel python3-lhapdf lhapdf lhapdf-devel \
                       texlive texlive-bibtex  'tex(multirow.sty)'   'tex(longtable.sty)' 'tex(helvet.sty)' 'tex(epsfig.sty)'  'tex(subfigure.sty)'   'tex(amsmath.sty)' 'tex(setspace.sty)' \
                      gcc-gfortran gcc-c++ bzip2 flang pythia8 pythia8-devel pythia8-lhapdf pythia8-data pythia8-doc HepMC3 HepMC3-devel openssl-devel openssl && yum -y clean all &&  lhapdf update &&  lhapdf install MSHT20qed_nnlo cteq6l1 && \
git clone https://github.com/scarrazza/apfel.git  && cd  apfel && ./configure --disable-pywrap --prefix=/usr && make && make install && cd  ../ && rm -rf apfel && \
wget https://superchic.hepforge.org/SF_MSHT20qed_nnlo.tar.gz && tar zxvf SF_MSHT20qed_nnlo.tar.gz && mv SF_MSHT20qed_nnlo /usr/share/LHAPDF/
