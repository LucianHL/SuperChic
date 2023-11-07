# SuperChic MC event generator for central exclusive production

Version 4.22

See manual in `doc/` folder for further instructions!


The latest version of Superchic can be compiled with `CMake` build system.
Build requirements: 
 - CMake (https://cmake.org/) > 3.16
 - Fortran compiller: GNU, Intel, flang, NVFortran were tested
 - APFEL, https://apfel.hepforge.org/
 - LHAPDF, see https://lhapdf.hepforge.org/
 - Internet connection if the installation of PDFs was requested
 - C++ compiler, HepMC3 (https://ep-dep-sft.web.cern.ch/project/hepmc) and Pythia8 (https://pythia.org/) for the tests
Runtime requirements:
 - APFEL, https://apfel.hepforge.org/
 - LHAPDF, see https://lhapdf.hepforge.org/


## Obtaining the required dependencies

### CentOS8

    Most of the dependencies are present in the default repositories
    
    ```
    yum -y install epel-release
    dnf config-manager --set-enabled PowerTools
    yum -y install  gcc gcc-c++ gcc-gfortran make cmake cmake-data git
    yum -y install HepMC3 HepMC3-devel HepMC3-search lhapdf lhapdf-devel python-lhapdf pythia8 pythia8-devel
    ```
    
    But `apfel` should be installed from the sources.

    ```
    git clone https://github.com/scarrazza/apfel.git  && cd  apfel && ./configure --disable-pywrap --prefix=/where_you_want && make && make install && cd  ../ && rm -rf apfel
    ```

### Fedora

    Most of the dependencies are present in the default repositories
    ```
    yum -y install  gcc gcc-c++ gcc-gfortran make cmake cmake-data git
    yum -y install HepMC3 HepMC3-devel HepMC3-search lhapdf lhapdf-devel python-lhapdf pythia8 pythia8-devel
    ```

    But `apfel` should be installed from the sources.

    ```
    git clone https://github.com/scarrazza/apfel.git  && cd  apfel && ./configure --disable-pywrap --prefix=/where_you_want && make && make install && cd  ../ && rm -rf apfel
    ```

### MacOS X
    All the dependencies are present in the `homebrew`(https://brew.sh/) repositories
    ```
    brew tap davidchall/hep
    brew install gcc cmake hepmc3 apfel lhapdf
    ```

## Building SuperChic

```
cmake -S . -B BUILD <extra flags>
cmake --build BUILD
cmake --install BUILD
```

and optionally

```
ctest --test-dir BUILD
```
to run all the tests.


The extra flags might be:
- generic CMake flags, e.g. `-DCMAKE_INSTALL_PREFIX=/my/home/dir`, `-DCMAKE_Fortran_COMPILER=ifort`, `-DCMAKE_Fortran_FLAGS="-O2 -g"`, etc.
- flags pointing to the dependencies, `-DLHAPDF_DIR=/where/the/LHAPDF/is`, `-DAPFEL_DIR=/where/the/APFEL/is`
- flags that regulate the compilation. In the current version there are the following flags: 
 - `-DSUPERCHIC_ENABLE_TESTS=ON/OFF`     Enables building of tests.
 - `-SUPERCHIC_ENABLE_ALL_TESTS=ON/OFF`  Enables building of much more tests.
 - `-DSUPERCHIC_DOWNLOAD_PDFS=ON/OFF`    Enables downloading of PDFs for tests. Makes sense only when the testing is enabled.
 - `-DSUPERCHIC_ENABLE_FPES=ON/OFF`      Enables floating point exceptions in the code. 

## Running SuperChic

The options to run SuperChic, i.e. executables `superchic` and `init` are desribed in the manual.

### Environment variables
SuperChic needs a set of data files to run.
In the source and build tree those are located in `share/SuperChic` directory.
After the installation in the runtime the location of the files is `${CMAKE_INSTALL_DATADIR}/SuperChic`. Typically that is 
`/usr/share/SuperChic` or `/usr/local/share/SuperChic` directory.
The location of those files in runtime can be overwitten with an environment variable `SUPERCHIC_DATA_PATH`.


The location of the PDFs can be set with the environment variable `LHAPDF_DATA_PATH`, see the LHAPDF socumentation. 
This is done automatically in the st suite.


## Debugging SuperChic

For the debug purposes it is useful to compile SeuperChic with the following options (assuming GNU compilers)
`-DSUPERCHIC_ENABLE_FPES=ON -DSUPERCHIC_ENABLE_TESTS=ON -SUPERCHIC_ENABLE_ALL_TESTS=ON -DSUPERCHIC_DOWNLOAD_PDFS=ON -DCMAKE_Fortran_FLAGS="-O2 -g" `

