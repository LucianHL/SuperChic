# SuperChic MC event generator for central exclusive production

Version 4.22

See manual in `doc/` folder for further instructions!


The latest version of Superchic can be compiled with `CMake` build system.
Requirements: 
 - CMake > 3.16
 - Fortran compiller: GNU, Intel, flang, NVFortran were tested
 - APFEL
 - LHAPDF
 - Internet connection if the installation of PDFs was requested
 - C++ compiler, HepMC3 and Pythia8 for the tests

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


The extra flags might be:
- generic CMake flags, e.g. `-DCMAKE_INSTALL_PREFIX=/my/home/dir`, `-DCMAKE_Fortran_COMPILER=ifort`, etc.
- flags pointing to the dependencies, `-DLHAPDF_DIR=/where/the/LHAPDF/is`, `-DAPFEL_DIR=/where/the/APFEL/is`
- flags that regulate the compilation. In the current version there are two flags: `-DSUPERCHIC_ENABLE_TESTS=ON/OFF` (Enables building of tests) and
`-DSUPERCHIC_DOWNLOAD_PDFS=ON/OFF`(Enables downloading of PDFs for tests). 

## Running SuperChic

SuperChic needs a set of data files to run.
In the source and build tree those are located in `share/SuperChic` directory.
After the installation in the runtime the location of the files is `${CMAKE_INSTALL_DATADIR}/SuperChic`. Typically that is 
`/usr/share/SuperChic` or `/usr/local/share/SuperChic` directory.
The location of those files in runtime can be overwitten with an environment variable `SUPERCHIC_DATA_PATH`.
 
