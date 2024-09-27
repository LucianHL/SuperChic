# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.20.0)

# We name the project and the target for the ExternalProject_Add() call
# to something that will highlight to the user what we are working on if
# something goes wrong and an error message is produced.

project(msht20qed_nnlo-populate NONE)



include(ExternalProject)
ExternalProject_Add(msht20qed_nnlo-populate
                     "UPDATE_DISCONNECTED" "False" "URL" "https://lhapdfsets.web.cern.ch/current/MSHT20qed_nnlo.tar.gz" "URL_HASH" "MD5=cce676b0bc21bb985c40fceac67e97ce"
                    SOURCE_DIR          "/unix/theory2/lang/git_sc/tautau_final/BUILD/_deps/msht20qed_nnlo-src"
                    BINARY_DIR          "/unix/theory2/lang/git_sc/tautau_final/BUILD/_deps/msht20qed_nnlo-build"
                    CONFIGURE_COMMAND   ""
                    BUILD_COMMAND       ""
                    INSTALL_COMMAND     ""
                    TEST_COMMAND        ""
                    USES_TERMINAL_DOWNLOAD  YES
                    USES_TERMINAL_UPDATE    YES
)


