# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.20.0)

# We name the project and the target for the ExternalProject_Add() call
# to something that will highlight to the user what we are working on if
# something goes wrong and an error message is produced.

project(cteq6l1-populate NONE)



include(ExternalProject)
ExternalProject_Add(cteq6l1-populate
                     "UPDATE_DISCONNECTED" "False" "URL" "https://lhapdfsets.web.cern.ch/current/cteq6l1.tar.gz" "URL_HASH" "MD5=5611f1e9235151d9f67254aeb13bb65f"
                    SOURCE_DIR          "/unix/theory2/lang/git_sc/tautau_final/BUILD/_deps/cteq6l1-src"
                    BINARY_DIR          "/unix/theory2/lang/git_sc/tautau_final/BUILD/_deps/cteq6l1-build"
                    CONFIGURE_COMMAND   ""
                    BUILD_COMMAND       ""
                    INSTALL_COMMAND     ""
                    TEST_COMMAND        ""
                    USES_TERMINAL_DOWNLOAD  YES
                    USES_TERMINAL_UPDATE    YES
)


