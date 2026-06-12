# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file LICENSE.rst or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION ${CMAKE_VERSION}) # this file comes with cmake

# If CMAKE_DISABLE_SOURCE_CHANGES is set to true and the source directory is an
# existing directory in our source tree, calling file(MAKE_DIRECTORY) on it
# would cause a fatal error, even though it would be a no-op.
if(NOT EXISTS "/home/hsmorin/Documents/RCC/cpp/main/build/_deps/treehh-src")
  file(MAKE_DIRECTORY "/home/hsmorin/Documents/RCC/cpp/main/build/_deps/treehh-src")
endif()
file(MAKE_DIRECTORY
  "/home/hsmorin/Documents/RCC/cpp/main/build/_deps/treehh-build"
  "/home/hsmorin/Documents/RCC/cpp/main/build/_deps/treehh-subbuild/treehh-populate-prefix"
  "/home/hsmorin/Documents/RCC/cpp/main/build/_deps/treehh-subbuild/treehh-populate-prefix/tmp"
  "/home/hsmorin/Documents/RCC/cpp/main/build/_deps/treehh-subbuild/treehh-populate-prefix/src/treehh-populate-stamp"
  "/home/hsmorin/Documents/RCC/cpp/main/build/_deps/treehh-subbuild/treehh-populate-prefix/src"
  "/home/hsmorin/Documents/RCC/cpp/main/build/_deps/treehh-subbuild/treehh-populate-prefix/src/treehh-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/home/hsmorin/Documents/RCC/cpp/main/build/_deps/treehh-subbuild/treehh-populate-prefix/src/treehh-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/home/hsmorin/Documents/RCC/cpp/main/build/_deps/treehh-subbuild/treehh-populate-prefix/src/treehh-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
