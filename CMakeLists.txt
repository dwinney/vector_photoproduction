cmake_minimum_required (VERSION 2.9)

project (vector-photoproduction)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

include_directories("include")
include_directories("src")
include_directories("executables")

add_executable(test
executables/test.cpp
src/polarization_vector.cpp
src/two_body_state.cpp
src/dirac_spinor.cpp
src/pomeron_exchange.cpp
)
