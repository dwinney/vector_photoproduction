cmake_minimum_required (VERSION 2.9)

project (jpacPhoto)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_FLAGS "-fPIC -O3") 
set(CMAKE_BUILD_TYPE "Release")

# Make sure gcc version is atleast 6!
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 7.0)
    message(FATAL_ERROR "GCC version must be atleast 7.0!")
  endif()
endif()

# INSTALLATION SETTINGS
set( LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib )
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin )

# Make sure we can find ROOT
execute_process(COMMAND root-config --prefix 
                COMMAND tr -d \\n
                RESULT_VARIABLE RC_RESULT 
                OUTPUT_VARIABLE ROOT_PREFIX)

list(APPEND CMAKE_MODULE_PATH "${ROOT_PREFIX}/share/root/cmake")

## Connect ROOT
find_package(ROOT REQUIRED MathMore)
if (ROOT_FOUND)
    message("-- ROOT found!")
    include_directories(${ROOT_INCLUDE_DIRS})
    link_directories(${ROOT_LIBRARY_DIRS})
else()
    message(SEND_ERROR "-- ROOT not found!")
endif()

# BUILD LIBRARY FROM LOCAL FiLES
include_directories("include")
include_directories("src")
file(GLOB_RECURSE INC "include/*.hpp")
file(GLOB_RECURSE SRC "src/*.cpp")
add_library( jpacPhoto SHARED ${INC} ${SRC} )

# Find the jpacStyle library
find_library(JSTYLELIB NAMES jpacStyle libjpacStyle 
                       HINTS "$ENV{JPACSTYLE}/lib")
if (JSTYLELIB)
    message("-- jpacStyle found! (${JSTYLELIB})")
    include_directories("$ENV{JPACSTYLE}/include")
    link_libraries(${JSTYLELIB})
else()
    message("-- jpacStyle not found! Executables will not be installed.")
endif()

install(TARGETS jpacPhoto
  LIBRARY DESTINATION "${LIBRARY_OUTPUT_DIRECTORY}" )

# if Style is found
# complie all the executables in the bin folder
if (JSTYLELIB)
    include_directories("executables")
    file(GLOB_RECURSE EXE_FILES "executables/*.cpp")
    foreach( exefile ${EXE_FILES} )
        get_filename_component( exename ${exefile} NAME_WE)
        add_executable( ${exename} ${exefile} )
        target_link_libraries( ${exename} jpacPhoto)
        target_link_libraries( ${exename} ${ROOT_LIBRARIES})
        target_link_libraries( ${exename} ${JSTYLELIB})
    endforeach( exefile ${EXE_FILES} )
endif()