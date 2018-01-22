cmake_minimum_required(VERSION 2.8)
project(CHPSO CXX)

set (CMAKE_MODULE_PATH "../cmake;${CMAKE_MODULE_PATH}")

set(ARCH "native" CACHE STRING "architecture to use for -march flag to compiler")

if(${OLD_DWARF_FORMAT})
    set(DEBUG_FLAGS "-gdwarf-3")
else()
    set(DEBUG_FLAGS "-g")
endif()

set (OMP_FLAGS "${OpenMP_CXX_FLAGS}") #"-fsanitize=address")
set (DEFINES "-DR123_NO_SINCOS -Drestrict=__restrict__") # -DNDEBUG")
set (OPT     " -O3 -ffast-math -fno-finite-math-only -fno-omit-frame-pointer -march=${ARCH}")
set (WARN    "-pedantic -Wall -Wno-strict-overflow -Wno-unused-function -Wno-unknown-pragmas")

set (CMAKE_CXX_FLAGS "${OMP_FLAGS} ${CMAKE_CXX_FLAGS} ${DEBUG_FLAGS} -std=c++11 ${DEFINES} ${OPT} ${WARN}" )
set (CMAKE_LD_FLAGS  "${OMP_FLAGS} ${CMAKE_LD_FLAGS}" )
include_directories(SYSTEM "include")

add_executable (CHPSO
    chpso.cpp
    )
