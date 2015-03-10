# Copyright (c) 2014, Florent Hedin, Markus Meuwly, and the University of Basel
# All rights reserved.
# The 3-clause BSD license is applied to this software.
# see LICENSE.txt

cmake_minimum_required (VERSION 2.6)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

set(TGT "mc_LJ" CACHE STRING "The executable name.")

project(${TGT} C)

# list of include folders with .h files
include_directories(
include
dSFMT
# if later -DLUA_PLUGINS is used for enabling external lua scripts for energy please add here path where to find lua.h
#/usr/local/include/luajit-2.0
#/home/hedin/bin/luajit_last/include/luajit-2.0
)

# list all source files
set(
SRCS
src/ener.c
src/io.c
src/logger.c
src/main.c
src/MCclassic.c
src/MCspav.c
src/memory.c
src/minim.c
src/parsing.c
src/plugins_lua.c
src/rand.c
src/tools.c
dSFMT/dSFMT.c
)

# gnu compiler : should perform well on all cpus
set(CMAKE_C_FLAGS_RELEASE "-std=c99 -Wall -Wextra -msse2 -fno-strict-aliasing -O3")
set(CMAKE_C_FLAGS_DEBUG   "-std=c99 -Wall -Wextra -msse2 -fno-strict-aliasing -O0 -g")

# gnu compiler : optimise for current machine
#set(CMAKE_C_FLAGS_RELEASE "-std=c99 -Wall -Wextra -msse2 -fno-strict-aliasing -O3 -march=native")
#set(CMAKE_C_FLAGS_DEBUG   "-std=c99 -Wall -Wextra -msse2 -fno-strict-aliasing -O0 -g")

# when using intel compiler : optimise for modern intel cpus
#set(CMAKE_C_FLAGS_RELEASE "-std=c99 -Wall -Wextra -msse2 -ax,CORE-AVX2,CORE-AVX-I,AVX,SSE4.2,SSE4.1,SSSE3,SSE3 -fno-strict-aliasing -O2")
#set(CMAKE_C_FLAGS_DEBUG   "-std=c99 -Wall -Wextra -msse2 -ax,CORE-AVX2,CORE-AVX-I,AVX,SSE4.2,SSE4.1,SSSE3,SSE3 -fno-strict-aliasing -O0 -g")

# when using intel compiler : optimise for current machine
#set(CMAKE_C_FLAGS_DEBUG   "-std=c99 -Wall -Wextra -fno-strict-aliasing -O0 -g")
#set(CMAKE_C_FLAGS_RELEASE "-std=c99 -Wall -Wextra -fno-strict-aliasing -O3 -xHost")

# never remove -DHAVE_SSE2 -DDSFMT_MEXP=19937 as they are necessary for the dSFMT random numbers generator
# uncomment line with -DLUA_PLUGINS if you want to enable use of external lua scripts for energy evaluation
#add_definitions(-DHAVE_SSE2 -DDSFMT_MEXP=19937 -DLUA_PLUGINS)
add_definitions(-DHAVE_SSE2 -DDSFMT_MEXP=19937)

# if -DLUA_PLUGINS was added we need the path where to find liblua*.so of libluajit*.so
#link_directories(/home/hedin/bin/luajit_last/lib/)
#link_directories(/usr/local/include/luajit-2.0/libluajit-5.1.so)

add_executable(${TGT} ${SRCS})
target_link_libraries(${TGT} m rt dl)

#target_link_libraries(${TGT} m rt dl luajit-5.1)

