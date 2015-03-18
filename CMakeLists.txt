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
)

# if Lua scripts are required for custom energy evaluation, we may need here the path where to find lua.h
if (USE_LUA)
include_directories(
    #/usr/include
    #/usr/local/include/lua53
    /usr/local/include/luajit-2.0
    #/usr/local/lua-5.3.0/include
)
endif()

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

# never remove -DHAVE_SSE2 -DDSFMT_MEXP=19937 as they are necessary for the dSFMT random numbers generator
add_definitions(-DHAVE_SSE2 -DDSFMT_MEXP=19937)
# if you want to siable dSFMT and use the standard C random numbers (NOT RECOMMENDED), comment previous line and uncomment the following :
#add_definitions(-DSTDRAND)

if (USE_LUA)
    # if user wants to enable use of external lua scripts for energy evaluation
    add_definitions(-DLUA_PLUGINS)

    # if -DLUA_PLUGINS was added we may need the path where to find liblua*.so of libluajit*.so
    link_directories(/usr/local/lib)
    #link_directories(/usr/local/lua-5.3.0/lib)
    #link_directories(/usr/lib64/libluajit-5.1.so)
    #link_directories(/usr/lib64)
endif()

add_executable(${TGT} ${SRCS})

if ("${CMAKE_C_COMPILER_ID}" MATCHES "Clang")
    # clang compiler : should perform well on all cpus
    set(CMAKE_C_FLAGS_RELEASE "-std=c99 -Wall -Wextra -msse2 -fno-strict-aliasing -O3")
    set(CMAKE_C_FLAGS_DEBUG   "-std=c99 -Wall -Wextra -msse2 -fno-strict-aliasing -O0 -g")

    # clang compiler : optimise for current machine
    #set(CMAKE_C_FLAGS_RELEASE "-std=c99 -Wall -Wextra -msse2 -fno-strict-aliasing -O3 -march=native")
    #set(CMAKE_C_FLAGS_DEBUG   "-std=c99 -Wall -Wextra -msse2 -fno-strict-aliasing -O0 -g")

	# clang
	target_link_libraries(${TGT} m)

	if (USE_LUA)
        # clang with lua
        #target_link_libraries(${TGT} m lua-5.3)
        
        # clang with luajit
        target_link_libraries(${TGT} m luajit-5.1)
    endif()

elseif("${CMAKE_C_COMPILER_ID}" MATCHES "GNU")
    # gnu compiler : should perform well on all cpus
    set(CMAKE_C_FLAGS_RELEASE "-std=c99 -Wall -Wextra -msse2 -fno-strict-aliasing -O3")
    set(CMAKE_C_FLAGS_DEBUG   "-std=c99 -Wall -Wextra -msse2 -fno-strict-aliasing -O0 -g")

    # gnu compiler : optimise for current machine
    #set(CMAKE_C_FLAGS_RELEASE "-std=c99 -Wall -Wextra -msse2 -fno-strict-aliasing -O3 -march=native")
    #set(CMAKE_C_FLAGS_DEBUG   "-std=c99 -Wall -Wextra -msse2 -fno-strict-aliasing -O0 -g")

	# gcc
	target_link_libraries(${TGT} m)

	if (USE_LUA)
        # gcc with lua
        #target_link_libraries(${TGT} m lua)
        #target_link_libraries(${TGT} m lua-5.2)
        
        # gcc with luajit
        target_link_libraries(${TGT} m luajit-5.1)
	endif()
	
elseif ("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
    # when using intel compiler : optimise for modern intel cpus
    set(CMAKE_C_FLAGS_RELEASE "-std=c99 -Wall -Wextra -msse2 -ax,CORE-AVX2,CORE-AVX-I,AVX,SSE4.2,SSE4.1,SSSE3,SSE3 -fno-strict-aliasing -O2")
    set(CMAKE_C_FLAGS_DEBUG   "-std=c99 -Wall -Wextra -msse2 -ax,CORE-AVX2,CORE-AVX-I,AVX,SSE4.2,SSE4.1,SSSE3,SSE3 -fno-strict-aliasing -O0 -g")

    # when using intel compiler : optimise for current machine
    #set(CMAKE_C_FLAGS_DEBUG   "-std=c99 -Wall -Wextra -fno-strict-aliasing -O0 -g")
    #set(CMAKE_C_FLAGS_RELEASE "-std=c99 -Wall -Wextra -fno-strict-aliasing -O3 -xHost")
    
    # icc
    target_link_libraries(${TGT} m)
    
    
    if (USE_LUA)
        # icc with luajit
        target_link_libraries(${TGT} m luajit-5.1)

        # icc with lua
        #target_link_libraries(${TGT} m lua)
        #target_link_libraries(${TGT} m lua-5.2)
    endif()
    
endif()
