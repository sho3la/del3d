cmake_minimum_required(VERSION 3.16)
project(del3d VERSION 0.1 LANGUAGES CXX)

# del3d - 3D Delaunay triangulation.
#
# The library itself uses nothing but the C++17 standard library, and is built
# as a plain static library so it can be dropped into any project. So do the
# tests and the examples. The viewer is the exception - it needs Polyscope, and
# is therefore a separate target that is off by default.

option(DEL3D_BUILD_TESTS    "Build del3d's tests"                        OFF)
option(DEL3D_BUILD_EXAMPLES "Build del3d's examples"                     OFF)
option(DEL3D_BUILD_VIEWER   "Build the Polyscope viewer (needs OpenGL)"  OFF)
option(DEL3D_FETCH_POLYSCOPE
       "If the viewer is on and Polyscope is not installed, download it"  ON)
set(DEL3D_POLYSCOPE_TAG "v2.6.1" CACHE STRING
    "The Polyscope tag to fetch when DEL3D_FETCH_POLYSCOPE is used")

# The headers are listed as sources so that IDE project generators - the Visual
# Studio solution explorer, Xcode, and the file-api that VS Code and CLion read
# - show them under the target. CMake does not compile them; source_group()
# below mirrors the directory layout instead of dumping them into one folder.
set(DEL3D_PUBLIC_HEADERS
    include/del3d/delaunay.h
    include/del3d/predicates.h
    include/del3d/predicates2.h
    include/del3d/spatial_sort.h
)

set(DEL3D_PRIVATE_HEADERS
    src/compact_container.h
    src/expansion.h
    src/interval.h
    src/rand48.h
    src/tds.h
)

set(DEL3D_SOURCES
    src/predicates.cpp
    src/predicates2.cpp
    src/spatial_sort.cpp
    src/triangulation.cpp
)

add_library(del3d STATIC
    ${DEL3D_SOURCES}
    ${DEL3D_PUBLIC_HEADERS}
    ${DEL3D_PRIVATE_HEADERS}
)

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR}
    FILES ${DEL3D_SOURCES} ${DEL3D_PUBLIC_HEADERS} ${DEL3D_PRIVATE_HEADERS}
)

target_include_directories(del3d PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_include_directories(del3d PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

target_compile_features(del3d PUBLIC cxx_std_17)

# The predicates are exact only under strict IEEE-754 semantics. The error-free
# transformations in src/expansion.h recover the round-off of each operation by
# recomputing it from the rounded result, which requires that every operation be
# correctly rounded and that none of them be reassociated, contracted into a
# fused multiply-add, or evaluated in a wider register than the type asks for.
#
# src/interval.h additionally changes the FPU rounding mode at run time, so the
# compiler must not assume round-to-nearest when folding constants or hoisting
# arithmetic across that change.
#
# Building with -ffast-math (or /fp:fast) breaks both. Do not.
if(MSVC)
    target_compile_options(del3d PRIVATE /fp:strict)
else()
    target_compile_options(del3d PRIVATE -ffp-contract=off -frounding-math)
    # 32-bit x86 evaluates in 80-bit x87 registers unless told otherwise, which
    # keeps excess precision and invalidates the error-free transformations.
    if(CMAKE_SIZEOF_VOID_P EQUAL 4)
        target_compile_options(del3d PRIVATE -msse2 -mfpmath=sse)
    endif()
endif()

if(DEL3D_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

if(DEL3D_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

if(DEL3D_BUILD_VIEWER)
    add_subdirectory(viewer)
endif()
