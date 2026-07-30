# del3d's viewer. Built only when DEL3D_BUILD_VIEWER is ON.
#
# This is the one place in the tree with a third-party dependency. It is kept in
# its own target and off by default so that the library, the tests and the
# examples stay buildable with nothing but a C++17 compiler.
#
# Polyscope is looked for as an installed package first, so a vendored or
# system copy wins; only if that fails - and only if DEL3D_FETCH_POLYSCOPE is
# left ON - is it downloaded at configure time.
find_package(polyscope CONFIG QUIET)

if(NOT polyscope_FOUND)
    if(NOT DEL3D_FETCH_POLYSCOPE)
        message(FATAL_ERROR
            "DEL3D_BUILD_VIEWER is ON but Polyscope was not found and "
            "DEL3D_FETCH_POLYSCOPE is OFF. Install Polyscope, point "
            "CMAKE_PREFIX_PATH at it, or allow the download.")
    endif()
    message(STATUS "del3d: Polyscope not found, fetching ${DEL3D_POLYSCOPE_TAG}")
    include(FetchContent)
    FetchContent_Declare(polyscope
        GIT_REPOSITORY https://github.com/nmwsharp/polyscope.git
        GIT_TAG        ${DEL3D_POLYSCOPE_TAG}
        GIT_SHALLOW    TRUE
    )
    FetchContent_MakeAvailable(polyscope)
endif()

add_executable(del3d_viewer
    main.cpp
    mesh_io.cpp
    mesh_io.h
)
target_link_libraries(del3d_viewer PRIVATE del3d polyscope)

# --check drives the whole viewer against Polyscope's mock backend and exits
# without opening a window, so it is worth running wherever the tests run.
if(DEL3D_BUILD_TESTS)
    add_test(NAME del3d_viewer_check COMMAND del3d_viewer --check)
endif()

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/mesh_io.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/mesh_io.h
)
