# del3d tests.
#
# Self-contained, like the library: nothing here needs anything beyond the C++17
# standard library. The predicates are checked against tests/exact.h, an
# independent arbitrary-precision integer reference that shares none of the
# library's filtered arithmetic.

# Shared by every test but test_container, and listed on those targets so IDEs
# show them alongside the test they belong to.
set(DEL3D_TEST_HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/exact.h
    ${CMAKE_CURRENT_SOURCE_DIR}/pointsets.h
)

foreach(t container predicates predicates2 spatial_sort delaunay)
    if("${t}" STREQUAL "container")
        add_executable(del3d_test_${t} test_${t}.cpp)
    else()
        add_executable(del3d_test_${t} test_${t}.cpp ${DEL3D_TEST_HEADERS})
    endif()
    target_link_libraries(del3d_test_${t} PRIVATE del3d)
    target_include_directories(del3d_test_${t} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../src)
    source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR}
        FILES ${CMAKE_CURRENT_SOURCE_DIR}/test_${t}.cpp ${DEL3D_TEST_HEADERS})
    add_test(NAME del3d_${t} COMMAND del3d_test_${t})
endforeach()
