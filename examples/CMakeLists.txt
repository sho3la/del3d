# del3d examples. Built only when DEL3D_BUILD_EXAMPLES is ON.
#
# Everything here uses the public interface and nothing else, so an example is
# also a check that include/ is self-sufficient: it links del3d without ever
# adding src/ to the include path.
add_executable(del3d_example_simple simple.cpp)
target_link_libraries(del3d_example_simple PRIVATE del3d)

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${CMAKE_CURRENT_SOURCE_DIR}/simple.cpp)
