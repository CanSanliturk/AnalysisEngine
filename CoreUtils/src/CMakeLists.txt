# Notice name prefix of this variable, set by CMake according
# to value given with "project()" in the root CMakeLists.txt.
include_directories(${Engine_SOURCE_DIR}/CoreUtils/include)
add_library(CoreUtils UtilMethods.cpp Vector.cpp GeometryHelper.cpp XYData.cpp Shape.cpp)
