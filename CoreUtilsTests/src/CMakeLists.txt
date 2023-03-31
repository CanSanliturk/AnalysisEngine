# include headers
include_directories(${Engine_SOURCE_DIR}/CoreUtils/include)
include_directories(${Engine_SOURCE_DIR}/CoreUtilsTests/include)

link_directories(${Engine_SOURCE_DIR}/CoreUtils/src)
add_executable(CoreUtilsTests tests.cpp)
target_link_libraries(CoreUtilsTests CoreUtils)
