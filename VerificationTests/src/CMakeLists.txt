# include headers
include_directories(${Engine_SOURCE_DIR}/CoreUtils/include)
include_directories(${Engine_SOURCE_DIR}/Eigen)
include_directories(${Engine_SOURCE_DIR}/Core/include)

link_directories(${Engine_SOURCE_DIR}/CoreUtils/src)
add_executable(Verifications verifications.cpp)
target_link_libraries(Verifications CoreUtils Core)
