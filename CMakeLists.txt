# Select minimum version for cmake
cmake_minimum_required(VERSION 3.15)

# Set c++ standard to be used
set (CMAKE_CXX_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# Name the project
project(Engine)

# Add subdirectories
add_subdirectory(CoreUtils)
add_subdirectory(CoreUtilsTests)
add_subdirectory(Core)
add_subdirectory(VerificationTests)
