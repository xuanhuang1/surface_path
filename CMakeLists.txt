cmake_minimum_required(VERSION 3.5)
project(surface_path)

if (NOT WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
endif()

add_executable(surface_path
  surface_path.cpp)

