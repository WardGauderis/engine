cmake_minimum_required(VERSION 3.13)
project(engine_cmake)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)
include_directories(l_parser)
set(CMAKE_CXX_FLAGS "-O2")
add_executable(engine_cmake
        l_parser/l_parser.cc
        l_parser/l_parser.h
        vector/vector3d.cc
        vector/vector3d.h
        easy_image.cc
        easy_image.h
        engine.cc
        ini_configuration.cc
        ini_configuration.h
        Line2D.cpp
        Line2D.h Figure.cpp Figure.h "ZBuffer.cpp" "ZBuffer.h")
