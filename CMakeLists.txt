cmake_minimum_required(VERSION 3.27)
project(bezier_simplification)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ~/Documents/cartocrow-bins/)

include_directories("~/Documents/cartocrow-bins/include")
find_package(CartoCrow REQUIRED)
find_package(spiro)

set(CMAKE_AUTOMOC ON)

add_subdirectory(library)
add_subdirectory(demos)
#add_subdirectory(test)
