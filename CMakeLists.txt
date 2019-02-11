cmake_minimum_required(VERSION 3.9)
project(plastic_adjust CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB headers src/*.h src/*.hpp)
file(GLOB sources src/*.cpp)

set(main_SRCS ${sources} ${headers})
add_executable(plastic_adjust ${main_SRCS})

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
configure_file(input.json ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
