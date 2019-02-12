cmake_minimum_required(VERSION 3.9)
project(plastic_adjust CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB headers src/*.h src/*.hpp)
file(GLOB sources src/*.cpp)

find_package(NLopt REQUIRED)

set(main_SRCS ${sources} ${headers})
add_executable(plastic_adjust ${main_SRCS})

configure_file(input.json ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

if (WIN32)
    add_custom_command( TARGET plastic_adjust POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy_if_different
    	"${CMAKE_CURRENT_BINARY_DIR}/input.json" $<TARGET_FILE_DIR:plastic_adjust>
    )
endif()
