cmake_minimum_required(VERSION 3.9)
project(plastic_adjust CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB headers src/*.h src/*.hpp)
file(GLOB sources src/*.cpp)

set(CMAKE_CONFIGURATION_TYPES Debug Release CACHE TYPE  INTERNAL FORCE )

find_package(NLopt)

include_directories( ${NLopt_INCLUDE_DIR})

IF (WIN32)
   find_package( PZ PATHS "C:/Arquivos de Programas/PZ" REQUIRED NO_DEFAULT_PATH)
ELSE()
   find_package( PZ REQUIRED )
ENDIF (WIN32)

include_directories(${PZ_INCLUDE_DIRS})

if(WIN32)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_STRUCT_TIMESPEC")
   #add a define saying it is a VC compiler
   set(VC "#Will add a define VC on pz_config.h")
   add_definitions(-DVC)
   #define _USE_MATH_DEFINES for example: M_PI constant definitions.
   add_definitions(-D_USE_MATH_DEFINES)
   #disabling VC warnings
   add_definitions(-D_SCL_SECURE_NO_WARNINGS)
   add_definitions(-D_CRT_SECURE_NO_WARNINGS)

   #define use of pthread static lib.
   add_definitions(-DPTW32_STATIC_LIB)
   #enabling /bigobj
   add_definitions("/bigobj")

   FOREACH(FLAG_TYPE EXE MODULE SHARED)
       STRING (REPLACE "INCREMENTAL:YES" "INCREMENTAL" FLAG_TMP "${CMAKE_${FLAG_TYPE}_LINKER_FLAGS_DEBUG}")
       STRING (REPLACE "INCREMENTAL:NO" "INCREMENTAL" FLAG_TMP ${FLAG_TMP})
       STRING (REPLACE "INCREMENTAL" "INCREMENTAL:NO" FLAG_TMP ${FLAG_TMP})
       STRING (REPLACE "/EDITANDCONTINUE" "" FLAG_TMP ${FLAG_TMP})
       SET(CMAKE_${FLAG_TYPE}_LINKER_FLAGS_DEBUG "${FLAG_TMP}" CACHE STRING "Overriding default debug ${FLAG_TYPE} linker flags." FORCE)
       MARK_AS_ADVANCED (CMAKE_${FLAG_TYPE}_LINKER_FLAGS_DEBUG)
   ENDFOREACH ()
else()
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-narrowing")#this flag is for preventing errors when casting from double to float
endif()


set(main_SRCS ${sources} ${headers})
add_executable(plastic_adjust ${main_SRCS})

target_link_libraries(plastic_adjust ${NLopt_LIB} pz)

configure_file(input.json ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(input2.json ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
file(COPY Inputdata DESTINATION .)

if (WIN32 OR APPLE)
    add_custom_command( TARGET plastic_adjust POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy_if_different
    	"${CMAKE_CURRENT_BINARY_DIR}/input.json" $<TARGET_FILE_DIR:plastic_adjust>
    )
endif()


