## LabMeC, Feb/2019
## Assuming following installation types:
#
# Windows: download tgz with source, configuration and generation via cmake (must generate file nlopt_cxx.lib) and installation
# Linux: sudo apt-get install libnlopt-dev (TODO: check if really not necessary a nlopt_cxx.lib file, as there is none in this package (maybe generate Windows likewise)
# Mac OS: brew install nlopt
#
# Resulting variables:
#
# NLopt_INCLUDE_DIR contains path to header files of NLopt library
#
# NLopt_LIB contains path to object files of NLopt library

if(WIN32)
   SET(NLopt_DIR "C:/Arquivos de Programas/nlopt" CACHE PATH "Root directory where the NLopt library and headers are installed")
elseif(APPLE)
   SET(NLopt_DIR "/usr/local/Cellar/nlopt" CACHE PATH "Root directory where the NLopt library and headers are istalled")
else() 
   SET(NLopt_DIR "/usr" CACHE PATH "Root directory where the NLopt library and headers are  istalled")
endif()

if (EXISTS ${NLopt_DIR})
    find_path( NLopt_INCLUDE_DIR  NAMES nlopt.h
                HINTS "${NLopt_DIR}/include" "${NLopt_DIR}/*/include" 
    )
  
    find_library(NLopt_LIB  NAMES nlopt_cxx.lib nlopt.lib libnlopt.so libnlopt_cxx.dylib libnlopt.dylib 
	              HINTS "${NLopt_DIR}/lib" "${NLopt_DIR}/*/lib" "${NLopt_DIR}/lib/x86_64-linux-gnu" NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)

				  
    if (${NLopt_LIB-NOTFOUND})
        find_library(NLopt_LIB  NAMES nlopt_cxx.lib nlopt.lib libnlopt.so libnlopt_cxx.dylib libnlopt.dylib PATHS ${SEARCH_DIRS})
    endif()

endif()


mark_as_advanced (
  NLopt_LIB
  NLopt_INCLUDE_DIR
)

