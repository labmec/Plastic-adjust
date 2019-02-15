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
 
elseif(MACOS)
   SET(NLopt_DIR "/usr/local/Cellar/nlopt" CACHE PATH "Root directory where the NLopt library and headers are istalled")
   message("NLopt in MacOS")
   
else() 
   SET(NLopt_DIR "/usr" CACHE PATH "Root directory where the NLopt library and headers are  istalled")
   message("NLopt in Linux")
endif()

if (EXISTS ${NLopt_DIR})
    find_path( NLopt_INCLUDE_DIR  NAMES nlopt.h
                HINTS "${NLopt_DIR}/include" "${NLopt_DIR}/*/include" 
    )

    #set ( NLopt_INCLUDE_DIR ${NLopt_INCLUDE_DIR} )
  
    find_library(NLopt_LIB  NAMES nlopt_cxx.lib libnlopt.so libnlopt_cxx.dylib 
	              HINTS "${NLopt_DIR}/lib" "${NLopt_DIR}/*/lib" NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)

	message("NLopt_LIB: ${NLopt_LIB}")			  
				  
    if (${NLopt_LIB-NOTFOUND})
        find_library(NLopt_LIB  NAMES nlopt_cxx.lib libnlopt.so libnlopt_cxx.dylib PATHS ${SEARCH_DIRS})
    endif()

endif()


mark_as_advanced (
  NLopt_LIB
  NLopt_INCLUDE_DIR
)

