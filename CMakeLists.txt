cmake_minimum_required (VERSION 2.6)
project (SBNfit)

# The version number.
set (SBNfit_VERSION_MAJOR 3)
set (SBNfit_VERSION_MINOR 0)

if (CMAKE_VERSION VERSION_LESS "3.1")
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
          set (CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
    endif ()
    else ()
	    set (CMAKE_CXX_STANDARD 11)
endif ()

message("@main CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
message("")
unset (CMAKE_CXX_FLAGS)

option(USE_GPU "USE GPU" OFF)
if(USE_GPU)
	message("@main Using GPU")
	message("")
	add_compile_definitions(USE_GPU)
else()
	message("@main Not using GPU")
	message("")
endif()


FIND_PACKAGE(ROOT)

if(ROOT_FOUND)
	message("@main Found ROOT using CMAKE FIND_PACKAGE")
	message("")
	include(${ROOT_USE_FILE})
ELSE()
	message("@main ROOT not compiled using CMAKE, use root-config")
	message("")
	FIND_PROGRAM(ROOT_CONFIG root-config)
	EXEC_PROGRAM(${ROOT_CONFIG} ARGS --cflags --glibs OUTPUT_VARIABLE ROOT_CXX_FLAGS)
	message("@main ROOT cxx flags set to: ${ROOT_CXX_FLAGS}")
	message("")
ENDIF()

# Get our GSL libraries
FIND_PROGRAM(GSL_CONFIG gsl-config)
EXEC_PROGRAM(${GSL_CONFIG} ARGS --cflags OUTPUT_VARIABLE GSL_CXX_FLAGS)
EXEC_PROGRAM(${GSL_CONFIG} ARGS --libs OUTPUT_VARIABLE GSL_LIBS)

message("@main CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
message("")
unset (CMAKE_CXX_FLAGS)

message("@main GSL_LIBS: ${GSL_LIBS}")
message("")
message("@main GSL_CXX_FLAGS: ${GSL_CXX_FLAGS}")
message("")

if(USE_GPU)
    set(CMAKE_GPU_CXX_FLAGS "-noswitcherror -ta=tesla:nollvm -Mcudalib=curand -Minfo=accel -acc -w" )
ENDIF()

#gopt is forcing LLVM compilation
# -noswitcherror only for pgcc

set (CMAKE_SO_CXX_FLAGS "${ROOT_CXX_FLAGS} ${GSL_LIBS} ${GSL_CXX_FLAGS} ${CMAKE_CXX_FLAGS} -w -Wno-deprecated-declarations")
if(USE_GPU)
	message("@main Add -noswitcherror to CMAKE_SO_CXX_FLAGS")
	message("")
	set (CMAKE_SO_CXX_FLAGS "-noswitcherror ${ROOT_CXX_FLAGS} ${GSL_LIBS} ${CMAKE_CXX_FLAGS} -w -Wno-deprecated-declarations")
ENDIF()

set (CMAKE_ACC_CXX_FLAGS "${CMAKE_GPU_CXX_FLAGS} ${CMAKE_SO_CXX_FLAGS}")
set (CMAKE_PCC_CXX_FLAGS "-ta=multicore -Minfo=accel -acc -w ${CMAKE_SO_CXX_FLAGS}")
set (CMAKE_OMP_CXX_FLAGS "-mp -acc -w ${CMAKE_SO_CXX_FLAGS}")

set (CMAKE_CXX_FLAGS "${CMAKE_ACC_CXX_FLAGS}")

message("@main CMAKE_SO_CXX_FLAGS: ${CMAKE_SO_CXX_FLAGS}")
message("@main CMAKE_ACC_CXX_FLAGS: ${CMAKE_ACC_CXX_FLAGS}")
message("@main CMAKE_ACC_PXX_FLAGS: ${CMAKE_PCC_CXX_FLAGS}")
message("@main CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")

include_directories(${PROJECT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/inc/)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/tinyxml/)

add_subdirectory(libio)

unset(CMAKE_CXX_FLAGS)
set(CMAKE_CXX_FLAGS "${CMAKE_ACC_CXX_FLAGS}")
message("@main CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
message("")

add_subdirectory(tinyxml)
add_subdirectory(src)
add_subdirectory(examples)

unset(CMAKE_CXX_FLAGS)
set(CMAKE_CXX_FLAGS "${CMAKE_OMP_CXX_FLAGS}")
message("@main CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
add_subdirectory(lee)



