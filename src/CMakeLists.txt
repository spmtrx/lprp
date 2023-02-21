include(GNUInstallDirs)

#include(FetchContent)
#FetchContent_Declare(
	#spmtrx
	#GIT_REPOSITORY https:://github.com/spmtrx/spmtrx.git
	#GIT_TAG tags/v0.1
	#GIT_SHALLOW TRUE)
#FetchContent_MakeAvailable(spmtrx)
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

message(STATUS "BOOST_ROOT=${BOOST_ROOT}")
find_package(Boost COMPONENTS serialization REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LINK_DIR})

include_directories(${CMAKE_SOURCE_DIR}/include)
add_library(lprp SHARED lprp.cc)
find_package(OpenMP)
if (OPENMP_FOUND)
	message(STATUS "OpenMP_CXX_FOUND : ${OpenMP_CXX_FOUND}")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

install(TARGETS lprp DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/lprp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} FILES_MATCHING PATTERN "*.h")
