cmake_minimum_required( VERSION 2.6 )
project( test-png )

#enable_testing()

#Get lib[ng parameters
#execute_process( COMMAND libpng-config --libs OUTPUT_VARIABLE PNG_LIBS )
set( PNG_LIBS "png12" )
set(CMAKE_CXX_FLAGS "-g -Wall -o3")        ## Optimize

#include_directories( include )
#include_directories( ext/gtest/include )
#FIND_PACKAGE(png++)

add_executable( main main.cpp )
TARGET_LINK_LIBRARIES(main ${PNG_LIBS} )

