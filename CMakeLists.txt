cmake_minimum_required( VERSION 2.6 )
project( test-png )

#Get libpng parameters (don't know, how to make it work)
#execute_process( COMMAND libpng-config --libs OUTPUT_VARIABLE PNG_LIBS )
set( PNG_LIBS "png12" )
set(CMAKE_CXX_FLAGS "-Wall -o3")        ## Optimize

find_package( Boost 1.34 COMPONENTS program_options  REQUIRED )
link_directories ( ${Boost_LIBRARY_DIRS} )
include_directories ( ${Boost_INCLUDE_DIRS} )


add_executable( vizualizer main.cpp )
TARGET_LINK_LIBRARIES(vizualizer ${PNG_LIBS} ${Boost_LIBRARIES} )

