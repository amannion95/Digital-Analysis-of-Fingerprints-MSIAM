cmake_minimum_required(VERSION 2.8)

project( src )

find_package( OpenCV REQUIRED )

include_directories( ${OpenCV_INCLUDE_DIRS} )

add_executable( boximage main.cpp Picture.cpp Usefull_functions.cpp )
add_executable( DisplayImage DisplayImage.cpp )

target_link_libraries( boximage ${OpenCV_LIBS} )
target_link_libraries( DisplayImage ${OpenCV_LIBS} )
