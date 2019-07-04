cmake_minimum_required(VERSION 3.4)
project(contorno VERSION 1.0.0 LANGUAGES CXX )

set(CMAKE_CXX_COMPILER g++)

find_package( OpenCV REQUIRED )
include_directories( ${OpenCV_INCLUDE_DIRS} )


add_executable(contorno contorno.cpp)
add_executable(blur blur.cpp)
target_link_libraries(contorno ${OpenCV_LIBS})
target_link_libraries(blur ${OpenCV_LIBS})
# target_link_libraries(contorno OpenCV::opencv)
target_compile_options(contorno PRIVATE -Wall)
target_compile_options(blur PRIVATE -Wall)

install( TARGETS contorno blur DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin COMPONENT binaries)
