PROJECT(geotagger CXX)
set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -g")
cmake_minimum_required(VERSION 2.6)
find_package(Exiv2)

find_package(Boost COMPONENTS filesystem system REQUIRED)
include_directories(${EXIV2_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})

ADD_EXECUTABLE(geotagger main.cpp tinyxml2.cpp)
TARGET_LINK_LIBRARIES(geotagger ${EXIV2_LIBRARIES} ${Boost_LIBRARIES})
