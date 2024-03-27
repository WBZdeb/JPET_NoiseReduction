cmake_minimum_required(VERSION 3.1...3.16)

project(load_data CXX) # using only C++

add_definitions(-std=c++11 -Wall -Wunused-parameter -Wextra )

find_package(ROOT REQUIRED COMPONENTS 
  Hist
  HistPainter
  Physics
  RIO
  Thread
  Gpad
  Tree
  Minuit
  Html
  Core
  RooFit
  RooFitCore
  TreePlayer
  MathMore
  )
include_directories(${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_DEFINITIONS})
link_directories(${ROOT_LIBRARY_DIRS})

add_executable(load_data.exe ${CMAKE_CURRENT_SOURCE_DIR}/load_data.cpp)
target_link_libraries(load_data.exe ${ROOT_LIBRARIES})
