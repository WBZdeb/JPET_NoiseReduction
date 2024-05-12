cmake_minimum_required(VERSION 3.1...3.16)

project(DTW LANGUAGES CXX)

add_definitions(-std=c++11 -Wall -Wunused-parameter -Wextra  -g)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

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
  ROOTDataFrame
  MultiProc
  ROOTVecOps 
  )

include_directories(${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_DEFINITIONS})
link_directories(${ROOT_LIBRARY_DIRS})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(helper_tools SHARED helper_tools.cpp ${ROOT_LIBRARIES})

add_executable(DTW.exe DTW.cpp)
target_link_libraries(DTW.exe PRIVATE helper_tools ${ROOT_LIBRARIES})
