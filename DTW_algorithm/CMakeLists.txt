cmake_minimum_required(VERSION 3.12)

project(DTW LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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

add_library(helper_tools SHARED helper_tools.cpp)
target_link_libraries(helper_tools PRIVATE ${ROOT_LIBRARIES})

add_executable(DTW.exe DTW.cpp ../Accidental_Events_Monte_Carlo/testDTW.cpp)
target_include_directories(DTW.exe PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(DTW.exe PRIVATE helper_tools ${ROOT_LIBRARIES})

