cmake_minimum_required(VERSION 2.8)
project(DirectedDensestSubgraph)

set(CMAKE_CXX_STANDARD 11)

find_package(Boost 1.66.0)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
endif()

add_executable(DirectedDensestSubgraph main.cpp Graph.cpp Graph.h Args.cpp Args.h FlowNetwork.cpp FlowNetwork.h EdgeFN.cpp EdgeFN.h)

target_link_libraries(DirectedDensestSubgraph ${Boost_LIBRARIES})
