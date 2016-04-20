
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/tools")


# link_libraries(${Boost_LIBRARIES} ${RG_LIBRARIES}  ta rg  boost_filesystem boost_system)
link_libraries( armadillo ${Boost_LIBRARIES} ${RG_LIBRARIES}  ta rg  boost_filesystem boost_system)

add_executable(runNaive runNaive.cc)
add_executable(runLemp runLemp.cc)
add_executable(runTa runTa.cc)
add_executable(runSimpleLsh runSimpleLsh.cc)
add_executable(runPcaTree runPcaTree.cpp)
add_executable(compare compare.cc)