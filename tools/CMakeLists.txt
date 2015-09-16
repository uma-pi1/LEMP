
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/tools")


# link_libraries(${Boost_LIBRARIES} ${RG_LIBRARIES}  ta rg  boost_filesystem boost_system)
link_libraries( armadillo ${Boost_LIBRARIES} ${RG_LIBRARIES}  ta rg  boost_filesystem boost_system)

add_executable(runNaive runNaive.cc)
add_executable(runAlgoWithTuner runAlgoWithTuner.cc)
add_executable(runTA_all runTA_all.cc)
# add_executable(runLSH_all runLSH_all.cc)
add_executable(runPCATree runPCATree.cpp)