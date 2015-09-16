
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/tools")

# link_libraries(${Boost_LIBRARIES} ${RG_LIBRARIES}  ta)
link_libraries(${Boost_LIBRARIES} ${RG_LIBRARIES}  ${GKLIB_LIBRARIES} ta )

add_executable(testAlgoTuning testAlgoTuning.cc)
add_executable(testTA_all testTA_all.cc)
add_executable(testTANRA_all testTANRA_all.cpp)
add_executable(testPCATree testPCATree.cpp)