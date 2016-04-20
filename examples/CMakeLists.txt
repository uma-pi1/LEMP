
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/examples")

# link_libraries(${Boost_LIBRARIES} ${RG_LIBRARIES}  ta)
link_libraries(${Boost_LIBRARIES} ${RG_LIBRARIES}  ${GKLIB_LIBRARIES} ta )

add_executable(testLemp testLemp.cc)
add_executable(testTa testTa.cc)
add_executable(testTaNra testTaNra.cpp)
add_executable(testPcaTree testPcaTree.cpp)
add_executable(testSimpleLsh testSimpleLsh.cc)
add_executable(testNaive testNaive.cpp)