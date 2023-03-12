cmake_minimum_required(VERSION 3.10)

set (ARMATUS_MAIN Armatus.cpp ArmatusUtil.cpp ArmatusParams.cpp ArmatusDAG.cpp IntervalScheduling.cpp)

set(MPI_EXECUTABLE_SUFFIX ".mpich")
find_package(MPI)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

include_directories( 
${Boost_INCLUDE_DIR}
${TOP_SOURCE_DIR}/include
)

link_directories(${Boost_LIBRARY_DIR})

add_executable(armatus ${ARMATUS_MAIN})
target_link_libraries(armatus ${Boost_LIBRARIES} z)
target_link_libraries(armatus ${MPI_C_LIBRARIES} z)


install (TARGETS armatus DESTINATION bin)
