cmake_minimum_required (VERSION 3.10)

project(CWDT)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_BUILD_TYPE "Release")
add_compile_definitions(CWDT_VERBOSE) # for test verbose
add_compile_definitions(CWDT_DEBUG) # for test output
add_compile_definitions(OUTPUT_PATH="${PROJECT_SOURCE_DIR}/data/") # set output file path
# add_compile_options(/openmp:llvm)

# CMake module path
list(PREPEND CMAKE_MODULE_PATH
        ${CMAKE_CURRENT_LIST_DIR}/cmake)

# Eigen
find_package(Eigen3 REQUIRED)
IF(Eigen3_FOUND)
    message(STATUS "Eigen3 FOUNDED ${EIGEN3_VERSION_STRING} IN ${EIGEN3_INCLUDE_DIR}")
    include_directories(${EIGEN3_INCLUDE_DIR})
ELSE()
    message(FATAL_ERROR "Eigen3 NOT FOUNDED")
ENDIF()

# Boost
find_package(Boost REQUIRED)
IF(Boost_FOUND)
    message(STATUS "Boost FOUNDED ${Boost_VERSION_STRING} IN ${Boost_INCLUDE_DIRS}")
    include_directories(${Boost_INCLUDE_DIRS})
ELSE()
    message(FATAL_ERROR "Boost NOT FOUNDED")
ENDIF()

# CGAL
find_package(CGAL 5.6 REQUIRED)
IF(CGAL_FOUND)
    message(STATUS "CGAL FOUNDED ${CGAL_VERSION_STRING} IN ${CGAL_INCLUDE_DIRS}")
    include_directories(${CGAL_INCLUDE_DIRS})
ELSE()
    message(FATAL_ERROR "CGAL NOT FOUNDED")
ENDIF()

# Gurobi
find_package(Gurobi REQUIRED)
IF(GUROBI_FOUND)
    message(STATUS "Gurobi FOUNDED IN ${GUROBI_INCLUDE_DIR}")
    include_directories(${GUROBI_INCLUDE_DIR})
ELSE()
    message(FATAL_ERROR "Gurobi NOT FOUNDED")
ENDIF()

# libigl
find_package(libigl REQUIRED)
IF(LIBIGL_FOUND)
    message(STATUS "libigl FOUNDED IN ${LIBIGL_INCLUDE_DIRS}")
    include_directories(${LIBIGL_INCLUDE_DIRS})
ELSE()
    message(FATAL_ERROR "libigl NOT FOUNDED")
ENDIF()


include_directories(${PROJECT_SOURCE_DIR}/src)

add_executable(CWDT
        "${CMAKE_SOURCE_DIR}/src/cwdt/cwdt.cpp"
        "${CMAKE_SOURCE_DIR}/src/main.cpp")
target_link_libraries(CWDT ${GUROBI_LIBRARY})
target_link_libraries(CWDT ${GUROBI_CPP_LIBRARY})