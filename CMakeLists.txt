# Setup project
cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(BaseProject)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# set build type variable in preprocessor
if(CMAKE_BUILD_TYPE)
    string(TOUPPER "${CMAKE_BUILD_TYPE}" _upper_build_type)
    add_compile_definitions(BUILD_${_upper_build_type}=1)
endif()

add_compile_definitions(HAVE_INLINE=1)

# add include dir
include_directories(${PROJECT_SOURCE_DIR}/include)

# add libraries
link_directories(${PROJECT_SOURCE_DIR}/lib)

# add sources and headers
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)

# add executable
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_executable(anim_test anim_test.cpp ${sources} ${headers})
    add_executable(supertriangle supertriangle.cpp ${sources} ${headers})
    add_executable(circumcentre circumcentre.cpp ${sources} ${headers})
    add_executable(test_edge_set edge_set.cpp ${sources} ${headers})
    add_executable(delaunay delaunay.cpp ${sources} ${headers})
    add_executable(check_turns check_turns.cpp ${sources} ${headers})
    add_executable(convex_hull convex_hull.cpp ${sources} ${headers})
endif()
add_executable(problem1 problem1.cpp ${sources} ${headers})
add_executable(problem2 problem2.cpp ${sources} ${headers})

# add libraries
find_package(GSL REQUIRED)

# add Boost libraries
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.71.0 COMPONENTS iostreams system filesystem container)
include_directories(${Boost_INCLUDE_DIRS})
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_link_libraries(anim_test GSL::gsl GSL::gslcblas ${Boost_LIBRARIES})
    target_link_libraries(supertriangle GSL::gsl GSL::gslcblas ${Boost_LIBRARIES})
    target_link_libraries(circumcentre GSL::gsl GSL::gslcblas ${Boost_LIBRARIES})
    target_link_libraries(test_edge_set GSL::gsl GSL::gslcblas ${Boost_LIBRARIES})
    target_link_libraries(delaunay GSL::gsl GSL::gslcblas ${Boost_LIBRARIES})
    target_link_libraries(check_turns GSL::gsl GSL::gslcblas ${Boost_LIBRARIES})
    target_link_libraries(convex_hull GSL::gsl GSL::gslcblas ${Boost_LIBRARIES})
endif()
target_link_libraries(problem1 GSL::gsl GSL::gslcblas ${Boost_LIBRARIES})
target_link_libraries(problem2 GSL::gsl GSL::gslcblas ${Boost_LIBRARIES})
