cmake_minimum_required(VERSION 3.5)
project(sltp-features VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


# Uncomment for debug build. Alternatively, build with `cmake -DCMAKE_BUILD_TYPE=Debug .`
#set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra")
set(CMAKE_CXX_FLAGS "-O3")


set(SOURCE_FILES
        src/main.cpp
        src/sltp/features
        src/sltp/utils.hxx
)

add_executable(featuregen ${SOURCE_FILES})
set_target_properties(featuregen PROPERTIES LINKER_LANGUAGE CXX)


set(PROJECT_ROOT .)
#MESSAGE(WARNING "Your FS_ROOT home directory is: ${PROJECT_ROOT}" )

include_directories(${PROJECT_ROOT}/src)
include_directories(${PROJECT_ROOT}/vendor)



#find_package(Boost 1.68 COMPONENTS program_options REQUIRED)
#include_directories(${Boost_INCLUDE_DIR})
#target_link_libraries(genfeatures LINK_PUBLIC ${Boost_LIBRARIES})