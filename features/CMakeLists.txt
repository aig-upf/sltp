cmake_minimum_required(VERSION 3.5)
project(sltp-features VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)




set(SOURCE_FILES
        src/sltp/features
        src/main.cpp src/sltp/io.hxx src/sltp/io.cpp)

add_executable(featuregen ${SOURCE_FILES})
set_target_properties(featuregen PROPERTIES LINKER_LANGUAGE CXX)


set(PROJECT_ROOT .)
#MESSAGE(WARNING "Your FS_ROOT home directory is: ${PROJECT_ROOT}" )

include_directories(${PROJECT_ROOT}/src)
include_directories(${PROJECT_ROOT}/vendor)



#find_package(Boost 1.68 COMPONENTS program_options REQUIRED)
#include_directories(${Boost_INCLUDE_DIR})
#target_link_libraries(genfeatures LINK_PUBLIC ${Boost_LIBRARIES})