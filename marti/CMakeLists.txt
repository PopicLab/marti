cmake_minimum_required(VERSION 3.20)
project(marti CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_compile_options(-Wno-error)
option(MARTI_BUILD_TESTS "Build tests" OFF)

# addresses third-party library cmake warnings
set(CMAKE_WARN_DEPRECATED OFF CACHE BOOL "" FORCE)
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0063 NEW)

######################
####### OpenMP #######
######################
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

########################
####### YAML-CPP #######
########################
include(FetchContent)
message(STATUS "Fetching yaml-cpp...")
FetchContent_Declare(
        yaml-cpp
        GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
        GIT_TAG 0.8.0
)
set(YAML_CPP_BUILD_TESTS OFF)
FetchContent_MakeAvailable(yaml-cpp)

#####################
####### edlib #######
#####################
message(STATUS "Fetching edlib...")
FetchContent_Declare(
        edlib
        GIT_REPOSITORY https://github.com/Martinsos/edlib.git
        GIT_TAG v1.2.7
)
set(BUILD_TESTING OFF)
set(EDLIB_BUILD_EXAMPLES OFF)
set(EDLIB_BUILD_UTILITIES OFF)
FetchContent_MakeAvailable(edlib)

#####################
####### marti #######
#####################
include_directories("/usr/local/include/") # (htslib)
link_directories("/usr/local/lib/")
message(STATUS "Configuring marti...")
include_directories(include)
set(SOURCES src/config.cc src/classifier.cc src/search.cc
        src/classifier.cc)
add_library(${PROJECT_NAME}_lib STATIC ${SOURCES})
add_executable (${PROJECT_NAME} src/main.cc)
target_link_libraries(${PROJECT_NAME}_lib PUBLIC edlib::edlib yaml-cpp::yaml-cpp hts)
target_link_libraries (${PROJECT_NAME} PUBLIC ${PROJECT_NAME}_lib)
target_compile_features(${PROJECT_NAME}_lib PUBLIC cxx_std_20)
set_target_properties(${PROJECT_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

if(MARTI_BUILD_TESTS)
    message(STATUS "Building tests...")
    enable_testing()
    add_subdirectory(tests)
endif()