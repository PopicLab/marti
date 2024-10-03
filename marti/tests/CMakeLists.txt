include(FetchContent)
FetchContent_Declare(
        doctest
        GIT_REPOSITORY https://github.com/onqtam/doctest.git
        GIT_TAG v2.4.11
)
FetchContent_MakeAvailable(doctest)
include_directories(${edlib_SOURCE_DIR}/doctest)
enable_testing()

###################
# test_classifier #
###################

add_executable(test_classifier test_classifier.cc simulator.cc)
target_compile_features(test_classifier PRIVATE cxx_std_20)
target_link_libraries(test_classifier PRIVATE ${PROJECT_NAME}_lib doctest::doctest)
set_target_properties(test_classifier PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
add_test(NAME test_classifier COMMAND test_classifier)
set_property(TEST test_classifier PROPERTY ENVIRONMENT MARTI_TESTS_ROOT=${CMAKE_SOURCE_DIR})

###############
# test_config #
###############

add_executable(test_config test_config.cc)
target_compile_features(test_config PRIVATE cxx_std_20)
target_link_libraries(test_config PRIVATE ${PROJECT_NAME}_lib doctest::doctest)
set_target_properties(test_config PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
add_test(NAME test_config COMMAND test_config)
set_property(TEST test_config PROPERTY ENVIRONMENT MARTI_TESTS_ROOT=${CMAKE_SOURCE_DIR})

###############
# test_search #
###############

add_executable(test_search test_search.cc)
target_compile_features(test_search PRIVATE cxx_std_20)
target_link_libraries(test_search PRIVATE ${PROJECT_NAME}_lib doctest::doctest)
set_target_properties(test_search PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
add_test(NAME test_search COMMAND test_search)
set_property(TEST test_search PROPERTY ENVIRONMENT MARTI_TESTS_ROOT=${CMAKE_SOURCE_DIR})
