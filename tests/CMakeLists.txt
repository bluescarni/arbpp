# NOTE: a bit hackish, as in principle Boost could
# be compiled statically. Proper way would be to check
# this with a try_compile or similar.
add_definitions(-DBOOST_TEST_DYN_LINK)

macro(ADD_ARBPP_TESTCASE arg1)
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        add_executable(${arg1} ${arg1}.cpp)
        target_link_libraries(${arg1} ${MANDATORY_LIBRARIES})
        add_test(${arg1} ${arg1})
    endif()
endmacro(ADD_ARBPP_TESTCASE)

macro(ADD_ARBPP_PERFORMANCE_TESTCASE arg1)
    if(CMAKE_BUILD_TYPE STREQUAL "Release")
        add_executable("${arg1}_perf" ${arg1}.cpp)
        target_link_libraries("${arg1}_perf" ${MANDATORY_LIBRARIES})
        add_test("${arg1}_perf" "${arg1}_perf")
    endif()
endmacro(ADD_ARBPP_PERFORMANCE_TESTCASE)

ADD_ARBPP_TESTCASE(arb)
