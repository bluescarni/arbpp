project(arbpp)

enable_testing()

cmake_minimum_required(VERSION 2.8.10)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules")

message(STATUS "System name: ${CMAKE_SYSTEM_NAME}")

# Set default build type to "Release".
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
endif()

# Build option: enable test set.
option(BUILD_TESTS "Build test set." OFF)

# Compiler setup.
if(BUILD_TESTS)
    # Setup compiler.
    include(ArbppCompilerLinkerSettings)
    # GMP.
    find_package(GMP REQUIRED)
    message(STATUS "GMP library found.")
    message(STATUS "GMP include dir is: ${GMP_INCLUDE_DIR}")
    message(STATUS "GMP library is: ${GMP_LIBRARIES}")
    include_directories(${GMP_INCLUDE_DIR})
    # FLINT.
    find_package(FLINT REQUIRED)
    message(STATUS "FLINT library found.")
    message(STATUS "FLINT include dir is: ${FLINT_INCLUDE_DIR}")
    message(STATUS "FLINT library is: ${FLINT_LIBRARIES}")
    # NOTE: arb looks for flint.h without prefix, but we look for
    # the include dir as $INCLUDE_DIR/flint/flint.h. Include both paths.
    include_directories(${FLINT_INCLUDE_DIR})
    include_directories(${FLINT_INCLUDE_DIR}/flint)
    # MPFR.
    find_package(MPFR REQUIRED)
    message(STATUS "MPFR library found.")
    message(STATUS "MPFR include dir is: ${MPFR_INCLUDE_DIR}")
    message(STATUS "MPFR library is: ${MPFR_LIBRARIES}")
    include_directories(${MPFR_INCLUDE_DIR})
    # Arb.
    find_package(Arb REQUIRED)
    message(STATUS "Arb library found.")
    message(STATUS "Arb include dir is: ${Arb_INCLUDE_DIR}")
    message(STATUS "Arb library is: ${Arb_LIBRARIES}")
    include_directories(${Arb_INCLUDE_DIR})
    # Boost unit test library.
    find_package(Boost 1.48.0 REQUIRED COMPONENTS "unit_test_framework")
    include_directories(${Boost_INCLUDE_DIRS})
    # Assemble all libraries and add the tests subdirectory.
    set(MANDATORY_LIBRARIES ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
        ${Arb_LIBRARIES}
        ${MPFR_LIBRARIES}
        ${FLINT_LIBRARIES}
        ${GMP_LIBRARIES}
    )
    add_subdirectory("${CMAKE_SOURCE_DIR}/tests")
endif()

# Install the headers.
install(FILES src/arbpp.hpp DESTINATION include/arbpp)
