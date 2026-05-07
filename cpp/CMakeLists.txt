cmake_minimum_required(VERSION 3.20)
project(cpp_incremental_slam VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(ISLAM_BUILD_TESTS "Build tests" ON)
option(ISLAM_WITH_SUITESPARSE "Enable SuiteSparse / CHOLMOD when found" ON)
option(ISLAM_REQUIRE_SUITESPARSE_CCOLAMD "Require the SuiteSparse CCOLAMD direct API for byte-for-byte compatibility audits" OFF)
option(ISLAM_WITH_GTSAM "Build optional GTSAM Option-A IGG-SPO implementation" OFF)
option(ISLAM_REQUIRE_GTSAM "Fail configuration when ISLAM_WITH_GTSAM=ON and GTSAM is not found" OFF)
option(ISLAM_ALLOW_FETCHCONTENT "Allow CMake to fetch missing third-party dependencies" OFF)
set(ISLAM_EIGEN_INCLUDE_DIR "" CACHE PATH "Path to local Eigen headers for offline builds")

find_package(Eigen3 3.3 QUIET NO_MODULE)
if(NOT Eigen3_FOUND)
    if(ISLAM_EIGEN_INCLUDE_DIR AND EXISTS "${ISLAM_EIGEN_INCLUDE_DIR}/Eigen/Core")
        message(STATUS "Using Eigen from ISLAM_EIGEN_INCLUDE_DIR=${ISLAM_EIGEN_INCLUDE_DIR}")
        add_library(Eigen3Eigen INTERFACE)
        target_include_directories(Eigen3Eigen INTERFACE ${ISLAM_EIGEN_INCLUDE_DIR})
        add_library(Eigen3::Eigen ALIAS Eigen3Eigen)
    elseif(EXISTS "/usr/include/eigen3/Eigen/Core")
        message(STATUS "Eigen3 package config not found; using /usr/include/eigen3")
        add_library(Eigen3Eigen INTERFACE)
        target_include_directories(Eigen3Eigen INTERFACE /usr/include/eigen3)
        add_library(Eigen3::Eigen ALIAS Eigen3Eigen)
    elseif(ISLAM_ALLOW_FETCHCONTENT)
        include(FetchContent)
        message(STATUS "Eigen3 not found locally; fetching Eigen 3.4.0 via FetchContent")
        FetchContent_Declare(
            eigen_external
            URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
            DOWNLOAD_EXTRACT_TIMESTAMP TRUE
        )
        FetchContent_MakeAvailable(eigen_external)
        add_library(Eigen3Eigen INTERFACE)
        target_include_directories(Eigen3Eigen INTERFACE ${eigen_external_SOURCE_DIR})
        add_library(Eigen3::Eigen ALIAS Eigen3Eigen)
    else()
        message(FATAL_ERROR "Eigen3 not found. Set -DISLAM_EIGEN_INCLUDE_DIR=/path/to/eigen or enable -DISLAM_ALLOW_FETCHCONTENT=ON.")
    endif()
endif()

set(ISLAM_HAS_CHOLMOD FALSE)
set(ISLAM_HAS_CCOLAMD FALSE)
if(ISLAM_WITH_SUITESPARSE)
    find_package(PkgConfig QUIET)
    if(PkgConfig_FOUND)
        pkg_check_modules(CHOLMOD QUIET IMPORTED_TARGET cholmod)
        pkg_check_modules(CCOLAMD QUIET IMPORTED_TARGET ccolamd)
    endif()

    find_path(CHOLMOD_INCLUDE_DIR cholmod.h PATH_SUFFIXES suitesparse)
    find_path(CCOLAMD_INCLUDE_DIR ccolamd.h PATH_SUFFIXES suitesparse)
    find_library(CHOLMOD_LIBRARY cholmod)
    find_library(AMD_LIBRARY amd)
    find_library(CAMD_LIBRARY camd)
    find_library(COLAMD_LIBRARY colamd)
    find_library(CCOLAMD_LIBRARY ccolamd)
    find_library(SUITESPARSECONFIG_LIBRARY suitesparseconfig)
    find_package(BLAS QUIET)

    if(CHOLMOD_FOUND)
        set(ISLAM_HAS_CHOLMOD TRUE)
        set(ISLAM_SUITESPARSE_TARGET PkgConfig::CHOLMOD)
        message(STATUS "Found CHOLMOD via pkg-config")
    elseif(CHOLMOD_INCLUDE_DIR AND CHOLMOD_LIBRARY)
        add_library(ISLAM_CHOLMOD INTERFACE)
        target_include_directories(ISLAM_CHOLMOD INTERFACE ${CHOLMOD_INCLUDE_DIR})
        target_link_libraries(ISLAM_CHOLMOD INTERFACE ${CHOLMOD_LIBRARY})
        if(AMD_LIBRARY)
            target_link_libraries(ISLAM_CHOLMOD INTERFACE ${AMD_LIBRARY})
        endif()
        if(CAMD_LIBRARY)
            target_link_libraries(ISLAM_CHOLMOD INTERFACE ${CAMD_LIBRARY})
        endif()
        if(COLAMD_LIBRARY)
            target_link_libraries(ISLAM_CHOLMOD INTERFACE ${COLAMD_LIBRARY})
        endif()
        if(CCOLAMD_LIBRARY)
            target_link_libraries(ISLAM_CHOLMOD INTERFACE ${CCOLAMD_LIBRARY})
        endif()
        if(SUITESPARSECONFIG_LIBRARY)
            target_link_libraries(ISLAM_CHOLMOD INTERFACE ${SUITESPARSECONFIG_LIBRARY})
        endif()
        if(BLAS_FOUND)
            target_link_libraries(ISLAM_CHOLMOD INTERFACE BLAS::BLAS)
        endif()
        add_library(ISLAM::CHOLMOD ALIAS ISLAM_CHOLMOD)
        set(ISLAM_HAS_CHOLMOD TRUE)
        set(ISLAM_SUITESPARSE_TARGET ISLAM::CHOLMOD)
        message(STATUS "Found CHOLMOD via manual search")
    else()
        message(STATUS "CHOLMOD not found; using Eigen fallback")
    endif()

    if(CCOLAMD_FOUND)
        set(ISLAM_HAS_CCOLAMD TRUE)
        add_library(ISLAM::CCOLAMD ALIAS PkgConfig::CCOLAMD)
        message(STATUS "Found SuiteSparse CCOLAMD direct API via pkg-config")
    elseif(CCOLAMD_INCLUDE_DIR AND CCOLAMD_LIBRARY)
        add_library(ISLAM_CCOLAMD INTERFACE)
        target_include_directories(ISLAM_CCOLAMD INTERFACE ${CCOLAMD_INCLUDE_DIR})
        target_link_libraries(ISLAM_CCOLAMD INTERFACE ${CCOLAMD_LIBRARY})
        if(SUITESPARSECONFIG_LIBRARY)
            target_link_libraries(ISLAM_CCOLAMD INTERFACE ${SUITESPARSECONFIG_LIBRARY})
        endif()
        add_library(ISLAM::CCOLAMD ALIAS ISLAM_CCOLAMD)
        set(ISLAM_HAS_CCOLAMD TRUE)
        message(STATUS "Found SuiteSparse CCOLAMD direct API via manual search")
    else()
        message(STATUS "SuiteSparse CCOLAMD direct API not found")
    endif()
endif()


set(ISLAM_HAS_GTSAM FALSE)
if(ISLAM_WITH_GTSAM)
    # GTSAM installations are not completely uniform: source/build-tree installs
    # commonly provide GTSAMConfig.cmake, while some package managers expose a
    # lowercase gtsam config package. Try both before deciding that the optional
    # backend is unavailable.
    find_package(GTSAM QUIET CONFIG)
    if(NOT GTSAM_FOUND)
        find_package(gtsam QUIET CONFIG)
    endif()
    if(GTSAM_FOUND OR gtsam_FOUND)
        set(ISLAM_HAS_GTSAM TRUE)
        if(TARGET gtsam)
            set(ISLAM_GTSAM_TARGET gtsam)
        elseif(TARGET GTSAM::gtsam)
            set(ISLAM_GTSAM_TARGET GTSAM::gtsam)
        elseif(TARGET gtsam::gtsam)
            set(ISLAM_GTSAM_TARGET gtsam::gtsam)
        elseif(GTSAM_LIBRARIES)
            set(ISLAM_GTSAM_TARGET ${GTSAM_LIBRARIES})
        elseif(gtsam_LIBRARIES)
            set(ISLAM_GTSAM_TARGET ${gtsam_LIBRARIES})
        else()
            set(ISLAM_HAS_GTSAM FALSE)
        endif()
        message(STATUS "Found GTSAM for optional Option-A IGG-SPO backend")
    endif()
    if(ISLAM_WITH_GTSAM AND NOT ISLAM_HAS_GTSAM)
        if(ISLAM_REQUIRE_GTSAM)
            message(FATAL_ERROR "ISLAM_WITH_GTSAM=ON but GTSAM was not found. Set CMAKE_PREFIX_PATH or disable ISLAM_REQUIRE_GTSAM.")
        else()
            message(STATUS "GTSAM not found; skipping optional Option-A IGG-SPO backend")
        endif()
    endif()
endif()

if(ISLAM_REQUIRE_SUITESPARSE_CCOLAMD AND NOT ISLAM_HAS_CCOLAMD)
    message(FATAL_ERROR "ISLAM_REQUIRE_SUITESPARSE_CCOLAMD=ON but ccolamd.h/libccolamd were not found. Install SuiteSparse development headers/libraries or set CMAKE_PREFIX_PATH.")
endif()

set(ISLAM_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")
if(NOT EXISTS "${ISLAM_INCLUDE_DIR}/islam/factor_manager.hpp")
    set(ISLAM_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
endif()

set(ISLAM_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
if(NOT EXISTS "${ISLAM_SRC_DIR}/factor_manager.cpp")
    set(ISLAM_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
endif()

set(ISLAM_TEST_DIR "${CMAKE_CURRENT_SOURCE_DIR}/tests")
if(NOT EXISTS "${ISLAM_TEST_DIR}/test_smoke.cpp")
    set(ISLAM_TEST_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
endif()

set(ISLAM_CORE_SOURCES
    ${ISLAM_SRC_DIR}/se2.cpp
    ${ISLAM_SRC_DIR}/se3.cpp
    ${ISLAM_SRC_DIR}/io_g2o.cpp
    ${ISLAM_SRC_DIR}/io_toro.cpp
    ${ISLAM_SRC_DIR}/linearize.cpp
    ${ISLAM_SRC_DIR}/reorder_edges.cpp
    ${ISLAM_SRC_DIR}/update_graph.cpp
    ${ISLAM_SRC_DIR}/loop_closure.cpp
    ${ISLAM_SRC_DIR}/factor_manager.cpp
    ${ISLAM_SRC_DIR}/sparse_expanding_cholesky.cpp
    ${ISLAM_SRC_DIR}/affected.cpp
    ${ISLAM_SRC_DIR}/selective_solver.cpp
    ${ISLAM_SRC_DIR}/synthetic.cpp
    ${ISLAM_SRC_DIR}/benchmark.cpp
    ${ISLAM_SRC_DIR}/symbolic_engine.cpp
    ${ISLAM_SRC_DIR}/dynamic_ccolamd_engine.cpp
    ${ISLAM_SRC_DIR}/deterministic_batch_ccolamd.cpp
    ${ISLAM_SRC_DIR}/dynamic_exact_ccolamd.cpp
)
if(EXISTS "${ISLAM_SRC_DIR}/suitesparse_ccolamd_reference.cpp")
    list(APPEND ISLAM_CORE_SOURCES ${ISLAM_SRC_DIR}/suitesparse_ccolamd_reference.cpp)
endif()

add_library(islam_core ${ISLAM_CORE_SOURCES})

target_include_directories(islam_core
    PUBLIC
        ${ISLAM_INCLUDE_DIR}
)

target_link_libraries(islam_core
    PUBLIC
        Eigen3::Eigen
)

if(ISLAM_HAS_CHOLMOD)
    target_link_libraries(islam_core PUBLIC ${ISLAM_SUITESPARSE_TARGET})
    target_compile_definitions(islam_core PUBLIC ISLAM_HAS_CHOLMOD=1)
endif()

if(ISLAM_HAS_CCOLAMD)
    target_link_libraries(islam_core PUBLIC ISLAM::CCOLAMD)
    target_compile_definitions(islam_core PUBLIC ISLAM_HAS_CCOLAMD=1)
endif()

target_compile_features(islam_core PUBLIC cxx_std_17)

add_executable(islam_demo
    ${ISLAM_SRC_DIR}/main_demo.cpp
)

target_link_libraries(islam_demo PRIVATE islam_core)

add_executable(islam_run_benchmark
    ${ISLAM_SRC_DIR}/run_benchmark.cpp
)

target_link_libraries(islam_run_benchmark PRIVATE islam_core)


if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tools/suitesparse_compat_audit.cpp")
    add_executable(islam_suitesparse_compat_audit
        ${CMAKE_CURRENT_SOURCE_DIR}/tools/suitesparse_compat_audit.cpp
    )
    target_link_libraries(islam_suitesparse_compat_audit PRIVATE islam_core)
endif()


if(ISLAM_HAS_GTSAM)
    add_library(islam_gtsam_igg_spo
        ${ISLAM_SRC_DIR}/gtsam_igg_spo.cpp
    )
    target_link_libraries(islam_gtsam_igg_spo
        PUBLIC
            islam_core
            ${ISLAM_GTSAM_TARGET}
    )
    target_compile_definitions(islam_gtsam_igg_spo PUBLIC ISLAM_HAS_GTSAM=1)
    if(GTSAM_INCLUDE_DIRS)
        target_include_directories(islam_gtsam_igg_spo PUBLIC ${GTSAM_INCLUDE_DIRS})
    elseif(GTSAM_INCLUDE_DIR)
        target_include_directories(islam_gtsam_igg_spo PUBLIC ${GTSAM_INCLUDE_DIR})
    endif()
    target_compile_features(islam_gtsam_igg_spo PUBLIC cxx_std_17)

    add_executable(islam_gtsam_igg_spo_example
        ${CMAKE_CURRENT_SOURCE_DIR}/tools/gtsam_igg_spo_example.cpp
    )
    target_link_libraries(islam_gtsam_igg_spo_example PRIVATE islam_gtsam_igg_spo)
endif()


if(ISLAM_BUILD_TESTS)
    include(CTest)
    enable_testing()

    add_executable(islam_test_smoke
        ${ISLAM_TEST_DIR}/test_smoke.cpp
    )
    target_link_libraries(islam_test_smoke PRIVATE islam_core)
    add_test(NAME smoke COMMAND islam_test_smoke)

    if(EXISTS "${ISLAM_TEST_DIR}/test_sparse_expanding_cholesky_randomized.cpp")
        add_executable(islam_test_sparse_expanding_cholesky_randomized
            ${ISLAM_TEST_DIR}/test_sparse_expanding_cholesky_randomized.cpp
            ${ISLAM_SRC_DIR}/sparse_expanding_cholesky.cpp
        )
    elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/test_sparse_expanding_cholesky_randomized.cpp")
        add_executable(islam_test_sparse_expanding_cholesky_randomized
            ${CMAKE_CURRENT_SOURCE_DIR}/test_sparse_expanding_cholesky_randomized.cpp
            ${ISLAM_SRC_DIR}/sparse_expanding_cholesky.cpp
        )
    endif()
    if(TARGET islam_test_sparse_expanding_cholesky_randomized)
        target_link_libraries(islam_test_sparse_expanding_cholesky_randomized PRIVATE Eigen3::Eigen)
        target_include_directories(islam_test_sparse_expanding_cholesky_randomized PRIVATE ${ISLAM_INCLUDE_DIR})
        target_compile_features(islam_test_sparse_expanding_cholesky_randomized PRIVATE cxx_std_17)
        add_test(NAME sparse_expanding_cholesky_randomized
                 COMMAND islam_test_sparse_expanding_cholesky_randomized)
    endif()

    if(EXISTS "${ISLAM_TEST_DIR}/test_sparse_expanding_cholesky_trace_replay.cpp")
        add_executable(islam_test_sparse_expanding_cholesky_trace_replay
            ${ISLAM_TEST_DIR}/test_sparse_expanding_cholesky_trace_replay.cpp
            ${ISLAM_SRC_DIR}/sparse_expanding_cholesky.cpp
        )
    elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/test_sparse_expanding_cholesky_trace_replay.cpp")
        add_executable(islam_test_sparse_expanding_cholesky_trace_replay
            ${CMAKE_CURRENT_SOURCE_DIR}/test_sparse_expanding_cholesky_trace_replay.cpp
            ${ISLAM_SRC_DIR}/sparse_expanding_cholesky.cpp
        )
    endif()
    if(TARGET islam_test_sparse_expanding_cholesky_trace_replay)
        target_link_libraries(islam_test_sparse_expanding_cholesky_trace_replay PRIVATE Eigen3::Eigen)
        target_include_directories(islam_test_sparse_expanding_cholesky_trace_replay PRIVATE ${ISLAM_INCLUDE_DIR})
        target_compile_features(islam_test_sparse_expanding_cholesky_trace_replay PRIVATE cxx_std_17)
        add_test(NAME sparse_expanding_cholesky_trace_replay
                 COMMAND islam_test_sparse_expanding_cholesky_trace_replay)
    endif()
endif()

message(STATUS "ISLAM complete bundle configured. CHOLMOD enabled: ${ISLAM_HAS_CHOLMOD}; CCOLAMD direct API enabled: ${ISLAM_HAS_CCOLAMD}; GTSAM Option-A backend enabled: ${ISLAM_HAS_GTSAM}")
