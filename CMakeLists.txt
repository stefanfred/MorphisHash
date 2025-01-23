cmake_minimum_required (VERSION 3.25)
project(MorphisHash
        LANGUAGES CXX)

if(TARGET MorphisHash)
    return()
endif()

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif ()

if((CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo") AND PROJECT_IS_TOP_LEVEL)
    add_compile_options(-march=native)
endif()

# ---------------------------- Dependencies ----------------------------

if(NOT TARGET Sux)
    add_library(Sux INTERFACE)
    target_include_directories(Sux SYSTEM INTERFACE extlib/sux)
endif()

if(NOT TARGET tlx)
    set(TLX_INSTALL_INCLUDE_DIR tlx CACHE PATH "Workaround for TLX breaking the first cmake call")
    add_subdirectory(extlib/tlx SYSTEM EXCLUDE_FROM_ALL)
endif()

if(NOT TARGET sdsl)
    add_subdirectory(extlib/sdsl-lite/external SYSTEM EXCLUDE_FROM_ALL)
    add_subdirectory(extlib/sdsl-lite/include SYSTEM EXCLUDE_FROM_ALL)
    add_subdirectory(extlib/sdsl-lite/lib SYSTEM EXCLUDE_FROM_ALL)
    target_compile_options(sdsl PRIVATE -w)
    target_include_directories(sdsl SYSTEM INTERFACE extlib/sdsl-lite/include)
endif()


if(NOT TARGET ips2ra)
    set(IPS2RA_DISABLE_PARALLEL ON CACHE PATH "ips2ra's FindTBB greps a file that does not exist in recent TBB versions")
    add_subdirectory(extlib/ips2ra SYSTEM)
    find_package(TBB)
    target_compile_options(ips2ra INTERFACE -D_REENTRANT)
    target_link_libraries(ips2ra INTERFACE pthread atomic TBB::tbb)
endif()

if(NOT TARGET Ips2raMorphisHashSorter)
    add_library(Ips2raMorphisHashSorter SHARED src/Sorter.cpp)
    target_compile_features(Ips2raMorphisHashSorter PRIVATE cxx_std_20)
    target_include_directories(Ips2raMorphisHashSorter PRIVATE include)
    target_link_libraries(Ips2raMorphisHashSorter PUBLIC ips2ra tlx Sux)
    target_compile_options(Ips2raMorphisHashSorter PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-march=native>)
endif()

if(NOT TARGET ByteHamsterUtil)
    add_subdirectory(extlib/util EXCLUDE_FROM_ALL)
endif()

# ---------------------------- Library Setup ----------------------------

add_library(MorphisHash2Precompiled SHARED src/MorphisHash-precompiled.cpp)
target_compile_features(MorphisHash2Precompiled PRIVATE cxx_std_20)
target_include_directories(MorphisHash2Precompiled PRIVATE include)
target_link_libraries(MorphisHash2Precompiled PUBLIC ips2ra tlx Sux ByteHamsterUtil)
target_compile_options(MorphisHash2Precompiled PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-march=native>)

add_library(MorphisHash INTERFACE)
target_include_directories(MorphisHash INTERFACE include)
target_compile_features(MorphisHash INTERFACE cxx_std_20)
target_link_libraries(MorphisHash INTERFACE Sux Ips2raMorphisHashSorter ByteHamsterUtil MorphisHash2Precompiled sdsl)



# ---------------------------- Benchmarks ----------------------------
if(PROJECT_IS_TOP_LEVEL)
    add_library(BenchmarkUtils INTERFACE)
    target_include_directories(BenchmarkUtils INTERFACE benchmark)
    target_include_directories(BenchmarkUtils INTERFACE test)
    target_link_libraries(BenchmarkUtils INTERFACE tlx ByteHamster::Util)

    # Warnings if this is the main project
    target_compile_options(MorphisHash INTERFACE $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Wpedantic -Werror -Wno-error=stringop-overflow -frecord-gcc-switches>)

    add_executable(Benchmark benchmark/benchmark_construction.cpp)
    target_link_libraries(Benchmark PUBLIC BenchmarkUtils MorphisHash)

    add_executable(Playground benchmark/playground.cpp)
    target_link_libraries(Playground PUBLIC BenchmarkUtils MorphisHash)

    add_executable(NumHashEvals benchmark/numHashEvals.cpp)
    target_link_libraries(NumHashEvals PUBLIC BenchmarkUtils MorphisHash)

    add_executable(GolombMemoTuner benchmark/golombMemoTuner.cpp)
    target_link_libraries(GolombMemoTuner PUBLIC MorphisHash BenchmarkUtils)
endif()