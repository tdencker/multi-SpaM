add_subdirectory(rasbhari)
add_subdirectory(raxml)

# TODO: linux and sse3 is required

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")

add_executable(multi-SpaM
    component.cpp
    options.cpp
    multi-SpaM.cpp
    pattern.cpp
    quartetblock.cpp
    randommatchfinder.cpp
    sequence.cpp
    stats.cpp
    word.cpp
)

find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    target_compile_options(multi-SpaM PRIVATE ${OpenMP_CXX_FLAGS})
    target_link_libraries(multi-SpaM PRIVATE ${OpenMP_CXX_FLAGS})
endif()

target_link_libraries(multi-SpaM PRIVATE rasbhari raxml)