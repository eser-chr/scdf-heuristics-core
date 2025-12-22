add_library(core STATIC)
target_include_directories(core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_sources(core
PRIVATE
    src/beam_search.cpp
    src/clustering.cpp
    src/construction.cpp
    src/grasp.cpp
    src/instance.cpp
    src/local_search.cpp
    src/neighborhoods.cpp
    src/random.cpp
    src/sa.cpp
    src/solution.cpp
    src/utils.cpp
    src/vnd.cpp
)