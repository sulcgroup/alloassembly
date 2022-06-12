# Written 11 June 2022 by Josh Evans
# Largely based on Lorenzo Rovigetti's CMakeLists.txt file

# require romano as a dependency
# I hope this works!
INCLUDE_DIRECTORIES(../romano)

SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# Probably going to be more of these
ADD_CUSTOM_TARGET(evans DEPENDS AllostericPatchyInteraction AllostericPatchyParticle)

# without this bit clang-based compilation fails with some "undefined symbols" errors during the linking stage
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    link_libraries(${ARGV0} "-undefined dynamic_lookup")
endif()

# this removes the prefix from each library. I didn't find any other way of doing it for all
# libraries at once if the library type is MODULE
function(add_library_no_prefix)
    add_library(${ARGN})
    set_target_properties(${ARGV0} PROPERTIES PREFIX "")
endfunction()

add_library_no_prefix(AllostericPatchyInteraction MODULE EXCLUDE_FROM_ALL src/Interactions/AllostericPatchyInteraction.cpp)
add_library_no_prefix(AllostericPatchyParticle MODULE EXCLUDE_FROM_ALL src/Particles/AllostericPatchyParticle.cpp)

if (CUDA)
    SET(CMAKE_MODULE_LINKER_FLAGS  "${CMAKE_MODULE_LINKER_FLAGS} -rdynamic")

    # same thing but for CUDA libs
    function(cuda_add_library_no_prefix)
        cuda_add_library(${ARGN})
        set_target_properties(${ARGV0} PROPERTIES PREFIX "")
    endfunction()

    cuda_add_library_no_prefix(CUDAAllostericPatchyInteraction MODULE EXCLUDE_FROM_ALL src/Interactions/CUDAAllostericPatchyInteraction.cu src/Interactions/CUDAAllostericPatchyInteraction.cpp)

    ADD_DEPENDENCIES(evans CUDAAllostericPatchyInteraction)
endif(CUDA)