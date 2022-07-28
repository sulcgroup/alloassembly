# Written 11 June 2022 by Josh Evans
# Largely based on Lorenzo Rovigetti's CMakeLists.txt file


SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(../rovigatti)
INCLUDE_DIRECTORIES(../../src)

# Probably going to be more of these
ADD_CUSTOM_TARGET(evans DEPENDS AllostericPatchySwapInteraction)

# without this bit clang-based compilation fails with some "undefined symbols" errors during the linking stage
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    link_libraries(${ARGV0} "-undefined dynamic_lookup")
endif()

# this removes the prefix from each library. Lorenzo didn't find any other way of doing it for all
# libraries at once if the library type is MODULE
function(add_library_no_prefix target source)
    # ${ARGN} holds the list of arguments past the last expected argument, which is "source"
    add_library(${target} MODULE EXCLUDE_FROM_ALL ${source} ${ARGN})
    set_target_properties(${target} PROPERTIES PREFIX "")
endfunction()

add_library_no_prefix(AllostericPatchySwapInteraction src/Interactions/AllostericPatchySwapInteraction.cpp src/Particles/AllostericPatchyParticle.cpp src/Particles/AllostericPatch.cpp)
#add_library_no_prefix(AllostericPatchyParticle src/Particles/AllostericPatchyParticle.cpp)
#add_library_no_prefix(AllostericPatch src/Particles/AllostericPatch.cpp)

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