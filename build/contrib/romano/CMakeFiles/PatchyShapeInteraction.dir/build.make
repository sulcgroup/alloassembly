# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/josh/git/oxDNA_torsion

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/josh/git/oxDNA_torsion/build

# Include any dependencies generated for this target.
include contrib/romano/CMakeFiles/PatchyShapeInteraction.dir/depend.make

# Include the progress variables for this target.
include contrib/romano/CMakeFiles/PatchyShapeInteraction.dir/progress.make

# Include the compile flags for this target's objects.
include contrib/romano/CMakeFiles/PatchyShapeInteraction.dir/flags.make

contrib/romano/CMakeFiles/PatchyShapeInteraction.dir/src/Interactions/PatchyShapeInteraction.cpp.o: contrib/romano/CMakeFiles/PatchyShapeInteraction.dir/flags.make
contrib/romano/CMakeFiles/PatchyShapeInteraction.dir/src/Interactions/PatchyShapeInteraction.cpp.o: ../contrib/romano/src/Interactions/PatchyShapeInteraction.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/josh/git/oxDNA_torsion/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object contrib/romano/CMakeFiles/PatchyShapeInteraction.dir/src/Interactions/PatchyShapeInteraction.cpp.o"
	cd /home/josh/git/oxDNA_torsion/build/contrib/romano && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/PatchyShapeInteraction.dir/src/Interactions/PatchyShapeInteraction.cpp.o -c /home/josh/git/oxDNA_torsion/contrib/romano/src/Interactions/PatchyShapeInteraction.cpp

contrib/romano/CMakeFiles/PatchyShapeInteraction.dir/src/Interactions/PatchyShapeInteraction.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/PatchyShapeInteraction.dir/src/Interactions/PatchyShapeInteraction.cpp.i"
	cd /home/josh/git/oxDNA_torsion/build/contrib/romano && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/josh/git/oxDNA_torsion/contrib/romano/src/Interactions/PatchyShapeInteraction.cpp > CMakeFiles/PatchyShapeInteraction.dir/src/Interactions/PatchyShapeInteraction.cpp.i

contrib/romano/CMakeFiles/PatchyShapeInteraction.dir/src/Interactions/PatchyShapeInteraction.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/PatchyShapeInteraction.dir/src/Interactions/PatchyShapeInteraction.cpp.s"
	cd /home/josh/git/oxDNA_torsion/build/contrib/romano && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/josh/git/oxDNA_torsion/contrib/romano/src/Interactions/PatchyShapeInteraction.cpp -o CMakeFiles/PatchyShapeInteraction.dir/src/Interactions/PatchyShapeInteraction.cpp.s

contrib/romano/CMakeFiles/PatchyShapeInteraction.dir/src/Particles/PatchyShapeParticle.cpp.o: contrib/romano/CMakeFiles/PatchyShapeInteraction.dir/flags.make
contrib/romano/CMakeFiles/PatchyShapeInteraction.dir/src/Particles/PatchyShapeParticle.cpp.o: ../contrib/romano/src/Particles/PatchyShapeParticle.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/josh/git/oxDNA_torsion/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object contrib/romano/CMakeFiles/PatchyShapeInteraction.dir/src/Particles/PatchyShapeParticle.cpp.o"
	cd /home/josh/git/oxDNA_torsion/build/contrib/romano && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/PatchyShapeInteraction.dir/src/Particles/PatchyShapeParticle.cpp.o -c /home/josh/git/oxDNA_torsion/contrib/romano/src/Particles/PatchyShapeParticle.cpp

contrib/romano/CMakeFiles/PatchyShapeInteraction.dir/src/Particles/PatchyShapeParticle.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/PatchyShapeInteraction.dir/src/Particles/PatchyShapeParticle.cpp.i"
	cd /home/josh/git/oxDNA_torsion/build/contrib/romano && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/josh/git/oxDNA_torsion/contrib/romano/src/Particles/PatchyShapeParticle.cpp > CMakeFiles/PatchyShapeInteraction.dir/src/Particles/PatchyShapeParticle.cpp.i

contrib/romano/CMakeFiles/PatchyShapeInteraction.dir/src/Particles/PatchyShapeParticle.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/PatchyShapeInteraction.dir/src/Particles/PatchyShapeParticle.cpp.s"
	cd /home/josh/git/oxDNA_torsion/build/contrib/romano && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/josh/git/oxDNA_torsion/contrib/romano/src/Particles/PatchyShapeParticle.cpp -o CMakeFiles/PatchyShapeInteraction.dir/src/Particles/PatchyShapeParticle.cpp.s

# Object files for target PatchyShapeInteraction
PatchyShapeInteraction_OBJECTS = \
"CMakeFiles/PatchyShapeInteraction.dir/src/Interactions/PatchyShapeInteraction.cpp.o" \
"CMakeFiles/PatchyShapeInteraction.dir/src/Particles/PatchyShapeParticle.cpp.o"

# External object files for target PatchyShapeInteraction
PatchyShapeInteraction_EXTERNAL_OBJECTS =

../contrib/romano/PatchyShapeInteraction.so: contrib/romano/CMakeFiles/PatchyShapeInteraction.dir/src/Interactions/PatchyShapeInteraction.cpp.o
../contrib/romano/PatchyShapeInteraction.so: contrib/romano/CMakeFiles/PatchyShapeInteraction.dir/src/Particles/PatchyShapeParticle.cpp.o
../contrib/romano/PatchyShapeInteraction.so: contrib/romano/CMakeFiles/PatchyShapeInteraction.dir/build.make
../contrib/romano/PatchyShapeInteraction.so: contrib/romano/CMakeFiles/PatchyShapeInteraction.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/josh/git/oxDNA_torsion/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX shared library ../../../contrib/romano/PatchyShapeInteraction.so"
	cd /home/josh/git/oxDNA_torsion/build/contrib/romano && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/PatchyShapeInteraction.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
contrib/romano/CMakeFiles/PatchyShapeInteraction.dir/build: ../contrib/romano/PatchyShapeInteraction.so

.PHONY : contrib/romano/CMakeFiles/PatchyShapeInteraction.dir/build

contrib/romano/CMakeFiles/PatchyShapeInteraction.dir/clean:
	cd /home/josh/git/oxDNA_torsion/build/contrib/romano && $(CMAKE_COMMAND) -P CMakeFiles/PatchyShapeInteraction.dir/cmake_clean.cmake
.PHONY : contrib/romano/CMakeFiles/PatchyShapeInteraction.dir/clean

contrib/romano/CMakeFiles/PatchyShapeInteraction.dir/depend:
	cd /home/josh/git/oxDNA_torsion/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/josh/git/oxDNA_torsion /home/josh/git/oxDNA_torsion/contrib/romano /home/josh/git/oxDNA_torsion/build /home/josh/git/oxDNA_torsion/build/contrib/romano /home/josh/git/oxDNA_torsion/build/contrib/romano/CMakeFiles/PatchyShapeInteraction.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : contrib/romano/CMakeFiles/PatchyShapeInteraction.dir/depend

