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
include contrib/rovigatti/CMakeFiles/PatchyToMgl.dir/depend.make

# Include the progress variables for this target.
include contrib/rovigatti/CMakeFiles/PatchyToMgl.dir/progress.make

# Include the compile flags for this target's objects.
include contrib/rovigatti/CMakeFiles/PatchyToMgl.dir/flags.make

contrib/rovigatti/CMakeFiles/PatchyToMgl.dir/src/Observables/PatchyToMgl.cpp.o: contrib/rovigatti/CMakeFiles/PatchyToMgl.dir/flags.make
contrib/rovigatti/CMakeFiles/PatchyToMgl.dir/src/Observables/PatchyToMgl.cpp.o: ../contrib/rovigatti/src/Observables/PatchyToMgl.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/josh/git/oxDNA_torsion/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object contrib/rovigatti/CMakeFiles/PatchyToMgl.dir/src/Observables/PatchyToMgl.cpp.o"
	cd /home/josh/git/oxDNA_torsion/build/contrib/rovigatti && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/PatchyToMgl.dir/src/Observables/PatchyToMgl.cpp.o -c /home/josh/git/oxDNA_torsion/contrib/rovigatti/src/Observables/PatchyToMgl.cpp

contrib/rovigatti/CMakeFiles/PatchyToMgl.dir/src/Observables/PatchyToMgl.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/PatchyToMgl.dir/src/Observables/PatchyToMgl.cpp.i"
	cd /home/josh/git/oxDNA_torsion/build/contrib/rovigatti && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/josh/git/oxDNA_torsion/contrib/rovigatti/src/Observables/PatchyToMgl.cpp > CMakeFiles/PatchyToMgl.dir/src/Observables/PatchyToMgl.cpp.i

contrib/rovigatti/CMakeFiles/PatchyToMgl.dir/src/Observables/PatchyToMgl.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/PatchyToMgl.dir/src/Observables/PatchyToMgl.cpp.s"
	cd /home/josh/git/oxDNA_torsion/build/contrib/rovigatti && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/josh/git/oxDNA_torsion/contrib/rovigatti/src/Observables/PatchyToMgl.cpp -o CMakeFiles/PatchyToMgl.dir/src/Observables/PatchyToMgl.cpp.s

# Object files for target PatchyToMgl
PatchyToMgl_OBJECTS = \
"CMakeFiles/PatchyToMgl.dir/src/Observables/PatchyToMgl.cpp.o"

# External object files for target PatchyToMgl
PatchyToMgl_EXTERNAL_OBJECTS =

../contrib/rovigatti/PatchyToMgl.so: contrib/rovigatti/CMakeFiles/PatchyToMgl.dir/src/Observables/PatchyToMgl.cpp.o
../contrib/rovigatti/PatchyToMgl.so: contrib/rovigatti/CMakeFiles/PatchyToMgl.dir/build.make
../contrib/rovigatti/PatchyToMgl.so: contrib/rovigatti/CMakeFiles/PatchyToMgl.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/josh/git/oxDNA_torsion/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared library ../../../contrib/rovigatti/PatchyToMgl.so"
	cd /home/josh/git/oxDNA_torsion/build/contrib/rovigatti && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/PatchyToMgl.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
contrib/rovigatti/CMakeFiles/PatchyToMgl.dir/build: ../contrib/rovigatti/PatchyToMgl.so

.PHONY : contrib/rovigatti/CMakeFiles/PatchyToMgl.dir/build

contrib/rovigatti/CMakeFiles/PatchyToMgl.dir/clean:
	cd /home/josh/git/oxDNA_torsion/build/contrib/rovigatti && $(CMAKE_COMMAND) -P CMakeFiles/PatchyToMgl.dir/cmake_clean.cmake
.PHONY : contrib/rovigatti/CMakeFiles/PatchyToMgl.dir/clean

contrib/rovigatti/CMakeFiles/PatchyToMgl.dir/depend:
	cd /home/josh/git/oxDNA_torsion/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/josh/git/oxDNA_torsion /home/josh/git/oxDNA_torsion/contrib/rovigatti /home/josh/git/oxDNA_torsion/build /home/josh/git/oxDNA_torsion/build/contrib/rovigatti /home/josh/git/oxDNA_torsion/build/contrib/rovigatti/CMakeFiles/PatchyToMgl.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : contrib/rovigatti/CMakeFiles/PatchyToMgl.dir/depend

