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
include contrib/romano/CMakeFiles/DirkInteraction.dir/depend.make

# Include the progress variables for this target.
include contrib/romano/CMakeFiles/DirkInteraction.dir/progress.make

# Include the compile flags for this target's objects.
include contrib/romano/CMakeFiles/DirkInteraction.dir/flags.make

contrib/romano/CMakeFiles/DirkInteraction.dir/src/Interactions/DirkInteraction.cpp.o: contrib/romano/CMakeFiles/DirkInteraction.dir/flags.make
contrib/romano/CMakeFiles/DirkInteraction.dir/src/Interactions/DirkInteraction.cpp.o: ../contrib/romano/src/Interactions/DirkInteraction.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/josh/git/oxDNA_torsion/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object contrib/romano/CMakeFiles/DirkInteraction.dir/src/Interactions/DirkInteraction.cpp.o"
	cd /home/josh/git/oxDNA_torsion/build/contrib/romano && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/DirkInteraction.dir/src/Interactions/DirkInteraction.cpp.o -c /home/josh/git/oxDNA_torsion/contrib/romano/src/Interactions/DirkInteraction.cpp

contrib/romano/CMakeFiles/DirkInteraction.dir/src/Interactions/DirkInteraction.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DirkInteraction.dir/src/Interactions/DirkInteraction.cpp.i"
	cd /home/josh/git/oxDNA_torsion/build/contrib/romano && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/josh/git/oxDNA_torsion/contrib/romano/src/Interactions/DirkInteraction.cpp > CMakeFiles/DirkInteraction.dir/src/Interactions/DirkInteraction.cpp.i

contrib/romano/CMakeFiles/DirkInteraction.dir/src/Interactions/DirkInteraction.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DirkInteraction.dir/src/Interactions/DirkInteraction.cpp.s"
	cd /home/josh/git/oxDNA_torsion/build/contrib/romano && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/josh/git/oxDNA_torsion/contrib/romano/src/Interactions/DirkInteraction.cpp -o CMakeFiles/DirkInteraction.dir/src/Interactions/DirkInteraction.cpp.s

# Object files for target DirkInteraction
DirkInteraction_OBJECTS = \
"CMakeFiles/DirkInteraction.dir/src/Interactions/DirkInteraction.cpp.o"

# External object files for target DirkInteraction
DirkInteraction_EXTERNAL_OBJECTS =

../contrib/romano/DirkInteraction.so: contrib/romano/CMakeFiles/DirkInteraction.dir/src/Interactions/DirkInteraction.cpp.o
../contrib/romano/DirkInteraction.so: contrib/romano/CMakeFiles/DirkInteraction.dir/build.make
../contrib/romano/DirkInteraction.so: contrib/romano/CMakeFiles/DirkInteraction.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/josh/git/oxDNA_torsion/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared library ../../../contrib/romano/DirkInteraction.so"
	cd /home/josh/git/oxDNA_torsion/build/contrib/romano && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/DirkInteraction.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
contrib/romano/CMakeFiles/DirkInteraction.dir/build: ../contrib/romano/DirkInteraction.so

.PHONY : contrib/romano/CMakeFiles/DirkInteraction.dir/build

contrib/romano/CMakeFiles/DirkInteraction.dir/clean:
	cd /home/josh/git/oxDNA_torsion/build/contrib/romano && $(CMAKE_COMMAND) -P CMakeFiles/DirkInteraction.dir/cmake_clean.cmake
.PHONY : contrib/romano/CMakeFiles/DirkInteraction.dir/clean

contrib/romano/CMakeFiles/DirkInteraction.dir/depend:
	cd /home/josh/git/oxDNA_torsion/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/josh/git/oxDNA_torsion /home/josh/git/oxDNA_torsion/contrib/romano /home/josh/git/oxDNA_torsion/build /home/josh/git/oxDNA_torsion/build/contrib/romano /home/josh/git/oxDNA_torsion/build/contrib/romano/CMakeFiles/DirkInteraction.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : contrib/romano/CMakeFiles/DirkInteraction.dir/depend

