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
include contrib/rovigatti/CMakeFiles/GraftedInteraction.dir/depend.make

# Include the progress variables for this target.
include contrib/rovigatti/CMakeFiles/GraftedInteraction.dir/progress.make

# Include the compile flags for this target's objects.
include contrib/rovigatti/CMakeFiles/GraftedInteraction.dir/flags.make

contrib/rovigatti/CMakeFiles/GraftedInteraction.dir/src/Interactions/GraftedInteraction.cpp.o: contrib/rovigatti/CMakeFiles/GraftedInteraction.dir/flags.make
contrib/rovigatti/CMakeFiles/GraftedInteraction.dir/src/Interactions/GraftedInteraction.cpp.o: ../contrib/rovigatti/src/Interactions/GraftedInteraction.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/josh/git/oxDNA_torsion/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object contrib/rovigatti/CMakeFiles/GraftedInteraction.dir/src/Interactions/GraftedInteraction.cpp.o"
	cd /home/josh/git/oxDNA_torsion/build/contrib/rovigatti && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/GraftedInteraction.dir/src/Interactions/GraftedInteraction.cpp.o -c /home/josh/git/oxDNA_torsion/contrib/rovigatti/src/Interactions/GraftedInteraction.cpp

contrib/rovigatti/CMakeFiles/GraftedInteraction.dir/src/Interactions/GraftedInteraction.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/GraftedInteraction.dir/src/Interactions/GraftedInteraction.cpp.i"
	cd /home/josh/git/oxDNA_torsion/build/contrib/rovigatti && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/josh/git/oxDNA_torsion/contrib/rovigatti/src/Interactions/GraftedInteraction.cpp > CMakeFiles/GraftedInteraction.dir/src/Interactions/GraftedInteraction.cpp.i

contrib/rovigatti/CMakeFiles/GraftedInteraction.dir/src/Interactions/GraftedInteraction.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/GraftedInteraction.dir/src/Interactions/GraftedInteraction.cpp.s"
	cd /home/josh/git/oxDNA_torsion/build/contrib/rovigatti && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/josh/git/oxDNA_torsion/contrib/rovigatti/src/Interactions/GraftedInteraction.cpp -o CMakeFiles/GraftedInteraction.dir/src/Interactions/GraftedInteraction.cpp.s

# Object files for target GraftedInteraction
GraftedInteraction_OBJECTS = \
"CMakeFiles/GraftedInteraction.dir/src/Interactions/GraftedInteraction.cpp.o"

# External object files for target GraftedInteraction
GraftedInteraction_EXTERNAL_OBJECTS =

../contrib/rovigatti/GraftedInteraction.so: contrib/rovigatti/CMakeFiles/GraftedInteraction.dir/src/Interactions/GraftedInteraction.cpp.o
../contrib/rovigatti/GraftedInteraction.so: contrib/rovigatti/CMakeFiles/GraftedInteraction.dir/build.make
../contrib/rovigatti/GraftedInteraction.so: contrib/rovigatti/CMakeFiles/GraftedInteraction.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/josh/git/oxDNA_torsion/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared library ../../../contrib/rovigatti/GraftedInteraction.so"
	cd /home/josh/git/oxDNA_torsion/build/contrib/rovigatti && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/GraftedInteraction.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
contrib/rovigatti/CMakeFiles/GraftedInteraction.dir/build: ../contrib/rovigatti/GraftedInteraction.so

.PHONY : contrib/rovigatti/CMakeFiles/GraftedInteraction.dir/build

contrib/rovigatti/CMakeFiles/GraftedInteraction.dir/clean:
	cd /home/josh/git/oxDNA_torsion/build/contrib/rovigatti && $(CMAKE_COMMAND) -P CMakeFiles/GraftedInteraction.dir/cmake_clean.cmake
.PHONY : contrib/rovigatti/CMakeFiles/GraftedInteraction.dir/clean

contrib/rovigatti/CMakeFiles/GraftedInteraction.dir/depend:
	cd /home/josh/git/oxDNA_torsion/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/josh/git/oxDNA_torsion /home/josh/git/oxDNA_torsion/contrib/rovigatti /home/josh/git/oxDNA_torsion/build /home/josh/git/oxDNA_torsion/build/contrib/rovigatti /home/josh/git/oxDNA_torsion/build/contrib/rovigatti/CMakeFiles/GraftedInteraction.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : contrib/rovigatti/CMakeFiles/GraftedInteraction.dir/depend

