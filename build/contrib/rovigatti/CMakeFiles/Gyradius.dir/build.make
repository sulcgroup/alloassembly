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
include contrib/rovigatti/CMakeFiles/Gyradius.dir/depend.make

# Include the progress variables for this target.
include contrib/rovigatti/CMakeFiles/Gyradius.dir/progress.make

# Include the compile flags for this target's objects.
include contrib/rovigatti/CMakeFiles/Gyradius.dir/flags.make

contrib/rovigatti/CMakeFiles/Gyradius.dir/src/Observables/Gyradius.cpp.o: contrib/rovigatti/CMakeFiles/Gyradius.dir/flags.make
contrib/rovigatti/CMakeFiles/Gyradius.dir/src/Observables/Gyradius.cpp.o: ../contrib/rovigatti/src/Observables/Gyradius.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/josh/git/oxDNA_torsion/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object contrib/rovigatti/CMakeFiles/Gyradius.dir/src/Observables/Gyradius.cpp.o"
	cd /home/josh/git/oxDNA_torsion/build/contrib/rovigatti && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Gyradius.dir/src/Observables/Gyradius.cpp.o -c /home/josh/git/oxDNA_torsion/contrib/rovigatti/src/Observables/Gyradius.cpp

contrib/rovigatti/CMakeFiles/Gyradius.dir/src/Observables/Gyradius.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Gyradius.dir/src/Observables/Gyradius.cpp.i"
	cd /home/josh/git/oxDNA_torsion/build/contrib/rovigatti && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/josh/git/oxDNA_torsion/contrib/rovigatti/src/Observables/Gyradius.cpp > CMakeFiles/Gyradius.dir/src/Observables/Gyradius.cpp.i

contrib/rovigatti/CMakeFiles/Gyradius.dir/src/Observables/Gyradius.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Gyradius.dir/src/Observables/Gyradius.cpp.s"
	cd /home/josh/git/oxDNA_torsion/build/contrib/rovigatti && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/josh/git/oxDNA_torsion/contrib/rovigatti/src/Observables/Gyradius.cpp -o CMakeFiles/Gyradius.dir/src/Observables/Gyradius.cpp.s

# Object files for target Gyradius
Gyradius_OBJECTS = \
"CMakeFiles/Gyradius.dir/src/Observables/Gyradius.cpp.o"

# External object files for target Gyradius
Gyradius_EXTERNAL_OBJECTS =

../contrib/rovigatti/Gyradius.so: contrib/rovigatti/CMakeFiles/Gyradius.dir/src/Observables/Gyradius.cpp.o
../contrib/rovigatti/Gyradius.so: contrib/rovigatti/CMakeFiles/Gyradius.dir/build.make
../contrib/rovigatti/Gyradius.so: contrib/rovigatti/CMakeFiles/Gyradius.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/josh/git/oxDNA_torsion/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared library ../../../contrib/rovigatti/Gyradius.so"
	cd /home/josh/git/oxDNA_torsion/build/contrib/rovigatti && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Gyradius.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
contrib/rovigatti/CMakeFiles/Gyradius.dir/build: ../contrib/rovigatti/Gyradius.so

.PHONY : contrib/rovigatti/CMakeFiles/Gyradius.dir/build

contrib/rovigatti/CMakeFiles/Gyradius.dir/clean:
	cd /home/josh/git/oxDNA_torsion/build/contrib/rovigatti && $(CMAKE_COMMAND) -P CMakeFiles/Gyradius.dir/cmake_clean.cmake
.PHONY : contrib/rovigatti/CMakeFiles/Gyradius.dir/clean

contrib/rovigatti/CMakeFiles/Gyradius.dir/depend:
	cd /home/josh/git/oxDNA_torsion/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/josh/git/oxDNA_torsion /home/josh/git/oxDNA_torsion/contrib/rovigatti /home/josh/git/oxDNA_torsion/build /home/josh/git/oxDNA_torsion/build/contrib/rovigatti /home/josh/git/oxDNA_torsion/build/contrib/rovigatti/CMakeFiles/Gyradius.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : contrib/rovigatti/CMakeFiles/Gyradius.dir/depend

