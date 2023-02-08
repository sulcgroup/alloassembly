# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.23

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /snap/clion/203/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /snap/clion/203/bin/cmake/linux/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/josh/git/alloassembly

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/josh/git/alloassembly/build

# Include any dependencies generated for this target.
include contrib/romano/CMakeFiles/CutVolume.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include contrib/romano/CMakeFiles/CutVolume.dir/compiler_depend.make

# Include the progress variables for this target.
include contrib/romano/CMakeFiles/CutVolume.dir/progress.make

# Include the compile flags for this target's objects.
include contrib/romano/CMakeFiles/CutVolume.dir/flags.make

contrib/romano/CMakeFiles/CutVolume.dir/src/Backends/MCMoves/CutVolume.cpp.o: contrib/romano/CMakeFiles/CutVolume.dir/flags.make
contrib/romano/CMakeFiles/CutVolume.dir/src/Backends/MCMoves/CutVolume.cpp.o: ../contrib/romano/src/Backends/MCMoves/CutVolume.cpp
contrib/romano/CMakeFiles/CutVolume.dir/src/Backends/MCMoves/CutVolume.cpp.o: contrib/romano/CMakeFiles/CutVolume.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/josh/git/alloassembly/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object contrib/romano/CMakeFiles/CutVolume.dir/src/Backends/MCMoves/CutVolume.cpp.o"
	cd /home/josh/git/alloassembly/build/contrib/romano && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT contrib/romano/CMakeFiles/CutVolume.dir/src/Backends/MCMoves/CutVolume.cpp.o -MF CMakeFiles/CutVolume.dir/src/Backends/MCMoves/CutVolume.cpp.o.d -o CMakeFiles/CutVolume.dir/src/Backends/MCMoves/CutVolume.cpp.o -c /home/josh/git/alloassembly/contrib/romano/src/Backends/MCMoves/CutVolume.cpp

contrib/romano/CMakeFiles/CutVolume.dir/src/Backends/MCMoves/CutVolume.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CutVolume.dir/src/Backends/MCMoves/CutVolume.cpp.i"
	cd /home/josh/git/alloassembly/build/contrib/romano && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/josh/git/alloassembly/contrib/romano/src/Backends/MCMoves/CutVolume.cpp > CMakeFiles/CutVolume.dir/src/Backends/MCMoves/CutVolume.cpp.i

contrib/romano/CMakeFiles/CutVolume.dir/src/Backends/MCMoves/CutVolume.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CutVolume.dir/src/Backends/MCMoves/CutVolume.cpp.s"
	cd /home/josh/git/alloassembly/build/contrib/romano && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/josh/git/alloassembly/contrib/romano/src/Backends/MCMoves/CutVolume.cpp -o CMakeFiles/CutVolume.dir/src/Backends/MCMoves/CutVolume.cpp.s

# Object files for target CutVolume
CutVolume_OBJECTS = \
"CMakeFiles/CutVolume.dir/src/Backends/MCMoves/CutVolume.cpp.o"

# External object files for target CutVolume
CutVolume_EXTERNAL_OBJECTS =

../contrib/romano/CutVolume.so: contrib/romano/CMakeFiles/CutVolume.dir/src/Backends/MCMoves/CutVolume.cpp.o
../contrib/romano/CutVolume.so: contrib/romano/CMakeFiles/CutVolume.dir/build.make
../contrib/romano/CutVolume.so: contrib/romano/CMakeFiles/CutVolume.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/josh/git/alloassembly/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared library ../../../contrib/romano/CutVolume.so"
	cd /home/josh/git/alloassembly/build/contrib/romano && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/CutVolume.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
contrib/romano/CMakeFiles/CutVolume.dir/build: ../contrib/romano/CutVolume.so
.PHONY : contrib/romano/CMakeFiles/CutVolume.dir/build

contrib/romano/CMakeFiles/CutVolume.dir/clean:
	cd /home/josh/git/alloassembly/build/contrib/romano && $(CMAKE_COMMAND) -P CMakeFiles/CutVolume.dir/cmake_clean.cmake
.PHONY : contrib/romano/CMakeFiles/CutVolume.dir/clean

contrib/romano/CMakeFiles/CutVolume.dir/depend:
	cd /home/josh/git/alloassembly/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/josh/git/alloassembly /home/josh/git/alloassembly/contrib/romano /home/josh/git/alloassembly/build /home/josh/git/alloassembly/build/contrib/romano /home/josh/git/alloassembly/build/contrib/romano/CMakeFiles/CutVolume.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : contrib/romano/CMakeFiles/CutVolume.dir/depend

