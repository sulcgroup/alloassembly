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
include contrib/rovigatti/CMakeFiles/DetailedPatchySwapInteraction.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include contrib/rovigatti/CMakeFiles/DetailedPatchySwapInteraction.dir/compiler_depend.make

# Include the progress variables for this target.
include contrib/rovigatti/CMakeFiles/DetailedPatchySwapInteraction.dir/progress.make

# Include the compile flags for this target's objects.
include contrib/rovigatti/CMakeFiles/DetailedPatchySwapInteraction.dir/flags.make

contrib/rovigatti/CMakeFiles/DetailedPatchySwapInteraction.dir/src/Interactions/DetailedPatchySwapInteraction.cpp.o: contrib/rovigatti/CMakeFiles/DetailedPatchySwapInteraction.dir/flags.make
contrib/rovigatti/CMakeFiles/DetailedPatchySwapInteraction.dir/src/Interactions/DetailedPatchySwapInteraction.cpp.o: ../contrib/rovigatti/src/Interactions/DetailedPatchySwapInteraction.cpp
contrib/rovigatti/CMakeFiles/DetailedPatchySwapInteraction.dir/src/Interactions/DetailedPatchySwapInteraction.cpp.o: contrib/rovigatti/CMakeFiles/DetailedPatchySwapInteraction.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/josh/git/alloassembly/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object contrib/rovigatti/CMakeFiles/DetailedPatchySwapInteraction.dir/src/Interactions/DetailedPatchySwapInteraction.cpp.o"
	cd /home/josh/git/alloassembly/build/contrib/rovigatti && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT contrib/rovigatti/CMakeFiles/DetailedPatchySwapInteraction.dir/src/Interactions/DetailedPatchySwapInteraction.cpp.o -MF CMakeFiles/DetailedPatchySwapInteraction.dir/src/Interactions/DetailedPatchySwapInteraction.cpp.o.d -o CMakeFiles/DetailedPatchySwapInteraction.dir/src/Interactions/DetailedPatchySwapInteraction.cpp.o -c /home/josh/git/alloassembly/contrib/rovigatti/src/Interactions/DetailedPatchySwapInteraction.cpp

contrib/rovigatti/CMakeFiles/DetailedPatchySwapInteraction.dir/src/Interactions/DetailedPatchySwapInteraction.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DetailedPatchySwapInteraction.dir/src/Interactions/DetailedPatchySwapInteraction.cpp.i"
	cd /home/josh/git/alloassembly/build/contrib/rovigatti && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/josh/git/alloassembly/contrib/rovigatti/src/Interactions/DetailedPatchySwapInteraction.cpp > CMakeFiles/DetailedPatchySwapInteraction.dir/src/Interactions/DetailedPatchySwapInteraction.cpp.i

contrib/rovigatti/CMakeFiles/DetailedPatchySwapInteraction.dir/src/Interactions/DetailedPatchySwapInteraction.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DetailedPatchySwapInteraction.dir/src/Interactions/DetailedPatchySwapInteraction.cpp.s"
	cd /home/josh/git/alloassembly/build/contrib/rovigatti && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/josh/git/alloassembly/contrib/rovigatti/src/Interactions/DetailedPatchySwapInteraction.cpp -o CMakeFiles/DetailedPatchySwapInteraction.dir/src/Interactions/DetailedPatchySwapInteraction.cpp.s

# Object files for target DetailedPatchySwapInteraction
DetailedPatchySwapInteraction_OBJECTS = \
"CMakeFiles/DetailedPatchySwapInteraction.dir/src/Interactions/DetailedPatchySwapInteraction.cpp.o"

# External object files for target DetailedPatchySwapInteraction
DetailedPatchySwapInteraction_EXTERNAL_OBJECTS =

../contrib/rovigatti/DetailedPatchySwapInteraction.so: contrib/rovigatti/CMakeFiles/DetailedPatchySwapInteraction.dir/src/Interactions/DetailedPatchySwapInteraction.cpp.o
../contrib/rovigatti/DetailedPatchySwapInteraction.so: contrib/rovigatti/CMakeFiles/DetailedPatchySwapInteraction.dir/build.make
../contrib/rovigatti/DetailedPatchySwapInteraction.so: contrib/rovigatti/CMakeFiles/DetailedPatchySwapInteraction.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/josh/git/alloassembly/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared module ../../../contrib/rovigatti/DetailedPatchySwapInteraction.so"
	cd /home/josh/git/alloassembly/build/contrib/rovigatti && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/DetailedPatchySwapInteraction.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
contrib/rovigatti/CMakeFiles/DetailedPatchySwapInteraction.dir/build: ../contrib/rovigatti/DetailedPatchySwapInteraction.so
.PHONY : contrib/rovigatti/CMakeFiles/DetailedPatchySwapInteraction.dir/build

contrib/rovigatti/CMakeFiles/DetailedPatchySwapInteraction.dir/clean:
	cd /home/josh/git/alloassembly/build/contrib/rovigatti && $(CMAKE_COMMAND) -P CMakeFiles/DetailedPatchySwapInteraction.dir/cmake_clean.cmake
.PHONY : contrib/rovigatti/CMakeFiles/DetailedPatchySwapInteraction.dir/clean

contrib/rovigatti/CMakeFiles/DetailedPatchySwapInteraction.dir/depend:
	cd /home/josh/git/alloassembly/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/josh/git/alloassembly /home/josh/git/alloassembly/contrib/rovigatti /home/josh/git/alloassembly/build /home/josh/git/alloassembly/build/contrib/rovigatti /home/josh/git/alloassembly/build/contrib/rovigatti/CMakeFiles/DetailedPatchySwapInteraction.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : contrib/rovigatti/CMakeFiles/DetailedPatchySwapInteraction.dir/depend

