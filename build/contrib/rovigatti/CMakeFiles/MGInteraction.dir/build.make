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
include contrib/rovigatti/CMakeFiles/MGInteraction.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include contrib/rovigatti/CMakeFiles/MGInteraction.dir/compiler_depend.make

# Include the progress variables for this target.
include contrib/rovigatti/CMakeFiles/MGInteraction.dir/progress.make

# Include the compile flags for this target's objects.
include contrib/rovigatti/CMakeFiles/MGInteraction.dir/flags.make

contrib/rovigatti/CMakeFiles/MGInteraction.dir/src/Interactions/MGInteraction.cpp.o: contrib/rovigatti/CMakeFiles/MGInteraction.dir/flags.make
contrib/rovigatti/CMakeFiles/MGInteraction.dir/src/Interactions/MGInteraction.cpp.o: ../contrib/rovigatti/src/Interactions/MGInteraction.cpp
contrib/rovigatti/CMakeFiles/MGInteraction.dir/src/Interactions/MGInteraction.cpp.o: contrib/rovigatti/CMakeFiles/MGInteraction.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/josh/git/alloassembly/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object contrib/rovigatti/CMakeFiles/MGInteraction.dir/src/Interactions/MGInteraction.cpp.o"
	cd /home/josh/git/alloassembly/build/contrib/rovigatti && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT contrib/rovigatti/CMakeFiles/MGInteraction.dir/src/Interactions/MGInteraction.cpp.o -MF CMakeFiles/MGInteraction.dir/src/Interactions/MGInteraction.cpp.o.d -o CMakeFiles/MGInteraction.dir/src/Interactions/MGInteraction.cpp.o -c /home/josh/git/alloassembly/contrib/rovigatti/src/Interactions/MGInteraction.cpp

contrib/rovigatti/CMakeFiles/MGInteraction.dir/src/Interactions/MGInteraction.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/MGInteraction.dir/src/Interactions/MGInteraction.cpp.i"
	cd /home/josh/git/alloassembly/build/contrib/rovigatti && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/josh/git/alloassembly/contrib/rovigatti/src/Interactions/MGInteraction.cpp > CMakeFiles/MGInteraction.dir/src/Interactions/MGInteraction.cpp.i

contrib/rovigatti/CMakeFiles/MGInteraction.dir/src/Interactions/MGInteraction.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/MGInteraction.dir/src/Interactions/MGInteraction.cpp.s"
	cd /home/josh/git/alloassembly/build/contrib/rovigatti && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/josh/git/alloassembly/contrib/rovigatti/src/Interactions/MGInteraction.cpp -o CMakeFiles/MGInteraction.dir/src/Interactions/MGInteraction.cpp.s

# Object files for target MGInteraction
MGInteraction_OBJECTS = \
"CMakeFiles/MGInteraction.dir/src/Interactions/MGInteraction.cpp.o"

# External object files for target MGInteraction
MGInteraction_EXTERNAL_OBJECTS =

../contrib/rovigatti/MGInteraction.so: contrib/rovigatti/CMakeFiles/MGInteraction.dir/src/Interactions/MGInteraction.cpp.o
../contrib/rovigatti/MGInteraction.so: contrib/rovigatti/CMakeFiles/MGInteraction.dir/build.make
../contrib/rovigatti/MGInteraction.so: contrib/rovigatti/CMakeFiles/MGInteraction.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/josh/git/alloassembly/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared module ../../../contrib/rovigatti/MGInteraction.so"
	cd /home/josh/git/alloassembly/build/contrib/rovigatti && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/MGInteraction.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
contrib/rovigatti/CMakeFiles/MGInteraction.dir/build: ../contrib/rovigatti/MGInteraction.so
.PHONY : contrib/rovigatti/CMakeFiles/MGInteraction.dir/build

contrib/rovigatti/CMakeFiles/MGInteraction.dir/clean:
	cd /home/josh/git/alloassembly/build/contrib/rovigatti && $(CMAKE_COMMAND) -P CMakeFiles/MGInteraction.dir/cmake_clean.cmake
.PHONY : contrib/rovigatti/CMakeFiles/MGInteraction.dir/clean

contrib/rovigatti/CMakeFiles/MGInteraction.dir/depend:
	cd /home/josh/git/alloassembly/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/josh/git/alloassembly /home/josh/git/alloassembly/contrib/rovigatti /home/josh/git/alloassembly/build /home/josh/git/alloassembly/build/contrib/rovigatti /home/josh/git/alloassembly/build/contrib/rovigatti/CMakeFiles/MGInteraction.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : contrib/rovigatti/CMakeFiles/MGInteraction.dir/depend

