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
include contrib/rovigatti/CMakeFiles/MGAnalysis.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include contrib/rovigatti/CMakeFiles/MGAnalysis.dir/compiler_depend.make

# Include the progress variables for this target.
include contrib/rovigatti/CMakeFiles/MGAnalysis.dir/progress.make

# Include the compile flags for this target's objects.
include contrib/rovigatti/CMakeFiles/MGAnalysis.dir/flags.make

contrib/rovigatti/CMakeFiles/MGAnalysis.dir/src/Observables/MGAnalysis.cpp.o: contrib/rovigatti/CMakeFiles/MGAnalysis.dir/flags.make
contrib/rovigatti/CMakeFiles/MGAnalysis.dir/src/Observables/MGAnalysis.cpp.o: ../contrib/rovigatti/src/Observables/MGAnalysis.cpp
contrib/rovigatti/CMakeFiles/MGAnalysis.dir/src/Observables/MGAnalysis.cpp.o: contrib/rovigatti/CMakeFiles/MGAnalysis.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/josh/git/alloassembly/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object contrib/rovigatti/CMakeFiles/MGAnalysis.dir/src/Observables/MGAnalysis.cpp.o"
	cd /home/josh/git/alloassembly/build/contrib/rovigatti && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT contrib/rovigatti/CMakeFiles/MGAnalysis.dir/src/Observables/MGAnalysis.cpp.o -MF CMakeFiles/MGAnalysis.dir/src/Observables/MGAnalysis.cpp.o.d -o CMakeFiles/MGAnalysis.dir/src/Observables/MGAnalysis.cpp.o -c /home/josh/git/alloassembly/contrib/rovigatti/src/Observables/MGAnalysis.cpp

contrib/rovigatti/CMakeFiles/MGAnalysis.dir/src/Observables/MGAnalysis.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/MGAnalysis.dir/src/Observables/MGAnalysis.cpp.i"
	cd /home/josh/git/alloassembly/build/contrib/rovigatti && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/josh/git/alloassembly/contrib/rovigatti/src/Observables/MGAnalysis.cpp > CMakeFiles/MGAnalysis.dir/src/Observables/MGAnalysis.cpp.i

contrib/rovigatti/CMakeFiles/MGAnalysis.dir/src/Observables/MGAnalysis.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/MGAnalysis.dir/src/Observables/MGAnalysis.cpp.s"
	cd /home/josh/git/alloassembly/build/contrib/rovigatti && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/josh/git/alloassembly/contrib/rovigatti/src/Observables/MGAnalysis.cpp -o CMakeFiles/MGAnalysis.dir/src/Observables/MGAnalysis.cpp.s

# Object files for target MGAnalysis
MGAnalysis_OBJECTS = \
"CMakeFiles/MGAnalysis.dir/src/Observables/MGAnalysis.cpp.o"

# External object files for target MGAnalysis
MGAnalysis_EXTERNAL_OBJECTS =

../contrib/rovigatti/MGAnalysis.so: contrib/rovigatti/CMakeFiles/MGAnalysis.dir/src/Observables/MGAnalysis.cpp.o
../contrib/rovigatti/MGAnalysis.so: contrib/rovigatti/CMakeFiles/MGAnalysis.dir/build.make
../contrib/rovigatti/MGAnalysis.so: contrib/rovigatti/CMakeFiles/MGAnalysis.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/josh/git/alloassembly/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared module ../../../contrib/rovigatti/MGAnalysis.so"
	cd /home/josh/git/alloassembly/build/contrib/rovigatti && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/MGAnalysis.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
contrib/rovigatti/CMakeFiles/MGAnalysis.dir/build: ../contrib/rovigatti/MGAnalysis.so
.PHONY : contrib/rovigatti/CMakeFiles/MGAnalysis.dir/build

contrib/rovigatti/CMakeFiles/MGAnalysis.dir/clean:
	cd /home/josh/git/alloassembly/build/contrib/rovigatti && $(CMAKE_COMMAND) -P CMakeFiles/MGAnalysis.dir/cmake_clean.cmake
.PHONY : contrib/rovigatti/CMakeFiles/MGAnalysis.dir/clean

contrib/rovigatti/CMakeFiles/MGAnalysis.dir/depend:
	cd /home/josh/git/alloassembly/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/josh/git/alloassembly /home/josh/git/alloassembly/contrib/rovigatti /home/josh/git/alloassembly/build /home/josh/git/alloassembly/build/contrib/rovigatti /home/josh/git/alloassembly/build/contrib/rovigatti/CMakeFiles/MGAnalysis.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : contrib/rovigatti/CMakeFiles/MGAnalysis.dir/depend

