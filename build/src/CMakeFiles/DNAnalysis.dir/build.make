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
include src/CMakeFiles/DNAnalysis.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/CMakeFiles/DNAnalysis.dir/compiler_depend.make

# Include the progress variables for this target.
include src/CMakeFiles/DNAnalysis.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/DNAnalysis.dir/flags.make

src/CMakeFiles/DNAnalysis.dir/DNAnalysis.cpp.o: src/CMakeFiles/DNAnalysis.dir/flags.make
src/CMakeFiles/DNAnalysis.dir/DNAnalysis.cpp.o: ../src/DNAnalysis.cpp
src/CMakeFiles/DNAnalysis.dir/DNAnalysis.cpp.o: src/CMakeFiles/DNAnalysis.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/josh/git/alloassembly/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/DNAnalysis.dir/DNAnalysis.cpp.o"
	cd /home/josh/git/alloassembly/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/CMakeFiles/DNAnalysis.dir/DNAnalysis.cpp.o -MF CMakeFiles/DNAnalysis.dir/DNAnalysis.cpp.o.d -o CMakeFiles/DNAnalysis.dir/DNAnalysis.cpp.o -c /home/josh/git/alloassembly/src/DNAnalysis.cpp

src/CMakeFiles/DNAnalysis.dir/DNAnalysis.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DNAnalysis.dir/DNAnalysis.cpp.i"
	cd /home/josh/git/alloassembly/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/josh/git/alloassembly/src/DNAnalysis.cpp > CMakeFiles/DNAnalysis.dir/DNAnalysis.cpp.i

src/CMakeFiles/DNAnalysis.dir/DNAnalysis.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DNAnalysis.dir/DNAnalysis.cpp.s"
	cd /home/josh/git/alloassembly/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/josh/git/alloassembly/src/DNAnalysis.cpp -o CMakeFiles/DNAnalysis.dir/DNAnalysis.cpp.s

src/CMakeFiles/DNAnalysis.dir/Managers/AnalysisManager.cpp.o: src/CMakeFiles/DNAnalysis.dir/flags.make
src/CMakeFiles/DNAnalysis.dir/Managers/AnalysisManager.cpp.o: ../src/Managers/AnalysisManager.cpp
src/CMakeFiles/DNAnalysis.dir/Managers/AnalysisManager.cpp.o: src/CMakeFiles/DNAnalysis.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/josh/git/alloassembly/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object src/CMakeFiles/DNAnalysis.dir/Managers/AnalysisManager.cpp.o"
	cd /home/josh/git/alloassembly/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/CMakeFiles/DNAnalysis.dir/Managers/AnalysisManager.cpp.o -MF CMakeFiles/DNAnalysis.dir/Managers/AnalysisManager.cpp.o.d -o CMakeFiles/DNAnalysis.dir/Managers/AnalysisManager.cpp.o -c /home/josh/git/alloassembly/src/Managers/AnalysisManager.cpp

src/CMakeFiles/DNAnalysis.dir/Managers/AnalysisManager.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DNAnalysis.dir/Managers/AnalysisManager.cpp.i"
	cd /home/josh/git/alloassembly/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/josh/git/alloassembly/src/Managers/AnalysisManager.cpp > CMakeFiles/DNAnalysis.dir/Managers/AnalysisManager.cpp.i

src/CMakeFiles/DNAnalysis.dir/Managers/AnalysisManager.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DNAnalysis.dir/Managers/AnalysisManager.cpp.s"
	cd /home/josh/git/alloassembly/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/josh/git/alloassembly/src/Managers/AnalysisManager.cpp -o CMakeFiles/DNAnalysis.dir/Managers/AnalysisManager.cpp.s

src/CMakeFiles/DNAnalysis.dir/Utilities/Timings.cpp.o: src/CMakeFiles/DNAnalysis.dir/flags.make
src/CMakeFiles/DNAnalysis.dir/Utilities/Timings.cpp.o: ../src/Utilities/Timings.cpp
src/CMakeFiles/DNAnalysis.dir/Utilities/Timings.cpp.o: src/CMakeFiles/DNAnalysis.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/josh/git/alloassembly/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object src/CMakeFiles/DNAnalysis.dir/Utilities/Timings.cpp.o"
	cd /home/josh/git/alloassembly/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/CMakeFiles/DNAnalysis.dir/Utilities/Timings.cpp.o -MF CMakeFiles/DNAnalysis.dir/Utilities/Timings.cpp.o.d -o CMakeFiles/DNAnalysis.dir/Utilities/Timings.cpp.o -c /home/josh/git/alloassembly/src/Utilities/Timings.cpp

src/CMakeFiles/DNAnalysis.dir/Utilities/Timings.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DNAnalysis.dir/Utilities/Timings.cpp.i"
	cd /home/josh/git/alloassembly/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/josh/git/alloassembly/src/Utilities/Timings.cpp > CMakeFiles/DNAnalysis.dir/Utilities/Timings.cpp.i

src/CMakeFiles/DNAnalysis.dir/Utilities/Timings.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DNAnalysis.dir/Utilities/Timings.cpp.s"
	cd /home/josh/git/alloassembly/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/josh/git/alloassembly/src/Utilities/Timings.cpp -o CMakeFiles/DNAnalysis.dir/Utilities/Timings.cpp.s

# Object files for target DNAnalysis
DNAnalysis_OBJECTS = \
"CMakeFiles/DNAnalysis.dir/DNAnalysis.cpp.o" \
"CMakeFiles/DNAnalysis.dir/Managers/AnalysisManager.cpp.o" \
"CMakeFiles/DNAnalysis.dir/Utilities/Timings.cpp.o"

# External object files for target DNAnalysis
DNAnalysis_EXTERNAL_OBJECTS =

bin/DNAnalysis: src/CMakeFiles/DNAnalysis.dir/DNAnalysis.cpp.o
bin/DNAnalysis: src/CMakeFiles/DNAnalysis.dir/Managers/AnalysisManager.cpp.o
bin/DNAnalysis: src/CMakeFiles/DNAnalysis.dir/Utilities/Timings.cpp.o
bin/DNAnalysis: src/CMakeFiles/DNAnalysis.dir/build.make
bin/DNAnalysis: /usr/local/cuda/lib64/libcudart_static.a
bin/DNAnalysis: /usr/lib/x86_64-linux-gnu/librt.so
bin/DNAnalysis: src/libcommon_debug.so
bin/DNAnalysis: /usr/local/cuda/lib64/libcudart_static.a
bin/DNAnalysis: /usr/lib/x86_64-linux-gnu/librt.so
bin/DNAnalysis: src/CMakeFiles/DNAnalysis.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/josh/git/alloassembly/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX executable ../bin/DNAnalysis"
	cd /home/josh/git/alloassembly/build/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/DNAnalysis.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/DNAnalysis.dir/build: bin/DNAnalysis
.PHONY : src/CMakeFiles/DNAnalysis.dir/build

src/CMakeFiles/DNAnalysis.dir/clean:
	cd /home/josh/git/alloassembly/build/src && $(CMAKE_COMMAND) -P CMakeFiles/DNAnalysis.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/DNAnalysis.dir/clean

src/CMakeFiles/DNAnalysis.dir/depend:
	cd /home/josh/git/alloassembly/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/josh/git/alloassembly /home/josh/git/alloassembly/src /home/josh/git/alloassembly/build /home/josh/git/alloassembly/build/src /home/josh/git/alloassembly/build/src/CMakeFiles/DNAnalysis.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/DNAnalysis.dir/depend

