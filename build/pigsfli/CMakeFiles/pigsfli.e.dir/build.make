# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.19

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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.19.1/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.19.1/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/ecasiano/Desktop/pigsfli/src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/ecasiano/Desktop/pigsfli/build

# Include any dependencies generated for this target.
include pigsfli/CMakeFiles/pigsfli.e.dir/depend.make

# Include the progress variables for this target.
include pigsfli/CMakeFiles/pigsfli.e.dir/progress.make

# Include the compile flags for this target's objects.
include pigsfli/CMakeFiles/pigsfli.e.dir/flags.make

pigsfli/CMakeFiles/pigsfli.e.dir/src/RNG.cpp.o: pigsfli/CMakeFiles/pigsfli.e.dir/flags.make
pigsfli/CMakeFiles/pigsfli.e.dir/src/RNG.cpp.o: /Users/ecasiano/Desktop/pigsfli/src/pigsfli/src/RNG.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/ecasiano/Desktop/pigsfli/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object pigsfli/CMakeFiles/pigsfli.e.dir/src/RNG.cpp.o"
	cd /Users/ecasiano/Desktop/pigsfli/build/pigsfli && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/pigsfli.e.dir/src/RNG.cpp.o -c /Users/ecasiano/Desktop/pigsfli/src/pigsfli/src/RNG.cpp

pigsfli/CMakeFiles/pigsfli.e.dir/src/RNG.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/pigsfli.e.dir/src/RNG.cpp.i"
	cd /Users/ecasiano/Desktop/pigsfli/build/pigsfli && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/ecasiano/Desktop/pigsfli/src/pigsfli/src/RNG.cpp > CMakeFiles/pigsfli.e.dir/src/RNG.cpp.i

pigsfli/CMakeFiles/pigsfli.e.dir/src/RNG.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/pigsfli.e.dir/src/RNG.cpp.s"
	cd /Users/ecasiano/Desktop/pigsfli/build/pigsfli && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/ecasiano/Desktop/pigsfli/src/pigsfli/src/RNG.cpp -o CMakeFiles/pigsfli.e.dir/src/RNG.cpp.s

pigsfli/CMakeFiles/pigsfli.e.dir/src/pimc.cpp.o: pigsfli/CMakeFiles/pigsfli.e.dir/flags.make
pigsfli/CMakeFiles/pigsfli.e.dir/src/pimc.cpp.o: /Users/ecasiano/Desktop/pigsfli/src/pigsfli/src/pimc.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/ecasiano/Desktop/pigsfli/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object pigsfli/CMakeFiles/pigsfli.e.dir/src/pimc.cpp.o"
	cd /Users/ecasiano/Desktop/pigsfli/build/pigsfli && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/pigsfli.e.dir/src/pimc.cpp.o -c /Users/ecasiano/Desktop/pigsfli/src/pigsfli/src/pimc.cpp

pigsfli/CMakeFiles/pigsfli.e.dir/src/pimc.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/pigsfli.e.dir/src/pimc.cpp.i"
	cd /Users/ecasiano/Desktop/pigsfli/build/pigsfli && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/ecasiano/Desktop/pigsfli/src/pigsfli/src/pimc.cpp > CMakeFiles/pigsfli.e.dir/src/pimc.cpp.i

pigsfli/CMakeFiles/pigsfli.e.dir/src/pimc.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/pigsfli.e.dir/src/pimc.cpp.s"
	cd /Users/ecasiano/Desktop/pigsfli/build/pigsfli && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/ecasiano/Desktop/pigsfli/src/pigsfli/src/pimc.cpp -o CMakeFiles/pigsfli.e.dir/src/pimc.cpp.s

# Object files for target pigsfli.e
pigsfli_e_OBJECTS = \
"CMakeFiles/pigsfli.e.dir/src/RNG.cpp.o" \
"CMakeFiles/pigsfli.e.dir/src/pimc.cpp.o"

# External object files for target pigsfli.e
pigsfli_e_EXTERNAL_OBJECTS =

pigsfli/pigsfli.e: pigsfli/CMakeFiles/pigsfli.e.dir/src/RNG.cpp.o
pigsfli/pigsfli.e: pigsfli/CMakeFiles/pigsfli.e.dir/src/pimc.cpp.o
pigsfli/pigsfli.e: pigsfli/CMakeFiles/pigsfli.e.dir/build.make
pigsfli/pigsfli.e: pigsfli/CMakeFiles/pigsfli.e.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/ecasiano/Desktop/pigsfli/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable pigsfli.e"
	cd /Users/ecasiano/Desktop/pigsfli/build/pigsfli && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/pigsfli.e.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
pigsfli/CMakeFiles/pigsfli.e.dir/build: pigsfli/pigsfli.e

.PHONY : pigsfli/CMakeFiles/pigsfli.e.dir/build

pigsfli/CMakeFiles/pigsfli.e.dir/clean:
	cd /Users/ecasiano/Desktop/pigsfli/build/pigsfli && $(CMAKE_COMMAND) -P CMakeFiles/pigsfli.e.dir/cmake_clean.cmake
.PHONY : pigsfli/CMakeFiles/pigsfli.e.dir/clean

pigsfli/CMakeFiles/pigsfli.e.dir/depend:
	cd /Users/ecasiano/Desktop/pigsfli/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/ecasiano/Desktop/pigsfli/src /Users/ecasiano/Desktop/pigsfli/src/pigsfli /Users/ecasiano/Desktop/pigsfli/build /Users/ecasiano/Desktop/pigsfli/build/pigsfli /Users/ecasiano/Desktop/pigsfli/build/pigsfli/CMakeFiles/pigsfli.e.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : pigsfli/CMakeFiles/pigsfli.e.dir/depend

