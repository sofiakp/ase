# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canoncical targets will work.
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

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /media/fusion10/work/chromatinVariation/src/ase_cpp

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /media/fusion10/work/chromatinVariation/src/ase_cpp/build/debug

# Include any dependencies generated for this target.
include deps/swak/test/CMakeFiles/TestSwak.dir/depend.make

# Include the progress variables for this target.
include deps/swak/test/CMakeFiles/TestSwak.dir/progress.make

# Include the compile flags for this target's objects.
include deps/swak/test/CMakeFiles/TestSwak.dir/flags.make

deps/swak/test/CMakeFiles/TestSwak.dir/TestSwak.cpp.o: deps/swak/test/CMakeFiles/TestSwak.dir/flags.make
deps/swak/test/CMakeFiles/TestSwak.dir/TestSwak.cpp.o: ../../deps/swak/test/TestSwak.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /media/fusion10/work/chromatinVariation/src/ase_cpp/build/debug/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object deps/swak/test/CMakeFiles/TestSwak.dir/TestSwak.cpp.o"
	cd /media/fusion10/work/chromatinVariation/src/ase_cpp/build/debug/deps/swak/test && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/TestSwak.dir/TestSwak.cpp.o -c /media/fusion10/work/chromatinVariation/src/ase_cpp/deps/swak/test/TestSwak.cpp

deps/swak/test/CMakeFiles/TestSwak.dir/TestSwak.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/TestSwak.dir/TestSwak.cpp.i"
	cd /media/fusion10/work/chromatinVariation/src/ase_cpp/build/debug/deps/swak/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /media/fusion10/work/chromatinVariation/src/ase_cpp/deps/swak/test/TestSwak.cpp > CMakeFiles/TestSwak.dir/TestSwak.cpp.i

deps/swak/test/CMakeFiles/TestSwak.dir/TestSwak.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/TestSwak.dir/TestSwak.cpp.s"
	cd /media/fusion10/work/chromatinVariation/src/ase_cpp/build/debug/deps/swak/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /media/fusion10/work/chromatinVariation/src/ase_cpp/deps/swak/test/TestSwak.cpp -o CMakeFiles/TestSwak.dir/TestSwak.cpp.s

deps/swak/test/CMakeFiles/TestSwak.dir/TestSwak.cpp.o.requires:
.PHONY : deps/swak/test/CMakeFiles/TestSwak.dir/TestSwak.cpp.o.requires

deps/swak/test/CMakeFiles/TestSwak.dir/TestSwak.cpp.o.provides: deps/swak/test/CMakeFiles/TestSwak.dir/TestSwak.cpp.o.requires
	$(MAKE) -f deps/swak/test/CMakeFiles/TestSwak.dir/build.make deps/swak/test/CMakeFiles/TestSwak.dir/TestSwak.cpp.o.provides.build
.PHONY : deps/swak/test/CMakeFiles/TestSwak.dir/TestSwak.cpp.o.provides

deps/swak/test/CMakeFiles/TestSwak.dir/TestSwak.cpp.o.provides.build: deps/swak/test/CMakeFiles/TestSwak.dir/TestSwak.cpp.o
.PHONY : deps/swak/test/CMakeFiles/TestSwak.dir/TestSwak.cpp.o.provides.build

deps/swak/test/CMakeFiles/TestSwak.dir/TestHelpers.cpp.o: deps/swak/test/CMakeFiles/TestSwak.dir/flags.make
deps/swak/test/CMakeFiles/TestSwak.dir/TestHelpers.cpp.o: ../../deps/swak/test/TestHelpers.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /media/fusion10/work/chromatinVariation/src/ase_cpp/build/debug/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object deps/swak/test/CMakeFiles/TestSwak.dir/TestHelpers.cpp.o"
	cd /media/fusion10/work/chromatinVariation/src/ase_cpp/build/debug/deps/swak/test && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/TestSwak.dir/TestHelpers.cpp.o -c /media/fusion10/work/chromatinVariation/src/ase_cpp/deps/swak/test/TestHelpers.cpp

deps/swak/test/CMakeFiles/TestSwak.dir/TestHelpers.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/TestSwak.dir/TestHelpers.cpp.i"
	cd /media/fusion10/work/chromatinVariation/src/ase_cpp/build/debug/deps/swak/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /media/fusion10/work/chromatinVariation/src/ase_cpp/deps/swak/test/TestHelpers.cpp > CMakeFiles/TestSwak.dir/TestHelpers.cpp.i

deps/swak/test/CMakeFiles/TestSwak.dir/TestHelpers.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/TestSwak.dir/TestHelpers.cpp.s"
	cd /media/fusion10/work/chromatinVariation/src/ase_cpp/build/debug/deps/swak/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /media/fusion10/work/chromatinVariation/src/ase_cpp/deps/swak/test/TestHelpers.cpp -o CMakeFiles/TestSwak.dir/TestHelpers.cpp.s

deps/swak/test/CMakeFiles/TestSwak.dir/TestHelpers.cpp.o.requires:
.PHONY : deps/swak/test/CMakeFiles/TestSwak.dir/TestHelpers.cpp.o.requires

deps/swak/test/CMakeFiles/TestSwak.dir/TestHelpers.cpp.o.provides: deps/swak/test/CMakeFiles/TestSwak.dir/TestHelpers.cpp.o.requires
	$(MAKE) -f deps/swak/test/CMakeFiles/TestSwak.dir/build.make deps/swak/test/CMakeFiles/TestSwak.dir/TestHelpers.cpp.o.provides.build
.PHONY : deps/swak/test/CMakeFiles/TestSwak.dir/TestHelpers.cpp.o.provides

deps/swak/test/CMakeFiles/TestSwak.dir/TestHelpers.cpp.o.provides.build: deps/swak/test/CMakeFiles/TestSwak.dir/TestHelpers.cpp.o
.PHONY : deps/swak/test/CMakeFiles/TestSwak.dir/TestHelpers.cpp.o.provides.build

deps/swak/test/CMakeFiles/TestSwak.dir/TestSystem.cpp.o: deps/swak/test/CMakeFiles/TestSwak.dir/flags.make
deps/swak/test/CMakeFiles/TestSwak.dir/TestSystem.cpp.o: ../../deps/swak/test/TestSystem.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /media/fusion10/work/chromatinVariation/src/ase_cpp/build/debug/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object deps/swak/test/CMakeFiles/TestSwak.dir/TestSystem.cpp.o"
	cd /media/fusion10/work/chromatinVariation/src/ase_cpp/build/debug/deps/swak/test && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/TestSwak.dir/TestSystem.cpp.o -c /media/fusion10/work/chromatinVariation/src/ase_cpp/deps/swak/test/TestSystem.cpp

deps/swak/test/CMakeFiles/TestSwak.dir/TestSystem.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/TestSwak.dir/TestSystem.cpp.i"
	cd /media/fusion10/work/chromatinVariation/src/ase_cpp/build/debug/deps/swak/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /media/fusion10/work/chromatinVariation/src/ase_cpp/deps/swak/test/TestSystem.cpp > CMakeFiles/TestSwak.dir/TestSystem.cpp.i

deps/swak/test/CMakeFiles/TestSwak.dir/TestSystem.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/TestSwak.dir/TestSystem.cpp.s"
	cd /media/fusion10/work/chromatinVariation/src/ase_cpp/build/debug/deps/swak/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /media/fusion10/work/chromatinVariation/src/ase_cpp/deps/swak/test/TestSystem.cpp -o CMakeFiles/TestSwak.dir/TestSystem.cpp.s

deps/swak/test/CMakeFiles/TestSwak.dir/TestSystem.cpp.o.requires:
.PHONY : deps/swak/test/CMakeFiles/TestSwak.dir/TestSystem.cpp.o.requires

deps/swak/test/CMakeFiles/TestSwak.dir/TestSystem.cpp.o.provides: deps/swak/test/CMakeFiles/TestSwak.dir/TestSystem.cpp.o.requires
	$(MAKE) -f deps/swak/test/CMakeFiles/TestSwak.dir/build.make deps/swak/test/CMakeFiles/TestSwak.dir/TestSystem.cpp.o.provides.build
.PHONY : deps/swak/test/CMakeFiles/TestSwak.dir/TestSystem.cpp.o.provides

deps/swak/test/CMakeFiles/TestSwak.dir/TestSystem.cpp.o.provides.build: deps/swak/test/CMakeFiles/TestSwak.dir/TestSystem.cpp.o
.PHONY : deps/swak/test/CMakeFiles/TestSwak.dir/TestSystem.cpp.o.provides.build

deps/swak/test/CMakeFiles/TestSwak.dir/TestSafeVec.cpp.o: deps/swak/test/CMakeFiles/TestSwak.dir/flags.make
deps/swak/test/CMakeFiles/TestSwak.dir/TestSafeVec.cpp.o: ../../deps/swak/test/TestSafeVec.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /media/fusion10/work/chromatinVariation/src/ase_cpp/build/debug/CMakeFiles $(CMAKE_PROGRESS_4)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object deps/swak/test/CMakeFiles/TestSwak.dir/TestSafeVec.cpp.o"
	cd /media/fusion10/work/chromatinVariation/src/ase_cpp/build/debug/deps/swak/test && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/TestSwak.dir/TestSafeVec.cpp.o -c /media/fusion10/work/chromatinVariation/src/ase_cpp/deps/swak/test/TestSafeVec.cpp

deps/swak/test/CMakeFiles/TestSwak.dir/TestSafeVec.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/TestSwak.dir/TestSafeVec.cpp.i"
	cd /media/fusion10/work/chromatinVariation/src/ase_cpp/build/debug/deps/swak/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /media/fusion10/work/chromatinVariation/src/ase_cpp/deps/swak/test/TestSafeVec.cpp > CMakeFiles/TestSwak.dir/TestSafeVec.cpp.i

deps/swak/test/CMakeFiles/TestSwak.dir/TestSafeVec.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/TestSwak.dir/TestSafeVec.cpp.s"
	cd /media/fusion10/work/chromatinVariation/src/ase_cpp/build/debug/deps/swak/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /media/fusion10/work/chromatinVariation/src/ase_cpp/deps/swak/test/TestSafeVec.cpp -o CMakeFiles/TestSwak.dir/TestSafeVec.cpp.s

deps/swak/test/CMakeFiles/TestSwak.dir/TestSafeVec.cpp.o.requires:
.PHONY : deps/swak/test/CMakeFiles/TestSwak.dir/TestSafeVec.cpp.o.requires

deps/swak/test/CMakeFiles/TestSwak.dir/TestSafeVec.cpp.o.provides: deps/swak/test/CMakeFiles/TestSwak.dir/TestSafeVec.cpp.o.requires
	$(MAKE) -f deps/swak/test/CMakeFiles/TestSwak.dir/build.make deps/swak/test/CMakeFiles/TestSwak.dir/TestSafeVec.cpp.o.provides.build
.PHONY : deps/swak/test/CMakeFiles/TestSwak.dir/TestSafeVec.cpp.o.provides

deps/swak/test/CMakeFiles/TestSwak.dir/TestSafeVec.cpp.o.provides.build: deps/swak/test/CMakeFiles/TestSwak.dir/TestSafeVec.cpp.o
.PHONY : deps/swak/test/CMakeFiles/TestSwak.dir/TestSafeVec.cpp.o.provides.build

deps/swak/test/CMakeFiles/TestSwak.dir/TestExtractDigits.cpp.o: deps/swak/test/CMakeFiles/TestSwak.dir/flags.make
deps/swak/test/CMakeFiles/TestSwak.dir/TestExtractDigits.cpp.o: ../../deps/swak/test/TestExtractDigits.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /media/fusion10/work/chromatinVariation/src/ase_cpp/build/debug/CMakeFiles $(CMAKE_PROGRESS_5)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object deps/swak/test/CMakeFiles/TestSwak.dir/TestExtractDigits.cpp.o"
	cd /media/fusion10/work/chromatinVariation/src/ase_cpp/build/debug/deps/swak/test && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/TestSwak.dir/TestExtractDigits.cpp.o -c /media/fusion10/work/chromatinVariation/src/ase_cpp/deps/swak/test/TestExtractDigits.cpp

deps/swak/test/CMakeFiles/TestSwak.dir/TestExtractDigits.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/TestSwak.dir/TestExtractDigits.cpp.i"
	cd /media/fusion10/work/chromatinVariation/src/ase_cpp/build/debug/deps/swak/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /media/fusion10/work/chromatinVariation/src/ase_cpp/deps/swak/test/TestExtractDigits.cpp > CMakeFiles/TestSwak.dir/TestExtractDigits.cpp.i

deps/swak/test/CMakeFiles/TestSwak.dir/TestExtractDigits.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/TestSwak.dir/TestExtractDigits.cpp.s"
	cd /media/fusion10/work/chromatinVariation/src/ase_cpp/build/debug/deps/swak/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /media/fusion10/work/chromatinVariation/src/ase_cpp/deps/swak/test/TestExtractDigits.cpp -o CMakeFiles/TestSwak.dir/TestExtractDigits.cpp.s

deps/swak/test/CMakeFiles/TestSwak.dir/TestExtractDigits.cpp.o.requires:
.PHONY : deps/swak/test/CMakeFiles/TestSwak.dir/TestExtractDigits.cpp.o.requires

deps/swak/test/CMakeFiles/TestSwak.dir/TestExtractDigits.cpp.o.provides: deps/swak/test/CMakeFiles/TestSwak.dir/TestExtractDigits.cpp.o.requires
	$(MAKE) -f deps/swak/test/CMakeFiles/TestSwak.dir/build.make deps/swak/test/CMakeFiles/TestSwak.dir/TestExtractDigits.cpp.o.provides.build
.PHONY : deps/swak/test/CMakeFiles/TestSwak.dir/TestExtractDigits.cpp.o.provides

deps/swak/test/CMakeFiles/TestSwak.dir/TestExtractDigits.cpp.o.provides.build: deps/swak/test/CMakeFiles/TestSwak.dir/TestExtractDigits.cpp.o
.PHONY : deps/swak/test/CMakeFiles/TestSwak.dir/TestExtractDigits.cpp.o.provides.build

deps/swak/test/CMakeFiles/TestSwak.dir/TestBinaryIO.cpp.o: deps/swak/test/CMakeFiles/TestSwak.dir/flags.make
deps/swak/test/CMakeFiles/TestSwak.dir/TestBinaryIO.cpp.o: ../../deps/swak/test/TestBinaryIO.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /media/fusion10/work/chromatinVariation/src/ase_cpp/build/debug/CMakeFiles $(CMAKE_PROGRESS_6)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object deps/swak/test/CMakeFiles/TestSwak.dir/TestBinaryIO.cpp.o"
	cd /media/fusion10/work/chromatinVariation/src/ase_cpp/build/debug/deps/swak/test && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/TestSwak.dir/TestBinaryIO.cpp.o -c /media/fusion10/work/chromatinVariation/src/ase_cpp/deps/swak/test/TestBinaryIO.cpp

deps/swak/test/CMakeFiles/TestSwak.dir/TestBinaryIO.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/TestSwak.dir/TestBinaryIO.cpp.i"
	cd /media/fusion10/work/chromatinVariation/src/ase_cpp/build/debug/deps/swak/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /media/fusion10/work/chromatinVariation/src/ase_cpp/deps/swak/test/TestBinaryIO.cpp > CMakeFiles/TestSwak.dir/TestBinaryIO.cpp.i

deps/swak/test/CMakeFiles/TestSwak.dir/TestBinaryIO.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/TestSwak.dir/TestBinaryIO.cpp.s"
	cd /media/fusion10/work/chromatinVariation/src/ase_cpp/build/debug/deps/swak/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /media/fusion10/work/chromatinVariation/src/ase_cpp/deps/swak/test/TestBinaryIO.cpp -o CMakeFiles/TestSwak.dir/TestBinaryIO.cpp.s

deps/swak/test/CMakeFiles/TestSwak.dir/TestBinaryIO.cpp.o.requires:
.PHONY : deps/swak/test/CMakeFiles/TestSwak.dir/TestBinaryIO.cpp.o.requires

deps/swak/test/CMakeFiles/TestSwak.dir/TestBinaryIO.cpp.o.provides: deps/swak/test/CMakeFiles/TestSwak.dir/TestBinaryIO.cpp.o.requires
	$(MAKE) -f deps/swak/test/CMakeFiles/TestSwak.dir/build.make deps/swak/test/CMakeFiles/TestSwak.dir/TestBinaryIO.cpp.o.provides.build
.PHONY : deps/swak/test/CMakeFiles/TestSwak.dir/TestBinaryIO.cpp.o.provides

deps/swak/test/CMakeFiles/TestSwak.dir/TestBinaryIO.cpp.o.provides.build: deps/swak/test/CMakeFiles/TestSwak.dir/TestBinaryIO.cpp.o
.PHONY : deps/swak/test/CMakeFiles/TestSwak.dir/TestBinaryIO.cpp.o.provides.build

# Object files for target TestSwak
TestSwak_OBJECTS = \
"CMakeFiles/TestSwak.dir/TestSwak.cpp.o" \
"CMakeFiles/TestSwak.dir/TestHelpers.cpp.o" \
"CMakeFiles/TestSwak.dir/TestSystem.cpp.o" \
"CMakeFiles/TestSwak.dir/TestSafeVec.cpp.o" \
"CMakeFiles/TestSwak.dir/TestExtractDigits.cpp.o" \
"CMakeFiles/TestSwak.dir/TestBinaryIO.cpp.o"

# External object files for target TestSwak
TestSwak_EXTERNAL_OBJECTS =

deps/swak/test/TestSwak: deps/swak/test/CMakeFiles/TestSwak.dir/TestSwak.cpp.o
deps/swak/test/TestSwak: deps/swak/test/CMakeFiles/TestSwak.dir/TestHelpers.cpp.o
deps/swak/test/TestSwak: deps/swak/test/CMakeFiles/TestSwak.dir/TestSystem.cpp.o
deps/swak/test/TestSwak: deps/swak/test/CMakeFiles/TestSwak.dir/TestSafeVec.cpp.o
deps/swak/test/TestSwak: deps/swak/test/CMakeFiles/TestSwak.dir/TestExtractDigits.cpp.o
deps/swak/test/TestSwak: deps/swak/test/CMakeFiles/TestSwak.dir/TestBinaryIO.cpp.o
deps/swak/test/TestSwak: deps/swak/src/libswak.a
deps/swak/test/TestSwak: deps/swak/deps/yaml-cpp/libyaml-cpp.a
deps/swak/test/TestSwak: deps/swak/test/CMakeFiles/TestSwak.dir/build.make
deps/swak/test/TestSwak: deps/swak/test/CMakeFiles/TestSwak.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable TestSwak"
	cd /media/fusion10/work/chromatinVariation/src/ase_cpp/build/debug/deps/swak/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/TestSwak.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
deps/swak/test/CMakeFiles/TestSwak.dir/build: deps/swak/test/TestSwak
.PHONY : deps/swak/test/CMakeFiles/TestSwak.dir/build

deps/swak/test/CMakeFiles/TestSwak.dir/requires: deps/swak/test/CMakeFiles/TestSwak.dir/TestSwak.cpp.o.requires
deps/swak/test/CMakeFiles/TestSwak.dir/requires: deps/swak/test/CMakeFiles/TestSwak.dir/TestHelpers.cpp.o.requires
deps/swak/test/CMakeFiles/TestSwak.dir/requires: deps/swak/test/CMakeFiles/TestSwak.dir/TestSystem.cpp.o.requires
deps/swak/test/CMakeFiles/TestSwak.dir/requires: deps/swak/test/CMakeFiles/TestSwak.dir/TestSafeVec.cpp.o.requires
deps/swak/test/CMakeFiles/TestSwak.dir/requires: deps/swak/test/CMakeFiles/TestSwak.dir/TestExtractDigits.cpp.o.requires
deps/swak/test/CMakeFiles/TestSwak.dir/requires: deps/swak/test/CMakeFiles/TestSwak.dir/TestBinaryIO.cpp.o.requires
.PHONY : deps/swak/test/CMakeFiles/TestSwak.dir/requires

deps/swak/test/CMakeFiles/TestSwak.dir/clean:
	cd /media/fusion10/work/chromatinVariation/src/ase_cpp/build/debug/deps/swak/test && $(CMAKE_COMMAND) -P CMakeFiles/TestSwak.dir/cmake_clean.cmake
.PHONY : deps/swak/test/CMakeFiles/TestSwak.dir/clean

deps/swak/test/CMakeFiles/TestSwak.dir/depend:
	cd /media/fusion10/work/chromatinVariation/src/ase_cpp/build/debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /media/fusion10/work/chromatinVariation/src/ase_cpp /media/fusion10/work/chromatinVariation/src/ase_cpp/deps/swak/test /media/fusion10/work/chromatinVariation/src/ase_cpp/build/debug /media/fusion10/work/chromatinVariation/src/ase_cpp/build/debug/deps/swak/test /media/fusion10/work/chromatinVariation/src/ase_cpp/build/debug/deps/swak/test/CMakeFiles/TestSwak.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : deps/swak/test/CMakeFiles/TestSwak.dir/depend

