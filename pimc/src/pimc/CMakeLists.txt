cmake_minimum_required(VERSION 3.19.1)

project(pigsl.e)

# Set c++ standard
set(CMAKE_CXX_STANDARD 17)

# Check for static build
if(STATIC)
    message("Static build specified, setting library suffixes to ${CMAKE_STATIC_LIBRARY_SUFFIX}.")
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX})
    set(BUILD_SHARED_LIBS OFF)
    list(APPEND CMAKE_EXE_LINKER_FLAGS "-static")
    list(APPEND CMAKE_EXE_LINKER_FLAGS_DEBUG "-static")
    set(Boost_USE_STATIC_LIBS   ON)
endif()

# Determine the dimension of space (default = 3)
if (NOT NDIM)
    set(NDIM 3)
endif()
add_definitions(-D NDIM=${NDIM})

# Set default build flags
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -Wall -fno-math-errno -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fno-math-errno -O3")

# Check for <span> flags
if (USE_TCB_SPAN)
    add_definitions(-D USE_TCB_SPAN=${USE_TCB_SPAN})
endif()

# Determine executable name
set(exe pigsl.e)
if(CMAKE_BUILD_TYPE MATCHES Debug)
    set(exe pigsld.e)
elseif(CMAKE_BUILD_TYPE MATCHES PIMC)
    set(exe pimcl.e)
elseif(CMAKE_BUILD_TYPE MATCHES PIMCDebug)
    set(exe pimcld.e)
endif()

# Find source files ( better to list explicitly https://stackoverflow.com/questions/1027247/specify-source-files-globally-with-glob )
file(GLOB PIMC_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
    )
add_executable(${exe} ${PIMC_SRC})

# Define headers for target
target_include_directories(${exe} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    PRIVATE src)

# Set debug flags
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D PIMC_DEBUG")

# Add PIMC build mode
set(CMAKE_CXX_FLAGS_PIMC "${CMAKE_CXX_FLAGS} -DPIMC")
set(CMAKE_CXX_FLAGS_PIMC "${CMAKE_CXX_FLAGS_PIMC}" CACHE STRING
    "Flags used by the C++ compiler during PIMC builds."
    FORCE )
set(CMAKE_C_FLAGS_PIMC "" CACHE STRING
    "Flags used by the C compiler during PIMC builds."
    FORCE )
set(CMAKE_EXE_LINKER_FLAGS_PIMC
    "${CMAKE_EXE_LINKER_FLAGS}" CACHE STRING
    "Flags used for linking binaries during PIMC builds."
    FORCE )
set(CMAKE_SHARED_LINKER_FLAGS_PIMC
    "${CMAKE_SHARED_LINKER_FLAGS}" CACHE STRING
    "Flags used by the shared libraries linker during PIMC builds."
    FORCE )

# Add PIMCDEBUG build mode
set(CMAKE_CXX_FLAGS_PIMCDEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DPIMC")
set(CMAKE_CXX_FLAGS_PIMCDEBUG "${CMAKE_CXX_FLAGS_PIMCDEBUG}" CACHE STRING
    "Flags used by the C++ compiler during PIMC debug builds."
    FORCE )
set(CMAKE_C_FLAGS_PIMCDEBUG "${CMAKE_C_FLAGS_DEBUG}" CACHE STRING
    "Flags used by the C compiler during PIMC debug builds."
    FORCE )
set(CMAKE_EXE_LINKER_FLAGS_PIMCDEBUG
    "${CMAKE_EXE_LINKER_FLAGS_DEBUG}" CACHE STRING
    "Flags used for linking binaries during PIMC debug builds."
    FORCE )
set(CMAKE_SHARED_LINKER_FLAGS_PIMCDEBUG
    "${CMAKE_SHARED_LINKER_FLAGS_DEBUG}" CACHE STRING
    "Flags used by the shared libraries linker during PIMC debug builds."
    FORCE )

mark_as_advanced(
    CMAKE_CXX_FLAGS_PIMC
    CMAKE_C_FLAGS_PIMC
    CMAKE_EXE_LINKER_FLAGS_PIMC
    CMAKE_SHARED_LINKER_FLAGS_PIMC
    CMAKE_CXX_FLAGS_PIMCDEBUG
    CMAKE_C_FLAGS_PIMCDEBUG
    CMAKE_EXE_LINKER_FLAGS_PIMCDEBUG
    CMAKE_SHARED_LINKER_FLAGS_PIMCDEBUG )

if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE None
       CACHE STRING "Choose the type of build : None Debug Release PIMC PIMCDebug."
       FORCE)
endif()

# Find Boost
find_package( Boost 1.73.0 REQUIRED)

# Add include directories
include_directories( ${Boost_INCLUDE_DIRS} )

# Link libraries
target_link_libraries (${exe} ${Boost_LIBRARIES} )

# 'make install' to the correct locations (provided by GNUInstallDirs).
install(TARGETS ${exe} EXPORT ${exe}Config
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR})  # This is for Windows
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# This makes the project importable from the install directory
# Put config file in per-project dir (name MUST match), can also
# just go into 'cmake'.
install(EXPORT ${exe}Config DESTINATION share/${exe}/cmake)

# This makes the project importable from the build directory
export(TARGETS ${exe} FILE ${exe}Config.cmake)

# Unit tests
# We need to make some of these
