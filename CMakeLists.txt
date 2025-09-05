cmake_minimum_required(VERSION 3.10)
project(MyGraphProject)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Для Windows нужно линковать нужные библиотеки
if(WIN32)
    set(CMAKE_EXE_LINKER_FLAGS "-mwindows")
endif()

add_executable(MyGraphProject
    main.cpp
    MyGraph.cpp
    MyGraph.h
)

# Для Windows нужно линковать Gdi32 и User32
if(WIN32)
    target_link_libraries(MyGraphProject Gdi32 User32)
endif()