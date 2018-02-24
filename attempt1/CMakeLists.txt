cmake_minimum_required(VERSION 3.9)
project(attempt1)

set(CMAKE_CXX_STANDARD 11)

add_executable(attempt1
        check.sh
        datagen.c
        Lab3IO.c
        Lab3IO.h
        main.c
        Makefile
        ReadMe
        serialtester.c
        timer.h)
