cmake_minimum_required(VERSION 3.9)
project(ECE420_Lab3)

set(CMAKE_CXX_STANDARD 11)

add_executable(ECE420_Lab3
        attempt1/datagen.c
        attempt1/Lab3IO.c
        attempt1/Lab3IO.h
        attempt1/main.c
        attempt1/serialtester.c
        attempt1/timer.h
        datagen.c
        Lab3IO.c
        Lab3IO.h
        main.c
        serialtester.c
        timer.h)
