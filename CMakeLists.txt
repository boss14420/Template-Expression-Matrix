cmake_minimum_required (VERSION 2.6)
PROJECT (MATRIX)

INCLUDE_DIRECTORIES("${MATRIX}/src/classes")
INCLUDE_DIRECTORIES("${MATRIX}/src/trivial")

#FIND_LIBRARY(MATH m)
#FIND_LIBRARY(PTHREAD pthread)

#SET(CMAKE_CXX_FLAGS "-O3 -Wall -std=gnu++0x --param max-inline-recursive-depth-auto=32 --param max-inline-recursive-depth=32 -ftree-vectorize -msse2 -DNDEBUG_MT_PRODUCT_TIME -DNDEBUG")
#SET(CMAKE_CXX_FLAGS "-g -Wall -std=gnu++0x -DNDEBUG_MT_PRODUCT_TIME")
SET(CMAKE_CXX_FLAGS "-O2 -Wall -std=gnu++0x -fomit-frame-pointer -ftree-vectorize -march=native -DNDEBUG -DNDEBUG_MT_EXPRESSION_TIME")

ADD_EXECUTABLE(mattest
    src/test/mattest.cpp
  )

ADD_EXECUTABLE(matspeedtest
    src/test/matspeedtest.cpp
    src/trivial/matrix.hpp
  )

#ADD_EXECUTABLE(mattest2
#    src/test/mattest2.cpp
#  )

#TARGET_LINK_LIBRARIES(test_strassen_matrix ${MATH} ${PTHREAD})
#TARGET_LINK_LIBRARIES(test2 ${MATH} ${PTHREAD})
