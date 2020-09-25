#cmake_minimum_required(VERSION 3.0)
cmake_minimum_required(VERSION 2.8)
project(tubex-vnode)



################################################################################
# Looking for Ibex, Tubex
################################################################################

  set(IBEX_ROOT "/certis/3/neveub/RECHERCHE/ibex/ibex2-8/ibex-lib/") # Install path of ibex lib
  set(TUBEX_ROOT "/certis/3/neveub/RECHERCHE/tubex/tubex2/tubex-lib/")
  set(TUBEX_SOLVER_ROOT "/certis/3/neveub/RECHERCHE/tubex/tubex-solve/tubex-solve/")
 
  set(FILIB_LIB "/certis/3/neveub/RECHERCHE/ibex/ibex2-8/ibex-lib/__build__/3rd/filibsrc-3.0.2.2/libprim/.libs/")
  
  set(IBEX_INCLUDE_DIR  "${IBEX_ROOT}/include"
		        "${IBEX_ROOT}/include/ibex"
		        "${IBEX_ROOT}/include/ibex/3rd"
			"${IBEX_ROOT}/include/lib" 
			"${IBEX_ROOT}/include/lib/ibex/3rd")

  set(TUBEX_INCLUDE_DIR "${TUBEX_ROOT}/make/include"
			"${TUBEX_ROOT}/make/src/core"
			"${TUBEX_ROOT}/include/tubex"
			"${TUBEX_ROOT}/include/tubex-rob"
			"${TUBEX_ROOT}/include/lib")
  set(TUBEX_SOLVER_INCLUDE_DIR "${TUBEX_SOLVER_ROOT}/include/tubex-solve" 
    "${TUBEX_SOLVER_ROOT}/make/src" "${TUBEX_SOLVER_ROOT}/src")

  set(TUBEX_LIB_DIR "${TUBEX_ROOT}/lib")
  set(TUBEX_SOLVER_LIB_DIR "${TUBEX_SOLVER_ROOT}/lib")


  include_directories(${IBEX_INCLUDE_DIR})
  include_directories(${TUBEX_INCLUDE_DIR})
  include_directories(${TUBEX_SOLVER_INCLUDE_DIR})
	
  link_directories("/certis/3/neveub/RECHERCHE/lapack/lapack-3.9.0/build/lib"
    "/certis/3/neveub/RECHERCHE/lapack/lapack-3.9.0/build/lib"
                    "${IBEX_ROOT}/lib" 
               	   "${IBEX_ROOT}/lib/ibex/3rd" 
                   "${TUBEX_ROOT}/make/src/core"
                   "${TUBEX_ROOT}/make/src/robotics"
                   "${TUBEX_SOLVER_ROOT}/make/src"
                   "${FILIB_LIB}")
# set(TUBEX_LDFLAGS "tubex" "ibex" "gaol" "gdtoa" "ultim" "prim")
set(TUBEX_LDFLAGS "tubex-solver" "tubex" "ibex" "prim" "gfortran")

###########################################################"
# VNODE VARIABLES
set(VNODE_BUILD_DIR "/certis/3/neveub/RECHERCHE/vnode/")
set(VNODE_LIB_DIR "${VNODE_BUILD_DIR}/vnodelp/lib")

# VNODE INCLUDE/LINK DIRECTORIES
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -frounding-math")
add_definitions(-DFILIB_VNODE)
add_definitions(-DMAXORDER=100)
include_directories(${VNODE_BUILD_DIR}/vnodelp/include)
include_directories(${VNODE_BUILD_DIR}/vnodelp/FADBAD++)
link_directories(${VNODE_LIB_DIR})
set(VNODE_LIBS vnode lapack blas gfortran prim tubex-solve tubex tubex-rob ibex prim)

include_directories(vnode_tubex)


add_library(CtcVnodelp STATIC vnode_tubex/CtcVnodelp.cpp vnode_tubex/CtcVnodelp.h)
target_link_libraries(CtcVnodelp  ${VNODE_LIBS})
 

#Set C++ version
set(CMAKE_CXX_STANDARD 11)
#ivp
add_executable(01 examples/IVP_01.cpp)
add_executable(02 examples/IVP_02.cpp)
add_executable(03 examples/IVP_03.cpp)
add_executable(10 examples/IVP_10.cpp)
add_executable(11 examples/IVP_11.cpp)
add_executable(12 examples/IVP_12.cpp)
add_executable(13 examples/IVP_13.cpp)
add_executable(14 examples/IVP_14.cpp)
add_executable(19 examples/IVP_19.cpp)
add_executable(20 examples/IVP_20.cpp)
add_executable(21 examples/IVP_21.cpp)



add_executable(test examples/test.cpp)
target_link_libraries(test  ${VNODE_LIBS} CtcVnodelp)



target_link_libraries(01  ${VNODE_LIBS} CtcVnodelp)
target_link_libraries(02  ${VNODE_LIBS} CtcVnodelp)
target_link_libraries(03  ${VNODE_LIBS} CtcVnodelp)
target_link_libraries(10  ${VNODE_LIBS} CtcVnodelp)
target_link_libraries(11  ${VNODE_LIBS} CtcVnodelp)
target_link_libraries(12  ${VNODE_LIBS} CtcVnodelp)
target_link_libraries(13  ${VNODE_LIBS} CtcVnodelp)
target_link_libraries(14  ${VNODE_LIBS} CtcVnodelp)
target_link_libraries(19  ${VNODE_LIBS} CtcVnodelp)
target_link_libraries(20  ${VNODE_LIBS} CtcVnodelp)
target_link_libraries(21  ${VNODE_LIBS} CtcVnodelp)



#bvp

add_executable(b04 examples/BVP_04.cpp )
add_executable(b09 examples/BVP_09.cpp )
add_executable(b15 examples/BVP_15.cpp )
add_executable(b16 examples/BVP_16.cpp )
add_executable(b17 examples/BVP_17.cpp )
add_executable(b18 examples/BVP_18.cpp )
add_executable(b19 examples/BVP_19.cpp )

add_executable(bID examples/BVP_ID.cpp )
add_executable(b23 examples/BVP_23.cpp )
add_executable(b24 examples/BVP_24.cpp )
add_executable(b25 examples/BVP_25.cpp )



target_link_libraries(b04  ${VNODE_LIBS} CtcVnodelp)
target_link_libraries(b09  ${VNODE_LIBS} CtcVnodelp)
target_link_libraries(b15  ${VNODE_LIBS} CtcVnodelp)
target_link_libraries(b16  ${VNODE_LIBS} CtcVnodelp)
target_link_libraries(b17  ${VNODE_LIBS} CtcVnodelp)
target_link_libraries(b18  ${VNODE_LIBS} CtcVnodelp)
target_link_libraries(b19  ${VNODE_LIBS} CtcVnodelp)

target_link_libraries(bID  ${VNODE_LIBS} CtcVnodelp)
target_link_libraries(b23  ${VNODE_LIBS} CtcVnodelp)
target_link_libraries(b24  ${VNODE_LIBS} CtcVnodelp)
target_link_libraries(b25  ${VNODE_LIBS} CtcVnodelp)

#PA Add
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -msse3")
add_compile_options(-std=c++11)
add_definitions(-std=c++11) # CMake 2.8.11 or older
message("Coucou: ${CMAKE_CXX_FLAGS}")
#PA Add end

#bubibex

#add_executable(pendule bubibex/pendule.cpp vnode_tubex/CtcVnodelp.cpp vnode_tubex/CtcVnodelp.h)
#target_link_libraries(pendule  ${VNODE_LIBS})

#add_executable(dub bubibex/dubbins.cpp vnode_tubex/CtcVnodelp.cpp vnode_tubex/CtcVnodelp.h)
#target_link_libraries(dub  ${VNODE_LIBS})