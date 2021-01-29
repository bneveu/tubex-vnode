cmake_minimum_required(VERSION 3.0)

project(tubex-vnode)



################################################################################
# Looking for Ibex, Tubex
################################################################################

  set(IBEX_ROOT "/home/neveub/tubex/ibex-lib/install") # Install path of ibex lib
  set(TUBEX_ROOT "/home/neveub/tubex/tubex-lib/")
  set(TUBEX_SOLVER_ROOT "/home/neveub/tubex/tubex-solve/")
 
  set(IBEX_LIB_DIR "${IBEX_ROOT}/lib"
                   "${IBEX_ROOT}/lib/ibex/3rd")
  
  set(IBEX_INCLUDE_DIR  "${IBEX_ROOT}/include"
		        "${IBEX_ROOT}/include/ibex"
		        "${IBEX_ROOT}/include/ibex/3rd"
			"${IBEX_ROOT}/include/lib" 
			"${IBEX_ROOT}/include/lib/ibex/3rd")

  set(TUBEX_INCLUDE_DIR "${TUBEX_ROOT}/build/include"
			"${TUBEX_ROOT}/build/src/core"
			"${TUBEX_ROOT}/include/tubex"
			"${TUBEX_ROOT}/include/lib")
  set(TUBEX_SOLVER_INCLUDE_DIR "${TUBEX_SOLVER_ROOT}/src" 
    "${TUBEX_SOLVER_ROOT}/build/src")

  set(TUBEX_LIB_DIR "${TUBEX_ROOT}/build/src/core")
    #"${TUBEX_ROOT}/build/src/robotics")
  set(TUBEX_SOLVER_LIB_DIR "${TUBEX_SOLVER_ROOT}/build/src")

    


  include_directories(${IBEX_INCLUDE_DIR})
  include_directories(${TUBEX_INCLUDE_DIR})
  include_directories(${TUBEX_SOLVER_INCLUDE_DIR})

  link_directories(${IBEX_LIB_DIR})
  link_directories(${TUBEX_LIB_DIR})
  link_directories(${TUBEX_SOLVER_LIB_DIR})
  
  # set(TUBEX_LDFLAGS "tubex" "ibex" "gaol" "gdtoa" "ultim" "prim")
  set(TUBEX_LDFLAGS "tubex-solver" "tubex" "ibex" "prim" )

###########################################################"
# VNODE VARIABLES
set(VNODE_BUILD_DIR "/home/neveub/tubex/vnode/")
set(VNODE_LIB_DIR "${VNODE_BUILD_DIR}/vnodelp/lib")

# VNODE INCLUDE/LINK DIRECTORIES
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -frounding-math")
add_definitions(-DFILIB_VNODE)
add_definitions(-DMAXORDER=100)
include_directories(${VNODE_BUILD_DIR}/vnodelp/include)
include_directories(${VNODE_BUILD_DIR}/vnodelp/FADBAD++)
link_directories(${VNODE_LIB_DIR})

set(VNODE_LIBS vnode lapack blas tubex-solve tubex  ibex prim)

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
#add_executable(b17 examples/BVP_17.cpp )
#add_executable(b18 examples/BVP_18.cpp )
#add_executable(b19 examples/BVP_19.cpp )
#add_executable(b20 examples/BVP_20.cpp )

add_executable(bID examples/BVP_ID.cpp )
#add_executable(b23 examples/BVP_23.cpp )
#add_executable(b24 examples/BVP_24.cpp )
#add_executable(b25 examples/BVP_25.cpp )
#add_executable(b26 examples/BVP_26.cpp )
#add_executable(b32 examples/BVP_32.cpp )
#add_executable(b35 examples/BVP_35.cpp )
#add_executable(b40 examples/BVP_40.cpp )

add_executable(bvp01a examples/bvpsolve/bvp01a.cpp)
add_executable(bvp01b examples/bvpsolve/bvp01b.cpp)
add_executable(bvp01c examples/bvpsolve/bvp01c.cpp)
add_executable(bvp01d examples/bvpsolve/bvp01d.cpp)

add_executable(bvp02a examples/bvpsolve/bvp02a.cpp)
add_executable(bvp02b examples/bvpsolve/bvp02b.cpp)
add_executable(bvp02c examples/bvpsolve/bvp02c.cpp)
add_executable(bvp02d examples/bvpsolve/bvp02d.cpp)

add_executable(bvp22a examples/bvpsolve/bvp22a.cpp)
add_executable(bvp22b examples/bvpsolve/bvp22b.cpp)
add_executable(bvp22c examples/bvpsolve/bvp22c.cpp)

add_executable(bvp25a examples/bvpsolve/bvp25a.cpp)
add_executable(bvp25b examples/bvpsolve/bvp25b.cpp)
add_executable(bvp25c examples/bvpsolve/bvp25c.cpp)

add_executable(bvp32a examples/bvpsolve/bvp32a.cpp)
add_executable(bvp32b examples/bvpsolve/bvp32b.cpp)
add_executable(bvp32c examples/bvpsolve/bvp32c.cpp)
add_executable(bvp32d examples/bvpsolve/bvp32d.cpp)
add_executable(bvp32e examples/bvpsolve/bvp32e.cpp)


add_executable(bvp34a examples/bvpsolve/bvp34a.cpp)
add_executable(bvp34b examples/bvpsolve/bvp34b.cpp)
add_executable(bvp34c examples/bvpsolve/bvp34c.cpp)
add_executable(bvp34d examples/bvpsolve/bvp34d.cpp)

add_executable(bvp35a examples/bvpsolve/bvp35a.cpp)
add_executable(bvp35b examples/bvpsolve/bvp35b.cpp)
add_executable(bvp35c examples/bvpsolve/bvp35c.cpp)

add_executable(bvp40 examples/bvpsolve/bvp40.cpp)

target_link_libraries(b04  ${VNODE_LIBS} CtcVnodelp)
target_link_libraries(b09  ${VNODE_LIBS} CtcVnodelp)
target_link_libraries(b15  ${VNODE_LIBS} CtcVnodelp)
target_link_libraries(b16  ${VNODE_LIBS} CtcVnodelp)
#target_link_libraries(b17  ${VNODE_LIBS} CtcVnodelp)
#target_link_libraries(b18  ${VNODE_LIBS} CtcVnodelp)
#target_link_libraries(b19  ${VNODE_LIBS} CtcVnodelp)
#target_link_libraries(b20  ${VNODE_LIBS} CtcVnodelp)

target_link_libraries(bID  ${VNODE_LIBS} CtcVnodelp)
#target_link_libraries(b23  ${VNODE_LIBS} CtcVnodelp)
#target_link_libraries(b24  ${VNODE_LIBS} CtcVnodelp)
#target_link_libraries(b25  ${VNODE_LIBS} CtcVnodelp)
#target_link_libraries(b26  ${VNODE_LIBS} CtcVnodelp)
#target_link_libraries(b32  ${VNODE_LIBS} CtcVnodelp)
#target_link_libraries(b35  ${VNODE_LIBS} CtcVnodelp)
#target_link_libraries(b40  ${VNODE_LIBS} CtcVnodelp)

target_link_libraries(bvp01a ${VNODE_LIBS} CtcVnodelp)
target_link_libraries(bvp01b ${VNODE_LIBS} CtcVnodelp)
target_link_libraries(bvp01c ${VNODE_LIBS} CtcVnodelp)
target_link_libraries(bvp01d ${VNODE_LIBS} CtcVnodelp)

target_link_libraries(bvp02a ${VNODE_LIBS} CtcVnodelp)
target_link_libraries(bvp02b ${VNODE_LIBS} CtcVnodelp)
target_link_libraries(bvp02c ${VNODE_LIBS} CtcVnodelp)
target_link_libraries(bvp02d ${VNODE_LIBS} CtcVnodelp)

target_link_libraries(bvp22a ${VNODE_LIBS} CtcVnodelp)
target_link_libraries(bvp22b ${VNODE_LIBS} CtcVnodelp)
target_link_libraries(bvp22c ${VNODE_LIBS} CtcVnodelp)

target_link_libraries(bvp25a ${VNODE_LIBS} CtcVnodelp)
target_link_libraries(bvp25b ${VNODE_LIBS} CtcVnodelp)
target_link_libraries(bvp25c ${VNODE_LIBS} CtcVnodelp)

target_link_libraries(bvp32a ${VNODE_LIBS} CtcVnodelp)
target_link_libraries(bvp32b ${VNODE_LIBS} CtcVnodelp)
target_link_libraries(bvp32c ${VNODE_LIBS} CtcVnodelp)
target_link_libraries(bvp32d ${VNODE_LIBS} CtcVnodelp)
target_link_libraries(bvp32e ${VNODE_LIBS} CtcVnodelp)

target_link_libraries(bvp34a ${VNODE_LIBS} CtcVnodelp)
target_link_libraries(bvp34b ${VNODE_LIBS} CtcVnodelp)
target_link_libraries(bvp34c ${VNODE_LIBS} CtcVnodelp)
target_link_libraries(bvp34d ${VNODE_LIBS} CtcVnodelp)

target_link_libraries(bvp35a ${VNODE_LIBS} CtcVnodelp)
target_link_libraries(bvp35b ${VNODE_LIBS} CtcVnodelp)
target_link_libraries(bvp35c ${VNODE_LIBS} CtcVnodelp)

target_link_libraries(bvp40 ${VNODE_LIBS} CtcVnodelp)



#bubibex

#add_executable(pendule bubibex/pendule.cpp vnode_tubex/CtcVnodelp.cpp vnode_tubex/CtcVnodelp.h)
#target_link_libraries(pendule  ${VNODE_LIBS})

#add_executable(dub bubibex/dubbins.cpp vnode_tubex/CtcVnodelp.cpp vnode_tubex/CtcVnodelp.h)
#target_link_libraries(dub  ${VNODE_LIBS})
