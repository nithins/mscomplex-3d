project(pymscomplex3d)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${mscomplex3d_SOURCE_DIR}/cmake)

find_package(Boost COMPONENTS python REQUIRED)

find_package(PythonLibs REQUIRED)

set(${PROJECT_NAME}_SRCS ${${PROJECT_NAME}_SRCS}

  ${mscomplex3d_SOURCE_DIR}/dmsc/grid.h
  ${mscomplex3d_SOURCE_DIR}/dmsc/grid.cpp

  ${mscomplex3d_SOURCE_DIR}/dmsc/grid_dataset.h
  ${mscomplex3d_SOURCE_DIR}/dmsc/grid_dataset_ensure.h
  ${mscomplex3d_SOURCE_DIR}/dmsc/grid_dataset.cpp

  ${mscomplex3d_SOURCE_DIR}/dmsc/grid_mscomplex.h
  ${mscomplex3d_SOURCE_DIR}/dmsc/grid_mscomplex_ensure.h
  ${mscomplex3d_SOURCE_DIR}/dmsc/grid_mscomplex.cpp

  ${mscomplex3d_SOURCE_DIR}/dmsc/grid_datamanager.h
  ${mscomplex3d_SOURCE_DIR}/dmsc/grid_datamanager.cpp
  )

# in order to make this module independant, code from utls is compiled into it. 
# will change this to link to utls, if utls becomes installable

set(${PROJECT_NAME}_SRCS ${${PROJECT_NAME}_SRCS}

  ${mscomplex3d_SOURCE_DIR}/utls/include/cpputils.h
  ${mscomplex3d_SOURCE_DIR}/utls/src/cpputils.cpp

  ${mscomplex3d_SOURCE_DIR}/utls/include/n_vector.h
  ${mscomplex3d_SOURCE_DIR}/utls/src/n_vector.cpp
  ${mscomplex3d_SOURCE_DIR}/utls/include/aabb.h

  ${mscomplex3d_SOURCE_DIR}/utls/include/logutil.h
  ${mscomplex3d_SOURCE_DIR}/utls/src/logutil.cpp

  ${mscomplex3d_SOURCE_DIR}/utls/include/timer.h
  ${mscomplex3d_SOURCE_DIR}/utls/src/timer.cpp
  )

include_directories(
  ${mscomplex3d_SOURCE_DIR}/dmsc/
  ${mscomplex3d_SOURCE_DIR}/utls/include
  ${Boost_INCLUDE_DIRS}
  ${PYTHON_INCLUDE_DIRS}
  )

set(${PROJECT_NAME}_SRCS ${${PROJECT_NAME}_SRCS}
  ${CMAKE_CURRENT_SOURCE_DIR}/grid_python_interface.cpp
  )

set(${PROJECT_NAME}_LIBS ${${PROJECT_NAME}_LIBS} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})

add_library(${PROJECT_NAME} SHARED ${${PROJECT_NAME}_SRCS})

set_target_properties(${PROJECT_NAME} PROPERTIES PREFIX "")

target_link_libraries(${PROJECT_NAME} ${${PROJECT_NAME}_LIBS})
