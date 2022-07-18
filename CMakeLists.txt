###### Library ######
add_library(BayesianBlocks SHARED src/BayesianBlocks.cxx)
add_library(Fermitools::BayesianBlocks ALIAS BayesianBlocks)

target_include_directories(
  BayesianBlocks
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:>
)
target_link_libraries(BayesianBlocks PUBLIC astro st_facilities Likelihood irfInterface)

#################################################
### SWIG generates the python wrapper library ###
#################################################
include(UseSWIG)

set(UseSWIG_TARGET_NAME_PREFERENCE STANDARD)
set(CMAKE_SWIG_FLAGS -w362 -w389 -w401 -w402 -w503 -w509)

set_property(SOURCE src/BayesianBlocks.i PROPERTY CPLUSPLUS ON)
swig_add_library(py_BayesianBlocks LANGUAGE python SOURCES src/BayesianBlocks.i)
set_property(TARGET py_BayesianBlocks PROPERTY SWIG_USE_TARGET_INCLUDE_DIRECTORIES TRUE)
set_property(TARGET py_BayesianBlocks PROPERTY SWIG_MODULE_NAME _BayesianBlocks)

target_include_directories(
  py_BayesianBlocks
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:>
)

target_link_libraries(py_BayesianBlocks PUBLIC BayesianBlocks Python3::Python)

if(APPLE)
  target_link_options(py_BayesianBlocks PRIVATE -bundle -flat_namespace)
endif()


###############################################################
# Installation
###############################################################
install(DIRECTORY BayesianBlocks DESTINATION ${FERMI_INSTALL_INCLUDEDIR})
install(DIRECTORY python/ DESTINATION ${FERMI_INSTALL_PYTHON})
install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/BayesianBlocks.py
  DESTINATION ${FERMI_INSTALL_PYTHON}
)

install(
  TARGETS BayesianBlocks py_BayesianBlocks
  EXPORT fermiTargets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/BayesianBlocks.py
  DESTINATION ${FERMI_INSTALL_PYTHON}
)
