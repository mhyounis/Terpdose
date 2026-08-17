#  ==================
#    CMake settings  
#  ==================

CMAKE_MINIMUM_REQUIRED (VERSION 3.20)
PROJECT (Terpdose LANGUAGES C Fortran)

#  ==========
#    Python  
#  ==========

FIND_PACKAGE (Python REQUIRED COMPONENTS Interpreter Development.Module NumPy)

#  =======================================
#    Fortran Directory and Compatibility  
#  =======================================

SET (FORTRAN_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/Terpdose/fortran")

EXECUTE_PROCESS (
  COMMAND "${Python_EXECUTABLE}" -c "import numpy.f2py; print(numpy.f2py.get_include())"
  OUTPUT_VARIABLE F2PY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

#  ================
#    Fortran Code  
#  ================

SET (TERPDOSE_FORTRAN "_tdfort")

# Directories

SET (CORE_DIR "${FORTRAN_DIR}/core")
SET (PP_DIR   "${FORTRAN_DIR}/postprocessing")
SET (GRID_DIR "${FORTRAN_DIR}/grid")
SET (FEA_DIR  "${FORTRAN_DIR}/fea")

SET (FORTRAN_SOURCES
    "${CORE_DIR}/constants.f90"
    "${CORE_DIR}/math.f90"
    
    "${PP_DIR}/CalculateFluence.f90"
    "${PP_DIR}/TallyDeposition.f90"
    
    "${GRID_DIR}/LeafElementMask.f90"
    "${GRID_DIR}/FindElement.f90"
    
    "${FEA_DIR}/ShapeFunctions.f90"
    "${FEA_DIR}/MapToReference.f90"
    "${FEA_DIR}/FEAInterpolation.f90"
)

SET (F2PY_C "${CMAKE_CURRENT_BINARY_DIR}/${TERPDOSE_FORTRAN}module.c")
SET (F2PY_F77_WRAPS "${CMAKE_CURRENT_BINARY_DIR}/${TERPDOSE_FORTRAN}-f2pywrappers.f")
SET (F2PY_F90_WRAPS "${CMAKE_CURRENT_BINARY_DIR}/${TERPDOSE_FORTRAN}-f2pywrappers2.f90")

ADD_CUSTOM_COMMAND (
  OUTPUT "${F2PY_C}" "${F2PY_F77_WRAPS}" "${F2PY_F90_WRAPS}"
  COMMAND "${Python_EXECUTABLE}" -m numpy.f2py
          ${FORTRAN_SOURCES}
          -m "${TERPDOSE_FORTRAN}"
          --lower
  DEPENDS ${FORTRAN_SOURCES}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
  VERBATIM
)

ADD_CUSTOM_TARGET (gen_f2py_wrappers DEPENDS "${F2PY_C}" "${F2PY_F77_WRAPS}" "${F2PY_F90_WRAPS}")

Python_add_library(${TERPDOSE_FORTRAN} MODULE WITH_SOABI
  "${F2PY_C}"
  "${F2PY_F77_WRAPS}"
  "${F2PY_F90_WRAPS}"
  "${F2PY_INCLUDE_DIR}/fortranobject.c"
  ${FORTRAN_SOURCES}
)

ADD_DEPENDENCIES (${TERPDOSE_FORTRAN} gen_f2py_wrappers)
TARGET_INCLUDE_DIRECTORIES (${TERPDOSE_FORTRAN} PRIVATE "${F2PY_INCLUDE_DIR}")
TARGET_LINK_LIBRARIES (${TERPDOSE_FORTRAN} PRIVATE Python::NumPy)

INSTALL (TARGETS ${TERPDOSE_FORTRAN} LIBRARY DESTINATION Terpdose)

INSTALL (
    DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/terpdose/"
    DESTINATION Terpdose
    FILES_MATCHING PATTERN "*.py"
)