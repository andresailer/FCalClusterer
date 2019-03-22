

#FCal Utils Library
ADD_LIBRARY(FCalUtils SHARED src/RootUtils.cpp)
TARGET_INCLUDE_DIRECTORIES(FCalUtils PUBLIC include)
TARGET_INCLUDE_DIRECTORIES(FCalUtils SYSTEM PUBLIC ${ROOT_INCLUDE_DIR} ${ROOT_INCLUDE_DIRS})
TARGET_COMPILE_DEFINITIONS(FCalUtils PUBLIC ${ROOT_DEFINITIONS})
TARGET_LINK_LIBRARIES(FCalUtils ${ROOT_LIBRARIES})

#BCBackgroundPar Executable
ADD_EXECUTABLE( BCBackgroundPar src/BCBackgroundPar.cpp src/BackgroundFitter.cpp )
TARGET_LINK_LIBRARIES( BCBackgroundPar FCalUtils )


#Occupancies Executbale
IF( DD4hep_FOUND )
  ADD_EXECUTABLE(Occupancies src/Occupancies.cpp)
  ADD_EXECUTABLE(MergedOccupancies src/MergedOccupancies.cpp)
  ADD_EXECUTABLE(SafetyOccupancies src/SafetyOccupancies.cpp)
  TARGET_LINK_LIBRARIES(Occupancies BeamCalReco FCalUtils)
  TARGET_LINK_LIBRARIES(MergedOccupancies BeamCalReco FCalUtils)
  TARGET_LINK_LIBRARIES(SafetyOccupancies BeamCalReco FCalUtils)
  TARGET_INCLUDE_DIRECTORIES(Occupancies PUBLIC ${CMAKE_SOURCE_DIR}/source/BeamCalReco/include)
  TARGET_INCLUDE_DIRECTORIES(MergedOccupancies PUBLIC ${CMAKE_SOURCE_DIR}/source/BeamCalReco/include)
  TARGET_INCLUDE_DIRECTORIES(SafetyOccupancies PUBLIC ${CMAKE_SOURCE_DIR}/source/BeamCalReco/include)
  INSTALL(TARGETS Occupancies MergedOccupancies SafetyOccupancies DESTINATION bin)
ENDIF()

IF ( CPPCHECK_FOUND )
  ADD_CPPCHECK( FCalUtils UNUSED_FUNCTIONS STYLE POSSIBLE_ERRORS FORCE FAIL_ON_WARNINGS)
  ADD_CPPCHECK( BCBackgroundPar UNUSED_FUNCTIONS STYLE POSSIBLE_ERRORS FORCE FAIL_ON_WARNINGS)
ENDIF()

INSTALL( TARGETS
  FCalUtils
  LIBRARY DESTINATION lib)

INSTALL( TARGETS
  BCBackgroundPar
  DESTINATION bin)
