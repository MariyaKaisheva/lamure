############################################################
# CMake Build Script for the pvs_visibility_analyzer executable

link_directories(${SCHISM_LIBRARY_DIRS})

include_directories(${PVS_COMMON_INCLUDE_DIR}
                    ${COMMON_INCLUDE_DIR}
                    ${GLUT_INCLUDE_DIR}
                    ${FREEIMAGE_INCLUDE_DIR}
			        ${LAMURE_CONFIG_DIR})

include_directories(SYSTEM ${SCHISM_INCLUDE_DIRS}
						   ${Boost_INCLUDE_DIR})

link_directories(${SCHISM_LIBRARY_DIRS})

InitApp(${CMAKE_PROJECT_NAME}_pvs_visibility_analyzer)

############################################################
# Libraries

target_link_libraries(${PROJECT_NAME}
    ${PROJECT_LIBS}
    ${PVS_COMMON_LIBRARY}
    optimized ${SCHISM_CORE_LIBRARY} debug ${SCHISM_CORE_LIBRARY_DEBUG}
    )

add_dependencies(${PROJECT_NAME} lamure_pvs_common lamure_common)

MsvcPostBuild(${PROJECT_NAME})
