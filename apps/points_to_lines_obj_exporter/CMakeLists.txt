############################################################
# CMake Build Script for the points_to_lines_obj_exporter executable

link_directories(${SCHISM_LIBRARY_DIRS})

include_directories(${REND_INCLUDE_DIR} 
                    ${COMMON_INCLUDE_DIR}
					"/home/vajo3185/Programming/lamure/external/"
                    )

include_directories(SYSTEM ${SCHISM_INCLUDE_DIRS}
						   ${Boost_INCLUDE_DIR})


InitApp(${CMAKE_PROJECT_NAME}_points_to_lines_obj_exporter)

############################################################
# Libraries

target_link_libraries(${PROJECT_NAME}
    ${PROJECT_LIBS}
    ${REND_LIBRARY}
    ${OpenGL_LIBRARIES} 
    ${GLUT_LIBRARY}
    /home/vajo3185/Programming/lamure/external/alglib/cpp/src/libAlglib.so
    )

