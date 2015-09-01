add_definitions( -DFREEGLUT_EXPORTS -DFREEGLUT_SHARED -D_CRT_SECURE_NO_WARNINGS )

file (GLOB_RECURSE SOURCES "*.cpp")
file (GLOB_RECURSE MAIN "*main.cpp")
file (GLOB_RECURSE SIM "*caging.cpp")
list (REMOVE_ITEM SOURCES ${MAIN} ${SIM})

# create library
add_library (${CMAKE_PROJECT_NAME}_Core ${CAGING_LIB_TYPE} ${SOURCES})
target_link_libraries (${CMAKE_PROJECT_NAME}_Core ${DEPENDENCY_LIBS})

# create executable
add_executable(alpha main.cpp)
target_link_libraries(alpha ${CMAKE_PROJECT_NAME}_Core ${DEPENDENCY_LIBS})

# create executable
add_executable(compound compound_grippers_main.cpp)
target_link_libraries(compound ${CMAKE_PROJECT_NAME}_Core ${DEPENDENCY_LIBS})

