include_directories(${CMAKE_SOURCE_DIR}/inc)

#add_executable (MC_constraint_channels MC_constraint_channels.cxx)
#target_link_libraries (MC_constraint_channels SBNfitlib tinyxmllib EventWeight)

file(COPY ${CMAKE_SOURCE_DIR}/xml/example.xml DESTINATION ${CMAKE_SOURCE_DIR}/build/pelee_plotter)




