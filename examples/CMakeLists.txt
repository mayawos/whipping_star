include_directories(${CMAKE_SOURCE_DIR}/inc)

add_executable (example1 example_1.cxx)
target_link_libraries (example1 SBNfitlib tinyxmllib EventWeight)

add_executable (example2 example_2.cxx)
target_link_libraries (example2 SBNfitlib tinyxmllib EventWeight)

add_executable (example3 example_3.cxx)
target_link_libraries (example3 SBNfitlib tinyxmllib EventWeight)