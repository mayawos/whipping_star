include_directories(${CMAKE_SOURCE_DIR}/inc)

add_executable (DLSens DLSens.cxx)
add_executable (DLCovar DLCovariance.cxx)
add_executable (UnfoldScale UnfoldScale.cxx)
add_executable (DLUnfold DLUnfold.cxx)
add_executable (DLSensCLS DLSensCLS.cxx)

target_link_libraries (DLSens SBNfitlib tinyxmllib EventWeight)
target_link_libraries (DLCovar SBNfitlib tinyxmllib EventWeight)
target_link_libraries (UnfoldScale SBNfitlib tinyxmllib EventWeight)
target_link_libraries (DLUnfold SBNfitlib tinyxmllib EventWeight)
target_link_libraries (DLSensCLS SBNfitlib tinyxmllib EventWeight)

