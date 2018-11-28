include_directories(${CMAKE_SOURCE_DIR}/inc)

add_library (DLUtils SHARED DLUtil.cxx)

add_executable (DLSens DLSens.cxx)
add_executable (DLCovar DLCovariance.cxx)
add_executable (DLUnfold DLUnfold.cxx)
add_executable (DLCombinedFit DLCombinedFit.cxx)
add_executable (DLCombinedFitFake DLCombinedFitFake.cxx)
add_executable (DLCombinedFitDOF DLCombinedFitDOF.cxx)

target_link_libraries (DLSens SBNfitlib tinyxmllib EventWeight)
target_link_libraries (DLCovar SBNfitlib tinyxmllib EventWeight)
target_link_libraries (DLUnfold SBNfitlib tinyxmllib EventWeight)
target_link_libraries (DLCombinedFit SBNfitlib tinyxmllib EventWeight DLUtils)
target_link_libraries (DLCombinedFitFake SBNfitlib tinyxmllib EventWeight DLUtils)
target_link_libraries (DLCombinedFitDOF SBNfitlib tinyxmllib EventWeight DLUtils)

