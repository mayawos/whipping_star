add_library(SBNfit SHARED SBNchi.cxx  SBNconfig.cxx  SBNcovariance.cxx prob.cxx SBNgenerate.cxx SBNfeld.cxx  SBNosc.cxx  SBNspec.cxx tinyxmlparser.cpp tinyxmlerror.cpp tinystr.cpp tinyxml.cpp)
install(TARGETS SBNfit
    DESTINATION ${CMAKE_INSTALL_PREFIX/lib}
    LIBRARY DESTINATION lib  # <-- Add this line
   COMPONENT library
    )

install(FILES ${MY_PUBLIC_HEADERS}
    DESTINATION ${CMAKE_INSTALL_PREFIX}/include
    COMPONENT devel)

add_executable              (highfivewrite highfivewrite.cpp)
add_executable              (nopreosc nopreosc.cpp)

target_link_libraries (highfivewrite LINK_PUBLIC ${libraries} ${GSL_LIBS} SBNfit Eigen3::Eigen)
install(TARGETS highfivewrite
    DESTINATION ${CMAKE_INSTALL_PREFIX/bin}
    RUNTIME DESTINATION bin
    )

target_link_libraries (nopreosc LINK_PUBLIC ${libraries} ${GSL_LIBS} SBNfit Eigen3::Eigen)
install(TARGETS nopreosc
    DESTINATION ${CMAKE_INSTALL_PREFIX/bin}
    RUNTIME DESTINATION bin
    )
