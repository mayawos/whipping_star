add_library(SBNfit STATIC SBNchi.cxx  SBNconfig.cxx prob.cxx  SBNspec.cxx tinyxmlparser.cpp tinyxmlerror.cpp tinystr.cpp tinyxml.cpp)
install(TARGETS SBNfit
    DESTINATION ${CMAKE_INSTALL_PREFIX/lib}
    LIBRARY DESTINATION lib  # <-- Add this line
   COMPONENT library
    )

install(FILES ${MY_PUBLIC_HEADERS}
    DESTINATION ${CMAKE_INSTALL_PREFIX}/include
    COMPONENT devel)

add_executable              (pureeigen pureeigen.cpp)

target_link_libraries (pureeigen LINK_PUBLIC ${libraries} ${GSL_LIBS} SBNfit Eigen3::Eigen)
install(TARGETS pureeigen
    DESTINATION ${CMAKE_INSTALL_PREFIX/bin}
    RUNTIME DESTINATION bin
    )
