include_directories(swak ../deps/yaml-cpp/include)

add_library(swak
    Swak.cpp
    System.cpp
    Helpers.cpp
    PolyExe.cpp
    FieldReader.cpp
    StringUtil.cpp
    StlUtil.cpp
    FoldedOutStream.cpp
    FastaReader.cpp
    FastqReader.cpp
    OptionParser.cpp
    Bio.cpp
    SafeVec.cpp
    DensityEstimator.cpp
    Window.cpp
    ConfigReader.cpp
    gzstream.cpp
)

target_link_libraries(swak yaml-cpp z)

add_library(gzstream gzstream.cpp)
target_link_libraries(gzstream z)


# When we do make install... tell cmake where everything should go
install(TARGETS                                 swak
        RUNTIME             DESTINATION         bin
        LIBRARY             DESTINATION         lib
        ARCHIVE             DESTINATION         lib
)

file(GLOB header_files "swak/*.h")
install(FILES ${header_files} DESTINATION include/swak)
