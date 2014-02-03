include_directories(. ../deps/swak/src ../src ../deps/swak/deps/yaml-cpp/include)

add_executable(TestRna
    TestRna.cpp
)
target_link_libraries(TestRna swak gene)

add_executable(Test
    Test.cpp
    TestRodUtil.cpp
)
target_link_libraries(Test swak rod)
