include_directories(. ../src ../deps/yaml-cpp/include)

add_executable(TestSwak TestSwak.cpp TestHelpers.cpp TestSystem.cpp TestSafeVec.cpp TestExtractDigits.cpp TestBinaryIO.cpp)
target_link_libraries(TestSwak swak)

add_executable(parse-yaml parse.cpp)
target_link_libraries(parse-yaml yaml-cpp)
