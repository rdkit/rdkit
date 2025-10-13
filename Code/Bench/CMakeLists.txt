add_executable(bench EXCLUDE_FROM_ALL smiles.cpp)
target_link_libraries(bench rdkitCatch GraphMol SmilesParse)
