add_executable(main_a main_a.cpp)
target_link_libraries(main_a PRIVATE core)

add_executable(test_fairness test_fairness.cpp)
target_link_libraries(test_fairness PRIVATE core)

add_executable(tuning_ga tuning_ga.cpp)
target_link_libraries(tuning_ga PRIVATE core)

add_executable(tuning_ln tuning_ln.cpp)
target_link_libraries(tuning_ln PRIVATE core)

add_executable(test_algorithm test_algorithm.cpp)
target_link_libraries(test_algorithm PRIVATE core)

add_executable(competition competition.cpp)
target_link_libraries(competition PRIVATE core)
