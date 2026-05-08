# amflow_cli — JSON driver.

add_executable(amflow_cli
    main.cpp
)

target_link_libraries(amflow_cli
    PRIVATE
        amflow
)

target_compile_features(amflow_cli PRIVATE cxx_std_17)
