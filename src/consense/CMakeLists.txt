add_executable(consense
    cons.c
    consense.c
    phylip.c
)

target_link_libraries(consense PRIVATE m)