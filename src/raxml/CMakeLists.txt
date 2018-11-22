add_library(raxml
    avxLikelihood.c
    axml_extract.c
    evaluateGenericSpecial.c
    evaluatePartialGenericSpecial.c
    makenewzGenericSpecial.c
    mem_alloc.c
    newviewGenericSpecial.c
    optimizeModel.c
    searchAlgo.c
    models.c
)

if(UNIX)
# add avx stuff
endif(UNIX)

target_compile_definitions(raxml PRIVATE __SIM_SSE3 _GNU_SOURCE _FMA)
target_compile_options(raxml PRIVATE -march=native) # make better sse2 or whatever it actually needs

target_link_libraries(raxml m)