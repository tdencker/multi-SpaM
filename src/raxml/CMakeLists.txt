add_library(raxml
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

include(CheckCXXSourceRuns)
set(CMAKE_REQUIRED_FLAGS "-mfma")
check_cxx_source_runs("
        #include <immintrin.h>
        int main()
        {
        __m256 a, b, c, d;
        const float src[8] = { 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f };
        float dst[8];
        a = _mm256_loadu_ps( src );
        b = _mm256_loadu_ps( src );
        c = _mm256_loadu_ps( src );
        d = _mm256_fmadd_ps( a, b, c);
        _mm256_storeu_ps( dst, d );
        for( int i = 0; i < 8; i++ ){
            if( ( ( src[i] * src[i] ) + src[i] ) != dst[i] ){
            return -1;
            }
        }
        return 0;
        }"
HAVE_FMA)
if(HAVE_FMA)
    target_sources(raxml PRIVATE avxLikelihood.c)
    target_compile_definitions(raxml PRIVATE __AVX _FMA)
    target_compile_options(raxml PRIVATE -mfma)
endif()

target_compile_definitions(raxml PRIVATE __SIM_SSE3 _GNU_SOURCE)
target_compile_options(raxml PRIVATE -msse3)

target_link_libraries(raxml m)