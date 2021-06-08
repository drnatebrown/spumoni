FetchContent_GetProperties(r-index)
FetchContent_GetProperties(shaped_slp)
FetchContent_GetProperties(klib)
FetchContent_GetProperties(bigbwt)

set(FOLCA_SOURCE_DIR ${shaped_slp_SOURCE_DIR}/folca)
set(SUX_SOURCE_DIR ${shaped_slp_SOURCE_DIR}/external/sux/sux)

add_executable(rif_tests rif_tests.cpp)
target_link_libraries(rif_tests common sdsl divsufsort divsufsort64 malloc_count ri)
target_include_directories(rif_tests PUBLIC    "../../include/ms" 
                                        "../../include/common" 
                                        "${shaped_slp_SOURCE_DIR}" 
                                        "${FOLCA_SOURCE_DIR}" 
                                        "${SUX_SOURCE_DIR}/function" 
                                        "${SUX_SOURCE_DIR}/support"
                                        )
target_compile_options(rif_tests PUBLIC "-std=c++17")

add_executable(rle_string_tests rle_string_tests.cpp ${bigbwt_SOURCE_DIR}/xerrors.c)
target_link_libraries(rle_string_tests common sdsl divsufsort divsufsort64 malloc_count ri pthread)
target_include_directories(rle_string_tests PUBLIC    "../../include/ms" 
                                        "../../include/common" 
                                        "${bigbwt_SOURCE_DIR}"
                                        )
target_compile_options(rle_string_tests PUBLIC "-std=c++17")