cmake_minimum_required(VERSION 3.15)
project(nig CXX)

# ── C++ library ───────────────────────────────────────────────────────────────
set(NIG_SOURCES
    code/src/bessel_k0_scaled.cpp
    code/src/bessel_k1_scaled.cpp
    code/src/erfc.cpp
    code/src/nig_integration.cpp
    code/src/nig_series_beta_eq_zero.cpp
    code/src/nig_series_general.cpp
    code/src/nig_series_x_eq_mu.cpp
    code/src/nig.cpp
)

add_library(libnig STATIC ${NIG_SOURCES})
target_include_directories(libnig PUBLIC code/include)
target_compile_features(libnig PRIVATE cxx_std_17)
target_compile_options(libnig PRIVATE
    $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-O3 -march=native -mfma>
    $<$<CXX_COMPILER_ID:MSVC>:/O2>
)

# ── C++ executable ────────────────────────────────────────────────────────────
add_executable(nig_demo code/src/main.cpp)
target_link_libraries(nig_demo PRIVATE libnig)
target_compile_features(nig_demo PRIVATE cxx_std_17)

# ── Python extension (only when invoked via pip / scikit-build-core) ──────────
if(DEFINED SKBUILD)
    find_package(pybind11 CONFIG REQUIRED)
    pybind11_add_module(_nig code/src/nig_bindings.cpp)
    target_link_libraries(_nig PRIVATE libnig)
    install(TARGETS _nig DESTINATION nig)
endif()
