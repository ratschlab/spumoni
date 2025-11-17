#-------------------------------------------------------------------
# Builds the executable for compressing the dictionary
#-------------------------------------------------------------------
add_executable(compress_dictionary compress_dictionary.cpp)
target_link_libraries(compress_dictionary common_h sdsl malloc_count)
target_include_directories(compress_dictionary PUBLIC  "../include")

FetchContent_GetProperties(r-index)
FetchContent_GetProperties(shaped_slp)
FetchContent_GetProperties(klib)
FetchContent_GetProperties(bigbwt)

set(FOLCA_SOURCE_DIR ${shaped_slp_SOURCE_DIR}/folca)
set(SUX_SOURCE_DIR ${shaped_slp_SOURCE_DIR}/external/sux/sux)

#-------------------------------------------------------------------
# Builds the final executable 
#-------------------------------------------------------------------

## Declare some variables
set(bonsai_SOURCE_DIR "${CMAKE_BINARY_DIR}/_deps/bonsai-src")
set(zlib_SOURCE_DIR "${CMAKE_BINARY_DIR}/_deps/zlib-src")
set(klib_SOURCE_DIR "${CMAKE_BINARY_DIR}/_deps/klib-src")

set(CMAKE_VERBOSE_MAKEFILE ON)

add_executable(spumoni spumoni.cpp  compute_ms_pml.cpp doc_array.cpp 
                        refbuilder.cpp emp_null_database.cpp 
                        ks_test.cpp batch_loader.cpp)
target_link_libraries(spumoni PUBLIC sdsl common_h divsufsort divsufsort64 ri pthread zlib bonsai "-fopenmp")
target_include_directories(spumoni PUBLIC
                            "../include"
                            "${shaped_slp_SOURCE_DIR}" "${FOLCA_SOURCE_DIR}" "${SUX_SOURCE_DIR}/function" 
                            "${SUX_SOURCE_DIR}/support" 
                            "${bonsai_SOURCE_DIR}/include/bonsai"
                            "${bonsai_SOURCE_DIR}"
                            "${bonsai_SOURCE_DIR}/hll/include/"
                            "${bonsai_SOURCE_DIR}/include")
target_compile_options(spumoni PUBLIC "-std=c++17" "-fopenmp")

pybind11_add_module(_core
        pybindings.cpp
        spumoni.cpp
        compute_ms_pml.cpp
        doc_array.cpp
        refbuilder.cpp
        emp_null_database.cpp
        ks_test.cpp
        batch_loader.cpp
)

target_link_libraries(_core PUBLIC sdsl common_h divsufsort divsufsort64 ri pthread zlib bonsai "-fopenmp" "-fPIC")
target_include_directories(_core PUBLIC
        "../include"
        "../extern/parlaylib/include"
        "${shaped_slp_SOURCE_DIR}" "${FOLCA_SOURCE_DIR}" "${SUX_SOURCE_DIR}/function"
        "${SUX_SOURCE_DIR}/support"
        "${bonsai_SOURCE_DIR}/include/bonsai"
        "${bonsai_SOURCE_DIR}"
        "${bonsai_SOURCE_DIR}/hll/include/"
        "${bonsai_SOURCE_DIR}/include")
target_compile_options(_core PUBLIC "-std=c++17" "-fopenmp" "-fPIC")

# ===============================================
# Copy the built _core module into the Python package
# ===============================================

# Destination where the Python package lives in your source tree
set(PY_PACKAGE_DIR "${CMAKE_SOURCE_DIR}/src/pyspumoni")

# The compiled shared library (platform-independent pattern)
add_custom_command(
        TARGET _core POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "$<TARGET_FILE:_core>"
        "${PY_PACKAGE_DIR}/$<TARGET_FILE_NAME:_core>"
        COMMENT "📦 Copying built _core library into ${PY_PACKAGE_DIR}"
)

# ===============================================
# Optional: Generate Python stub files (.pyi)
# ===============================================
# Only run if pybind11-stubgen is available

#find_program(PYBIND11_STUBGEN pybind11-stubgen)
#if(PYBIND11_STUBGEN)
#    add_custom_command(
#            TARGET _core POST_BUILD
#            COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_SOURCE_DIR}/src ${PYBIND11_STUBGEN} pyspumoni._core
#            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
#            COMMENT "🧩 Generating .pyi stubs using pybind11-stubgen"
#    )
#
#    # Optionally copy the generated stub to your package
#    add_custom_command(
#            TARGET _core POST_BUILD
#            COMMAND ${CMAKE_COMMAND} -E copy_if_different
#            "${CMAKE_SOURCE_DIR}/stubs/pyspumoni/_core.pyi"
#            "${PY_PACKAGE_DIR}/_core.pyi"
#            COMMENT "📄 Copying _core.pyi stub into ${PY_PACKAGE_DIR}"
#            VERBATIM
#    )
#else()
#    message(WARNING "⚠️ pybind11-stubgen not found — skipping stub generation")
#endif()

# === Target 2: The C++ Test Harness Executable ===
add_executable(
        test_harness  # The name of your C++ executable
        test_harness.cpp         # The source file
)
target_compile_features(test_harness PUBLIC cxx_std_17)

# Link the executable against pybind and python
# pybind11::embed handles all the necessary flags for embedding
target_link_libraries(
        test_harness PRIVATE
        pybind11::embed  # This links pybind11
        ${Python_LIBRARIES}  # This links libpython
)