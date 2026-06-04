# PatchBonsaiVec.cmake
# Called via cmake -P with VEC_FILE set to the path of vec.h.
# On ARM64, sse2neon.h is used to provide SSE2 intrinsics (see thirdparty/CMakeLists.txt).
# This script only patches the scalar fallback cases in vec.h as a safety net for
# non-x86 platforms that lack sse2neon support.
#
# IMPORTANT: Do NOT add architecture guards (#if defined(__x86_64__)) here —
# those would prevent x86intrin.h from being included, blocking sse2neon.

cmake_minimum_required(VERSION 3.15)

if(NOT VEC_FILE OR NOT EXISTS "${VEC_FILE}")
    message(FATAL_ERROR "PatchBonsaiVec.cmake: VEC_FILE='${VEC_FILE}' does not exist")
endif()

file(READ "${VEC_FILE}" _vec)

# Scalar fallback for uint32_t (already patched if "using Type = uint32_t" present)
if(NOT "${_vec}" MATCHES "using Type = uint32_t")
    string(REPLACE
        "declare_all_int128_32(epi32,)\n#else\n#error(\"Need at least sse2\")"
        "declare_all_int128_32(epi32,)\n#else\n    using Type = uint32_t;"
        _vec "${_vec}")
endif()

# Scalar fallback for uint16_t
if(NOT "${_vec}" MATCHES "using Type = uint16_t")
    string(REPLACE
        "declare_all_int128_16(epi16,)\n#else\n#error(\"Need at least sse2\")"
        "declare_all_int128_16(epi16,)\n#else\n    using Type = uint16_t;"
        _vec "${_vec}")
endif()

# Scalar fallback for uint8_t
if(NOT "${_vec}" MATCHES "using Type = uint8_t")
    string(REPLACE
        "_mm_slli_epi32(a, imm8));\n    }\n#else\n#error(\"Need at least sse2\")"
        "_mm_slli_epi32(a, imm8));\n    }\n#else\n    using Type = uint8_t;"
        _vec "${_vec}")
endif()

# Scalar fallback for float
if(NOT "${_vec}" MATCHES "using Type = float")
    string(REPLACE
        "#  error(\"Need at least sse2\")\n#endif\n    static constexpr size_t ALN = sizeof(Type) / sizeof(char);\n    static constexpr size_t MASK = ALN - 1;\n    static constexpr size_t COUNT = sizeof(Type) / sizeof(ValueType);\n    template<typename T>\n    static constexpr bool aligned(T *ptr) {\n        return (reinterpret_cast<uint64_t>(ptr) & MASK) == 0;\n    }\n    using VType = UType<SIMDTypes<ValueType>>;\n};\n\ntemplate<>\nstruct SIMDTypes<double>{"
        "    using Type = float;\n#endif\n    static constexpr size_t ALN = sizeof(Type) / sizeof(char);\n    static constexpr size_t MASK = ALN - 1;\n    static constexpr size_t COUNT = sizeof(Type) / sizeof(ValueType);\n    template<typename T>\n    static constexpr bool aligned(T *ptr) {\n        return (reinterpret_cast<uint64_t>(ptr) & MASK) == 0;\n    }\n    using VType = UType<SIMDTypes<ValueType>>;\n};\n\ntemplate<>\nstruct SIMDTypes<double>{"
        _vec "${_vec}")
endif()

# Scalar fallback for double
if(NOT "${_vec}" MATCHES "using Type = double")
    string(REPLACE
        "#  error(\"Need at least sse2\")"
        "    using Type = double;"
        _vec "${_vec}")
endif()

file(WRITE "${VEC_FILE}" "${_vec}")
message(STATUS "Applied bonsai vec.h scalar fallback patches at ${VEC_FILE}")
