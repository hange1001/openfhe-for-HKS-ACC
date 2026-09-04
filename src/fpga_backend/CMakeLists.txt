if(OPENFHE_FPGA_ENABLE)
    message(FATAL_ERROR "Board-free oracle tests require OPENFHE_FPGA_ENABLE=OFF")
endif()
if(NOT BUILD_SHARED OR NOT NATIVE_SIZE EQUAL 64 OR WITH_REDUCED_NOISE)
    message(FATAL_ERROR "HKS C-model tests require BUILD_SHARED=ON, NATIVE_SIZE=64, WITH_REDUCED_NOISE=OFF")
endif()
set(HKS_HLS_ROOT "/tools/Xilinx/Vitis_HLS/2023.2" CACHE PATH "Vitis HLS installation (ap_int headers)")
if(NOT EXISTS "${HKS_HLS_ROOT}/include/ap_int.h")
    message(FATAL_ERROR "Set HKS_HLS_ROOT to a Vitis HLS installation containing include/ap_int.h")
endif()
set(HKS_KERNEL_SOURCES
    src/top.cpp src/load.cpp src/mod_add_kernel.cpp src/mod_mult_kernel.cpp
    src/mod_sub_kernel.cpp src/ntt_kernel.cpp src/arithmetic.cpp src/bconv.cpp
    src/interleave.cpp src/cg_ntt.cpp src/bconv_naive.cpp src/bconv_systolic.cpp)
add_library(hks_digit_cmodel STATIC ${HKS_KERNEL_SOURCES})
target_include_directories(hks_digit_cmodel PRIVATE include)
target_include_directories(hks_digit_cmodel SYSTEM PRIVATE "${HKS_HLS_ROOT}/include")
target_compile_options(hks_digit_cmodel PRIVATE -Wno-unknown-pragmas -Wno-unused-label -Wno-error=unused-but-set-variable)
target_link_libraries(hks_digit_cmodel PUBLIC gmp)
add_executable(hks-digit-openfhe-test testbench/hks_digit_openfhe_test.cpp)
target_include_directories(hks-digit-openfhe-test PRIVATE
    "${PROJECT_SOURCE_DIR}/src/core/include" "${PROJECT_BINARY_DIR}/src/core"
    "${PROJECT_SOURCE_DIR}/src/pke/include" "${PROJECT_SOURCE_DIR}/src/binfhe/include"
    "${PROJECT_SOURCE_DIR}/third-party/cereal/include")
target_link_libraries(hks-digit-openfhe-test PRIVATE hks_digit_cmodel OPENFHEpke OPENFHEcore OPENFHEbinfhe)
set_property(TARGET hks-digit-openfhe-test PROPERTY RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin/tests")
enable_testing()
add_test(NAME hks-digit-openfhe COMMAND hks-digit-openfhe-test)
set_tests_properties(hks-digit-openfhe PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2" TIMEOUT 180)
