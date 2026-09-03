# Isolated board-free HKS prototype projects; never reset baseline solutions.
set mode $::env(HKS_MODE)
if {$mode ni {csim csynth cosim baseline}} { error "Unsupported HKS_MODE: $mode" }
set proj_name "Solution/hks_digit_${mode}"
open_project $proj_name
set_top Top
set flags "-std=c++14 -I[file normalize ./include] -DHKS_DIGIT_TB"
if {$mode eq "baseline"} {
    set flags "-std=c++14 -I[file normalize ./include] -DHKS_DISABLE_FUSED_OPCODE"
}
foreach source $::env(HLS_SRC_FILES) {
    add_files [file normalize $source] -cflags $flags
}
foreach tb {testbench/top_tb.cpp testbench/hks_digit_tb.cpp} {
    add_files [file normalize $tb] -tb -cflags $flags
}
open_solution solution1
set hls_part [expr {[info exists ::env(HLS_PART)] ? $::env(HLS_PART) : "xcu55c-fsvh2892-2L-e"}]
set_part $hls_part
create_clock -period 6ns
set_clock_uncertainty 0.75ns
config_interface -m_axi_max_widen_bitwidth 512
if {$mode eq "csim"} {
    csim_design
} else {
    csynth_design
    if {$mode eq "cosim"} { cosim_design -rtl verilog }
}
exit
