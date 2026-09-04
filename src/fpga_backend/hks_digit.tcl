# Isolated board-free HKS prototype projects; never reset baseline solutions.
set mode $::env(HKS_MODE)
if {$mode ni {csim csynth cosim cosim-smoke cosim-perf baseline}} { error "Unsupported HKS_MODE: $mode" }
set source_root [pwd]
if {[info exists ::env(HKS_SOURCE_ROOT)]} {
    set source_root [file normalize $::env(HKS_SOURCE_ROOT)]
}
set proj_name "Solution/hks_digit_${mode}"
if {[info exists ::env(HKS_PROJECT_TAG)]} {
    if {![regexp {^[A-Za-z0-9_-]+$} $::env(HKS_PROJECT_TAG)]} {
        error "HKS_PROJECT_TAG must contain only letters, digits, _ or -"
    }
    append proj_name "_$::env(HKS_PROJECT_TAG)"
}
open_project $proj_name
set_top Top
set flags "-std=c++14 -I[file join $source_root include] -DHKS_DIGIT_TB"
set testbenches {testbench/top_tb.cpp testbench/hks_digit_tb.cpp}
if {$mode eq "cosim-smoke"} {
    append flags " -DHKS_RTL_SMOKE"
    set testbenches {testbench/hks_digit_tb.cpp}
}
if {$mode eq "cosim-perf"} {
    if {![info exists ::env(HKS_RTL_FIXTURE)]} { error "HKS_RTL_FIXTURE is required" }
    set fixture [file normalize $::env(HKS_RTL_FIXTURE)]
    if {[file tail $fixture] ne "openfhe_rtl_fixture.txt"} { error "Unexpected fixture basename" }
    set testbenches {testbench/hks_digit_perf_tb.cpp}
    add_files $fixture -tb
}
if {$mode eq "baseline"} {
    set flags "-std=c++14 -I[file join $source_root include] -DHKS_DISABLE_FUSED_OPCODE"
}
foreach source $::env(HLS_SRC_FILES) {
    add_files [file normalize [file join $source_root $source]] -cflags $flags
}
foreach tb $testbenches {
    add_files [file normalize $tb] -tb -cflags $flags
}
open_solution solution1
set hls_part [expr {[info exists ::env(HLS_PART)] ? $::env(HLS_PART) : "xcu55c-fsvh2892-2L-e"}]
set_part $hls_part
create_clock -period 6ns
set_clock_uncertainty 0.75ns
config_interface -m_axi_max_widen_bitwidth 256 -m_axi_alignment_byte_size 32
if {$mode eq "csim"} {
    csim_design
} else {
    csynth_design
    if {$mode in {cosim cosim-smoke cosim-perf}} { cosim_design -rtl verilog }
}
exit
