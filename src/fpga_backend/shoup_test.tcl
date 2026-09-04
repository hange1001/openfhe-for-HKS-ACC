# Independent primitive validation before integration into the BConv PE pool.
set mode [expr {[info exists ::env(SHOUP_MODE)] ? $::env(SHOUP_MODE) : "csim"}]
if {$mode ni {csim csynth cosim}} { error "Unsupported SHOUP_MODE" }
open_project Solution/shoup_primitive_$mode
set_top ShoupMul
add_files src/shoup.cpp -cflags "-std=c++14 -Iinclude"
add_files testbench/shoup_tb.cpp -tb -cflags "-std=c++14 -Iinclude"
open_solution solution1
set_part xcu55c-fsvh2892-2L-e
create_clock -period 6ns
set_clock_uncertainty 0.75ns
if {$mode eq "csim"} { csim_design } else {
    csynth_design
    if {$mode eq "cosim"} { cosim_design -rtl verilog }
}
exit
