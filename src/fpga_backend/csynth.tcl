# csynth.tcl
# 检查环境变量 TOP_FUNCTION，若未设置则使用默认值 "Top"
if { [info exists ::env(TOP_FUNCTION)] } {
    set top_func $::env(TOP_FUNCTION)
} else {
    set top_func "Top"   ;# 默认顶层
}
puts "INFO: Synthesizing top function: $top_func"

# 项目目录带上模块名后缀，方便区分
set proj_name "Solution_${top_func}"
open_project $proj_name

set_top $top_func

set my_cflags "-I./include -I/opt/xilinx/xrt/include"

add_files {
    ./src/top.cpp
    ./src/load.cpp
    ./src/arithmetic.cpp
    ./src/bconv.cpp
    ./src/auto.cpp
    ./src/mod_add_kernel.cpp
    ./src/mod_sub_kernel.cpp
    ./src/mod_mult_kernel.cpp
    ./src/ntt_kernel.cpp
    ./src/interleave.cpp
    ./src/bconv_naive.cpp
} -cflags $my_cflags

open_solution "solution1"
set_part xcu55c-fsvh2892-2L-e
create_clock -period 5ns

csynth_design