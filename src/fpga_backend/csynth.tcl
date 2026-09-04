# ===================================================================
# csynth.tcl - C Synthesis Script (由 Makefile 动态驱动)
# ===================================================================

# 1. 从环境变量获取 Makefile 传来的参数
set top_func  $::env(HLS_TOP_FUNC)
set src_files $::env(HLS_SRC_FILES)
set extra_cf  [expr {[info exists ::env(HLS_EXTRA_CFLAGS)] ? $::env(HLS_EXTRA_CFLAGS) : ""}]

puts "======================================="
puts "INFO: Running C Synthesis for: $top_func"
puts "INFO: Source files : $src_files"
puts "INFO: Extra cflags : $extra_cf"
puts "======================================="

# 2. 每个模块使用独立工程目录；实验可通过 HLS_PROJECT_NAME 覆盖，
#    避免覆盖基线 solution 或混入旧 RTL 产物。
file mkdir Solution
set proj_name [expr {[info exists ::env(HLS_PROJECT_NAME)] ? $::env(HLS_PROJECT_NAME) : "Solution/${top_func}"}]
open_project $proj_name

set_top $top_func
# 使用绝对路径，避免 HLS 把相对路径按工程子目录解析导致找不到源文件
set inc_dir [file normalize ./include]
set my_cflags "-I${inc_dir} -I/opt/xilinx/xrt/include $extra_cf"

# 3. 动态添加 Source 文件
foreach file $src_files {
    add_files [file normalize $file] -cflags $my_cflags
}

open_solution "solution1"
# 器件由 Makefile 的 HLS_PART 单点驱动，与 PLATFORM 保持同一块板
set hls_part [expr {[info exists ::env(HLS_PART)] ? $::env(HLS_PART) : "xcu55c-fsvh2892-2L-e"}]
set_part $hls_part
create_clock -period 6ns
set_clock_uncertainty 0.75ns

# 开启全局 AXI 位宽自动拓宽，上限 512 bit（Alveo U55C 原生宽度）
if {$top_func eq "Top"} {
    config_interface -m_axi_max_widen_bitwidth 256 -m_axi_alignment_byte_size 32
} else {
    config_interface -m_axi_max_widen_bitwidth 512
}

csynth_design

exit
