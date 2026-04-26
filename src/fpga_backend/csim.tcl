# ===================================================================
# csim.tcl - C Simulation Script (由 Makefile 动态驱动)
# ===================================================================

# 1. 从环境变量获取 Makefile 传来的参数
set top_func  $::env(HLS_TOP_FUNC)
set src_files $::env(HLS_SRC_FILES)
set tb_files  $::env(HLS_TB_FILES)
set extra_cf  [expr {[info exists ::env(HLS_EXTRA_CFLAGS)] ? $::env(HLS_EXTRA_CFLAGS) : ""}]

puts "======================================="
puts "INFO: Running C Simulation for: $top_func"
puts "INFO: Source files : $src_files"
puts "INFO: Testbench    : $tb_files"
puts "INFO: Extra cflags : $extra_cf"
puts "======================================="

# 2. 每个模块使用独立工程目录，避免互相覆盖
file mkdir Solution
set proj_name "Solution/${top_func}"
open_project $proj_name

set_top $top_func
# 使用绝对路径，避免 HLS 把相对路径按工程子目录解析导致找不到源文件
set inc_dir [file normalize ./include]
set my_cflags "-I${inc_dir} -I/opt/xilinx/xrt/include $extra_cf"

# 3. 动态添加 Source 文件
foreach file $src_files {
    add_files [file normalize $file] -cflags $my_cflags
}

# 4. 动态添加 Testbench 文件
if {[string length $tb_files] > 0} {
    foreach file $tb_files {
        add_files [file normalize $file] -cflags $my_cflags -tb
    }
} else {
    puts "WARNING: No testbench files provided for $top_func! csim_design will fail."
}

open_solution "solution1"
set_part xcu55c-fsvh2892-2L-e
create_clock -period 6ns

# 清理旧的编译缓存，防止遗留 obj 文件引发 multiple definition 错误
file delete -force ${proj_name}/solution1/csim

csim_design

exit
