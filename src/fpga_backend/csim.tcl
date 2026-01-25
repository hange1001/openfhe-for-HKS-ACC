open_project Solution

# 使用 Compute_BConv 作为顶层函数 (3D 数组接口)
set_top Compute_BConv

set my_cflags "-I./include -I/opt/xilinx/xrt/include"

# 添加源文件
add_files ./src/bconv.cpp -cflags $my_cflags
add_files ./src/arithmetic.cpp -cflags $my_cflags

# 添加测试文件
add_files ./bconv_test.cpp -cflags $my_cflags -tb

open_solution "solution1"
set_part xcu250-figd2104-2l-e
create_clock -period 4ns

# 运行 C 仿真
csim_design

exit