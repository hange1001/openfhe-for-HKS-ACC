<AutoPilot:project xmlns:AutoPilot="com.autoesl.autopilot.project" projectType="C/C++" name="Solution" ideType="classic" top="Compute_BConv">
    <Simulation argv="">
        <SimFlow name="csim" setup="false" optimizeCompile="false" clean="false" ldflags="" mflags=""/>
    </Simulation>
    <files>
        <file name="../../bconv_test.cpp" sc="0" tb="1" cflags="-I../../include -I/opt/xilinx/xrt/include -Wno-unknown-pragmas" csimflags="" blackbox="false"/>
        <file name="src/bconv.cpp" sc="0" tb="false" cflags="-I./include -I/opt/xilinx/xrt/include" csimflags="" blackbox="false"/>
        <file name="src/arithmetic.cpp" sc="0" tb="false" cflags="-I./include -I/opt/xilinx/xrt/include" csimflags="" blackbox="false"/>
        <file name="include/arithmetic.h" sc="0" tb="false" cflags="-I./include -I/opt/xilinx/xrt/include" csimflags="" blackbox="false"/>
        <file name="bconv.cpp" sc="0" tb="false" cflags="-I./include -I/opt/xilinx/xrt/include" csimflags="" blackbox="false"/>
    </files>
    <solutions>
        <solution name="solution1" status=""/>
    </solutions>
</AutoPilot:project>

