# Read the routed checkpoint without reopening/mutating the implementation run.
# Report both export defaults and an explicit 0.75 ns setup uncertainty.
if {![info exists ::env(HKS_IMPL_PROJECT)]} { error "HKS_IMPL_PROJECT is required" }
set project [file normalize $::env(HKS_IMPL_PROJECT)]
set impl [file join $project solution1 impl verilog]
set checkpoints [glob -nocomplain [file join $impl project.runs impl_1 *_routed.dcp]]
if {[llength $checkpoints] != 1} { error "Expected one completed routed checkpoint" }
set reports [file join $impl report postroute_margin]
file mkdir $reports
set_param general.maxThreads 4
open_checkpoint [lindex $checkpoints 0]
report_route_status -file [file join $reports route_status.rpt]
report_utilization -file [file join $reports utilization.rpt]
report_timing_summary -file [file join $reports timing_default.rpt]
if {[llength [get_clocks ap_clk]] != 1} { error "Expected the ap_clk clock" }
set_clock_uncertainty -setup 0.75 [get_clocks ap_clk]
report_timing_summary -report_unconstrained -file [file join $reports timing_setup075.rpt]
report_timing -max_paths 10 -file [file join $reports paths_setup075.rpt]
report_high_fanout_nets -file [file join $reports high_fanout.rpt]
report_design_analysis -congestion -file [file join $reports congestion.rpt]
# Vivado does not provide the generic Tcl redirect command. Keep verbose
# diagnostics in the captured run log; each timing summary also has check_timing.
check_timing -verbose
# Frequency sensitivity only: same routed cells/nets, no RTL change or reroute.
# This is an explicit internal-path timing scenario, not a board clock setting.
create_clock -name ap_clk -period 7.0 [get_ports ap_clk]
set_clock_uncertainty -setup 0.75 [get_clocks ap_clk]
report_timing_summary -report_unconstrained -file [file join $reports timing_period7_setup075.rpt]
close_design
exit
