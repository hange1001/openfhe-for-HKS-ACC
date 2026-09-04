# Resume OOC physical implementation after a resource-pressure interruption.
# Reuses the completed synthesized project, never resets a project or run.
if {![info exists ::env(HKS_IMPL_PROJECT)]} { error "HKS_IMPL_PROJECT is required" }
set project [file normalize $::env(HKS_IMPL_PROJECT)]
set xpr [file join $project solution1 impl verilog project.xpr]
if {![file exists $xpr]} { error "No existing Vivado project: $xpr" }
set_param general.maxThreads 4
open_project $xpr
if {[get_property PROGRESS [get_runs synth_1]] ne "100%"} {
    error "Synthesis must already be complete; no automatic reset/resynthesis"
}
if {[get_property PROGRESS [get_runs impl_1]] ne "100%"} {
    launch_runs impl_1 -to_step route_design -jobs 1
    wait_on_run impl_1
}
if {[get_property PROGRESS [get_runs impl_1]] ne "100%"} {
    error "Implementation incomplete: [get_property STATUS [get_runs impl_1]]"
}
open_run impl_1
set reports [file join $project solution1 impl verilog report]
file mkdir $reports
report_route_status -file [file join $reports resumed_route_status.rpt]
report_utilization -file [file join $reports resumed_utilization.rpt]
report_timing_summary -file [file join $reports resumed_timing.rpt]
close_project
exit
