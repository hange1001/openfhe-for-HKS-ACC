# Board-free out-of-context implementation of an already synthesized project.
# Run from a separate log directory; this does NOT link the U55C platform shell.
if {![info exists ::env(HKS_IMPL_PROJECT)]} { error "HKS_IMPL_PROJECT is required" }
set project [file normalize $::env(HKS_IMPL_PROJECT)]
if {![file isdirectory [file join $project solution1 syn verilog]]} {
    error "Expected an existing synthesized solution1: $project"
}
open_project $project
open_solution solution1
export_design -rtl verilog -format ip_catalog -flow impl
exit
