$(eval $(begin-module))

################################################################
# unit definitions
################################################################

# units
module_units_h += UNtoU3
module_units_cpp-h +=

# programs
module_programs_cpp +=

# test programs
module_programs_cpp_test += test_141 test_6114 test_input

################################################################
# library creation flag
################################################################

$(eval $(library))

################################################################
# special variable assignments, rules, and dependencies
################################################################

$(eval $(end-module))
