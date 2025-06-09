!< MOTIOn test: basic use module test.
program motion_use_module_basic_test
!< MOTIOn test: basic use module test.
use motion, only : xdmf_file_object

implicit none
type(xdmf_file_object) :: xdmf_file      !< A XDMF file.
logical                :: test_passed(1) !< List of passed tests.

test_passed = .true. ! nothing to test, just run
print "(A,L1)", new_line('a')//'Are all tests passed? ', all(test_passed)
endprogram motion_use_module_basic_test

