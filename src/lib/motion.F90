!< MOTIOn, Modular (HPC) Optimized Toolkit (for) IO (in fortra)n
module motion
!< MOTIOn, Modular (HPC) Optimized Toolkit (for) IO (in fortra)n

use motion_hdf5_file_object
use motion_xdmf_file_object
use motion_xh5f_file_object

implicit none
private
public :: hdf5_file_object, HDF5_PARAMETERS
public :: xdmf_file_object, XDMF_PARAMETERS
public :: xh5f_file_object, XH5F_PARAMETERS
endmodule motion
