!< MOTIOn, Modular (HPC) Optimized Toolkit (for) IO (in fortra)n
module motion
!< MOTIOn, Modular (HPC) Optimized Toolkit (for) IO (in fortra)n

use motion_hdf5_file_object
use motion_xdmf_file_object
use motion_xh5f_file_object

implicit none
private
public :: hdf5_file_object, HDF5_DATASPACE_TYPE_SIMPLE
public :: xdmf_file_object,                   &
          XDMF_ATTR_CENTER_NODE,              &
          XDMF_ATTR_CENTER_EDGE,              &
          XDMF_ATTR_CENTER_FACE,              &
          XDMF_ATTR_CENTER_CELL,              &
          XDMF_ATTR_CENTER_GRID,              &
          XDMF_ATTR_CENTER_OTHER,             &
          XDMF_ATTR_TYPE_SCALAR,              &
          XDMF_ATTR_TYPE_VECTOR,              &
          XDMF_ATTR_TYPE_TENSOR,              &
          XDMF_ATTR_TYPE_TENSOR6,             &
          XDMF_ATTR_TYPE_MATRIX,              &
          XDMF_ATTR_TYPE_GLOBALID,            &
          XDMF_ATTR_ELEM_FAMILY_CG,           &
          XDMF_ATTR_ELEM_FAMILY_DG,           &
          XDMF_ATTR_ELEM_FAMILY_Q,            &
          XDMF_ATTR_ELEM_FAMILY_DQ,           &
          XDMF_ATTR_ELEM_FAMILY_RT,           &
          XDMF_ATTR_ELEM_CELL_INTERVAL,       &
          XDMF_ATTR_ELEM_CELL_TRIANGLE,       &
          XDMF_ATTR_ELEM_CELL_TETRAHEDRON,    &
          XDMF_ATTR_ELEM_CELL_QUADRILATERAL,  &
          XDMF_ATTR_ELEM_CELL_HEXAHEDRON,     &
          XDMF_DATAITEM_ITEMTYPE_UNIFORM,     &
          XDMF_DATAITEM_ITEMTYPE_COLLECTION,  &
          XDMF_DATAITEM_ITEMTYPE_TREE,        &
          XDMF_DATAITEM_ITEMTYPE_HYPERSLAB,   &
          XDMF_DATAITEM_ITEMTYPE_COORDINATES, &
          XDMF_DATAITEM_ITEMTYPE_FUNCTION,    &
          XDMF_DATAITEM_NUMBER_FORMAT_HDF,    &
          XDMF_DATAITEM_NUMBER_FORMAT_XML,    &
          XDMF_DATAITEM_NUMBER_FORMAT_BINARY, &
          XDMF_DATAITEM_NUMBER_TYPE_FLOAT,    &
          XDMF_DATAITEM_NUMBER_TYPE_INT,      &
          XDMF_DATAITEM_NUMBER_TYPE_UINT,     &
          XDMF_DATAITEM_NUMBER_TYPE_CHAR,     &
          XDMF_DATAITEM_NUMBER_TYPE_UCHAR,    &
          XDMF_DATAITEM_NUMBER_PRECISION_1,   &
          XDMF_DATAITEM_NUMBER_PRECISION_2,   &
          XDMF_DATAITEM_NUMBER_PRECISION_4,   &
          XDMF_DATAITEM_NUMBER_PRECISION_8,   &
          XDMF_DATAITEM_ENDIAN_NATIVE,        &
          XDMF_DATAITEM_ENDIAN_BIG,           &
          XDMF_DATAITEM_ENDIAN_LITTLE,        &
          XDMF_DATAITEM_COMPRESSION_RAW,      &
          XDMF_DATAITEM_COMPRESSION_ZLIB,     &
          XDMF_DATAITEM_COMPRESSION_BZIP2,    &
          XDMF_GEOMETRY_TYPE_XYZ,             &
          XDMF_GEOMETRY_TYPE_XY,              &
          XDMF_GEOMETRY_TYPE_X_Y_Z,           &
          XDMF_GEOMETRY_TYPE_VXVYVZ,          &
          XDMF_GEOMETRY_TYPE_ODXYZ,           &
          XDMF_GEOMETRY_TYPE_ODXY,            &
          XDMF_GRID_TYPE_UNIFORM,             &
          XDMF_GRID_TYPE_COLLECTION,          &
          XDMF_GRID_TYPE_COLLECTION_ASYNC,    &
          XDMF_GRID_TYPE_TREE,                &
          XDMF_GRID_TYPE_SUBSET,              &
          XDMF_GRID_COLLECTION_TYPE_SPATIAL,  &
          XDMF_GRID_COLLECTION_TYPE_TEMPORAL, &
          XDMF_GRID_SECTION_DATAITEM,         &
          XDMF_GRID_SECTION_ALL,              &
          XDMF_TIME_TYPE_SINGLE,              &
          XDMF_TIME_TYPE_HYPERSLAB,           &
          XDMF_TIME_TYPE_LIST,                &
          XDMF_TIME_TYPE_RANGE,               &
          XDMF_TOPOLOGY_TYPE_3DSMESH,         &
          XDMF_TOPOLOGY_TYPE_3DRECTMESH,      &
          XDMF_TOPOLOGY_TYPE_3DCORECTMESH,    &
          XDMF_TOPOLOGY_TYPE_2DSMESH,         &
          XDMF_TOPOLOGY_TYPE_2DRECTMESH,      &
          XDMF_TOPOLOGY_TYPE_2DCORECTMESH,    &
          XDMF_TOPOLOGY_TYPE_POLYVERTEX,      &
          XDMF_TOPOLOGY_TYPE_POLYLINE,        &
          XDMF_TOPOLOGY_TYPE_POLYGON,         &
          XDMF_TOPOLOGY_TYPE_TRIANGLE,        &
          XDMF_TOPOLOGY_TYPE_QUADRILATERAL,   &
          XDMF_TOPOLOGY_TYPE_TETRAHEDRON,     &
          XDMF_TOPOLOGY_TYPE_PYRAMID,         &
          XDMF_TOPOLOGY_TYPE_WEDGE,           &
          XDMF_TOPOLOGY_TYPE_HEXAHEDRON,      &
          XDMF_TOPOLOGY_TYPE_EDGE_3,          &
          XDMF_TOPOLOGY_TYPE_TRIANGLE_6,      &
          XDMF_TOPOLOGY_TYPE_QUADRILATERAL_8, &
          XDMF_TOPOLOGY_TYPE_TETRAHEDRON_10,  &
          XDMF_TOPOLOGY_TYPE_PYRAMID_13,      &
          XDMF_TOPOLOGY_TYPE_WEDGE_15,        &
          XDMF_TOPOLOGY_TYPE_HEXAHEDRON_20,   &
          XDMF_TOPOLOGY_TYPE_MIXED
public :: xh5f_file_object,             &
          XH5F_BLOCK_CARTESIAN,         &
          XH5F_BLOCK_CARTESIAN_UNIFORM, &
          XH5F_BLOCK_CURVILINEAR
endmodule motion
