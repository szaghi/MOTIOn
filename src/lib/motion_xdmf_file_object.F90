!< MOTIOn, XDMF file object class.
module motion_xdmf_file_object
!< MOTIOn, XDMF file object class.

use foxy
use penf
use stringifor
use mpi

implicit none
private
public :: xdmf_file_object
public :: XDMF_PARAMETERS

character(1), parameter :: NL = new_line('a') !< New line (end record) character.

type :: xdmf_parameters_object
   !< Global named constants (paramters) class (container) of XDMF syntax.
   character(4)  :: XDMF_ATTR_CENTER_NODE              = 'Node'            !< XDMF attribute node  centered.
   character(4)  :: XDMF_ATTR_CENTER_EDGE              = 'Edge'            !< XDMF attribute edge  centered.
   character(4)  :: XDMF_ATTR_CENTER_FACE              = 'Face'            !< XDMF attribute face  centered.
   character(4)  :: XDMF_ATTR_CENTER_CELL              = 'Cell'            !< XDMF attribute cell  centered.
   character(4)  :: XDMF_ATTR_CENTER_GRID              = 'Grid'            !< XDMF attribute grid  centered.
   character(5)  :: XDMF_ATTR_CENTER_OTHER             = 'Other'           !< XDMF attribute other centered.
   character(6)  :: XDMF_ATTR_TYPE_SCALAR              = 'Scalar'          !< XDMF attribute type scalar.
   character(6)  :: XDMF_ATTR_TYPE_VECTOR              = 'Vector'          !< XDMF attribute type vector.
   character(6)  :: XDMF_ATTR_TYPE_TENSOR              = 'Tensor'          !< XDMF attribute type tensor.
   character(7)  :: XDMF_ATTR_TYPE_TENSOR6             = 'Tensor6'         !< XDMF attribute type tensor6.
   character(6)  :: XDMF_ATTR_TYPE_MATRIX              = 'Matrix'          !< XDMF attribute type matrix.
   character(8)  :: XDMF_ATTR_TYPE_GLOBALID            = 'GlobalID'        !< XDMF attribute type global ID.
   character(2)  :: XDMF_ATTR_ELEM_FAMILY_CG           = 'CG'              !< XDMF attribute element family CG.
   character(2)  :: XDMF_ATTR_ELEM_FAMILY_DG           = 'DG'              !< XDMF attribute element family DG.
   character(1)  :: XDMF_ATTR_ELEM_FAMILY_Q            = 'Q'               !< XDMF attribute element family Q.
   character(2)  :: XDMF_ATTR_ELEM_FAMILY_DQ           = 'DQ'              !< XDMF attribute element family DQ.
   character(2)  :: XDMF_ATTR_ELEM_FAMILY_RT           = 'RT'              !< XDMF attribute element family RT.
   character(8)  :: XDMF_ATTR_ELEM_CELL_INTERVAL       = 'interval'        !< XDMF attribute element cell interval.
   character(8)  :: XDMF_ATTR_ELEM_CELL_TRIANGLE       = 'triangle'        !< XDMF attribute element cell triangle.
   character(11) :: XDMF_ATTR_ELEM_CELL_TETRAHEDRON    = 'tetrahedron'     !< XDMF attribute element cell tetrahedron.
   character(13) :: XDMF_ATTR_ELEM_CELL_QUADRILATERAL  = 'quadrilateral'   !< XDMF attribute element cell quadrilateral.
   character(10) :: XDMF_ATTR_ELEM_CELL_HEXAHEDRON     = 'hexahedron'      !< XDMF attribute element cell hexahedron.
   character(7)  :: XDMF_DATAITEM_ITEMTYPE_UNIFORM     = 'Uniform'         !< XDMF dataitem item type Uniform
   character(10) :: XDMF_DATAITEM_ITEMTYPE_COLLECTION  = 'Collection'      !< XDMF dataitem item type collection
   character(4)  :: XDMF_DATAITEM_ITEMTYPE_TREE        = 'tree'            !< XDMF dataitem item type tree
   character(9)  :: XDMF_DATAITEM_ITEMTYPE_HYPERSLAB   = 'HyperSlab'       !< XDMF dataitem item type hyperSlab
   character(11) :: XDMF_DATAITEM_ITEMTYPE_COORDINATES = 'coordinates'     !< XDMF dataitem item type coordinates
   character(8)  :: XDMF_DATAITEM_ITEMTYPE_FUNCTION    = 'Function'        !< XDMF dataitem item type function
   character(3)  :: XDMF_DATAITEM_NUMBER_FORMAT_HDF    = 'HDF'             !< XDMF dataitem number format HDF.
   character(3)  :: XDMF_DATAITEM_NUMBER_FORMAT_XML    = 'XML'             !< XDMF dataitem number format XML.
   character(6)  :: XDMF_DATAITEM_NUMBER_FORMAT_BINARY = 'Binary'          !< XDMF dataitem number format Binary.
   character(5)  :: XDMF_DATAITEM_NUMBER_TYPE_FLOAT    = 'Float'           !< XDMF dataitem number type float.
   character(3)  :: XDMF_DATAITEM_NUMBER_TYPE_INT      = 'Int'             !< XDMF dataitem number type int.
   character(4)  :: XDMF_DATAITEM_NUMBER_TYPE_UINT     = 'UInt'            !< XDMF dataitem number type uInt.
   character(4)  :: XDMF_DATAITEM_NUMBER_TYPE_CHAR     = 'Char'            !< XDMF dataitem number type char.
   character(5)  :: XDMF_DATAITEM_NUMBER_TYPE_UCHAR    = 'UChar'           !< XDMF dataitem number type uChar.
   character(1)  :: XDMF_DATAITEM_NUMBER_PRECISION_1   = '1'               !< XDMF dataitem number precision 1.
   character(1)  :: XDMF_DATAITEM_NUMBER_PRECISION_2   = '2'               !< XDMF dataitem number precision 2.
   character(1)  :: XDMF_DATAITEM_NUMBER_PRECISION_4   = '4'               !< XDMF dataitem number precision 4.
   character(1)  :: XDMF_DATAITEM_NUMBER_PRECISION_8   = '8'               !< XDMF dataitem number precision 8.
   character(6)  :: XDMF_DATAITEM_ENDIAN_NATIVE        = 'Native'          !< XDMF dataitem number endian native.
   character(3)  :: XDMF_DATAITEM_ENDIAN_BIG           = 'Big'             !< XDMF dataitem number endian big.
   character(6)  :: XDMF_DATAITEM_ENDIAN_LITTLE        = 'Little'          !< XDMF dataitem number endian little.
   character(3)  :: XDMF_DATAITEM_COMPRESSION_RAW      = 'Raw'             !< XDMF dataitem number compression raw.
   character(4)  :: XDMF_DATAITEM_COMPRESSION_ZLIB     = 'Zlib'            !< XDMF dataitem number compression zlib.
   character(5)  :: XDMF_DATAITEM_COMPRESSION_BZIP2    = 'BZip2'           !< XDMF dataitem number compression bzip2.
   character(3)  :: XDMF_GEOMETRY_TYPE_XYZ             = 'XYZ'             !< XDMF geometry type xyz (interlaced values).
   character(2)  :: XDMF_GEOMETRY_TYPE_XY              = 'XY'              !< XDMF geometry type xy.
   character(5)  :: XDMF_GEOMETRY_TYPE_X_Y_Z           = 'X_Y_Z'           !< XDMF geometry type xyz (separated values).
   character(5)  :: XDMF_GEOMETRY_TYPE_VXVYVZ          = 'VXVYVZ'          !< XDMF geometry type xyz (3 arrays).
   character(13) :: XDMF_GEOMETRY_TYPE_ODXYZ           = 'ORIGIN_DXDYDZ'   !< XDMF geometry type origin-dxyz.
   character(11) :: XDMF_GEOMETRY_TYPE_ODXY            = 'ORIGIN_DXDY'     !< XDMF geometry type origin-dxy.
   character(7)  :: XDMF_GRID_TYPE_UNIFORM             = 'Uniform'         !< XDMF grid type uniform.
   character(10) :: XDMF_GRID_TYPE_COLLECTION          = 'Collection'      !< XDMF grid type collection.
   character(15) :: XDMF_GRID_TYPE_COLLECTION_ASYNC    = 'CollectionAsync' !< XDMF grid type collection.
   character(4)  :: XDMF_GRID_TYPE_TREE                = 'Tree'            !< XDMF grid type tree.
   character(6)  :: XDMF_GRID_TYPE_SUBSET              = 'Subset'          !< XDMF grid type subset.
   character(7)  :: XDMF_GRID_COLLECTION_TYPE_SPATIAL  = 'Spatial'         !< XDMF grid collection type spatial.
   character(8)  :: XDMF_GRID_COLLECTION_TYPE_TEMPORAL = 'Temporal'        !< XDMF grid collection type spatial.
   character(8)  :: XDMF_GRID_SECTION_DATAITEM         = 'DataItem'        !< XDMF grid section dataitem.
   character(3)  :: XDMF_GRID_SECTION_ALL              = 'All'             !< XDMF grid section all.
   character(6)  :: XDMF_TIME_TYPE_SINGLE              = 'Single'          !< XDMF time type single.
   character(9)  :: XDMF_TIME_TYPE_HYPERSLAB           = 'HyperSlab'       !< XDMF time type hyperslab.
   character(4)  :: XDMF_TIME_TYPE_LIST                = 'List'            !< XDMF time type list.
   character(5)  :: XDMF_TIME_TYPE_RANGE               = 'Range'           !< XDMF time type range.
   character(7)  :: XDMF_TOPOLOGY_TYPE_3DSMESH         = '3DSMesh'         !< XDMF topology type curvilinear mesh, 3D.
   character(10) :: XDMF_TOPOLOGY_TYPE_3DRECTMESH      = '3DRectMesh'      !< XDMF topology type cartesian mesh, 3D.
   character(12) :: XDMF_TOPOLOGY_TYPE_3DCORECTMESH    = '3DCoRectMesh'    !< XDMF topology type cart uniform mesh, 3D.
   character(7)  :: XDMF_TOPOLOGY_TYPE_2DSMESH         = '2DSMesh'         !< XDMF topology type curvilinear mesh, 2D.
   character(10) :: XDMF_TOPOLOGY_TYPE_2DRECTMESH      = '2DRectMesh'      !< XDMF topology type cartesian mesh, 2D.
   character(12) :: XDMF_TOPOLOGY_TYPE_2DCORECTMESH    = '2DCoRectMesh'    !< XDMF topology type cart uniform mesh, 2D.
   character(10) :: XDMF_TOPOLOGY_TYPE_POLYVERTEX      = 'Polyvertex'      !< XDMF topology type polyvertex.
   character(8)  :: XDMF_TOPOLOGY_TYPE_POLYLINE        = 'Polyline'        !< XDMF topology type polyline.
   character(7)  :: XDMF_TOPOLOGY_TYPE_POLYGON         = 'Polygon'         !< XDMF topology type polygon.
   character(8)  :: XDMF_TOPOLOGY_TYPE_TRIANGLE        = 'Triangle'        !< XDMF topology type triangle.
   character(13) :: XDMF_TOPOLOGY_TYPE_QUADRILATERAL   = 'Quadrilateral'   !< XDMF topology type quadrilateral.
   character(11) :: XDMF_TOPOLOGY_TYPE_TETRAHEDRON     = 'Tetrahedron'     !< XDMF topology type tetrahedron.
   character(7)  :: XDMF_TOPOLOGY_TYPE_PYRAMID         = 'Pyramid'         !< XDMF topology type pyramid.
   character(5)  :: XDMF_TOPOLOGY_TYPE_WEDGE           = 'Wedge'           !< XDMF topology type wedge.
   character(10) :: XDMF_TOPOLOGY_TYPE_HEXAHEDRON      = 'Hexahedron'      !< XDMF topology type hexahedron.
   character(6)  :: XDMF_TOPOLOGY_TYPE_EDGE_3          = 'Edge_3'          !< XDMF topology type edge_3.
   character(10) :: XDMF_TOPOLOGY_TYPE_TRIANGLE_6      = 'Triangle_6'      !< XDMF topology type triangle_6.
   character(15) :: XDMF_TOPOLOGY_TYPE_QUADRILATERAL_8 = 'Quadrilateral_8' !< XDMF topology type quadrilateral_8.
   character(14) :: XDMF_TOPOLOGY_TYPE_TETRAHEDRON_10  = 'Tetrahedron_10'  !< XDMF topology type tetrahedron_10.
   character(10) :: XDMF_TOPOLOGY_TYPE_PYRAMID_13      = 'Pyramid_13'      !< XDMF topology type pyramid_13.
   character(8)  :: XDMF_TOPOLOGY_TYPE_WEDGE_15        = 'Wedge_15'        !< XDMF topology type wedge_15.
   character(13) :: XDMF_TOPOLOGY_TYPE_HEXAHEDRON_20   = 'Hexahedron_20'   !< XDMF topology type hexahedron_20.
   character(5)  :: XDMF_TOPOLOGY_TYPE_MIXED           = 'Mixed'           !< XDMF topology type mixed.
endtype xdmf_parameters_object
type(xdmf_parameters_object), parameter :: XDMF_PARAMETERS=xdmf_parameters_object() !< List of XDMF named constants.

type :: xdmf_file_object
   !< XDMF file object class.
   type(string)  :: filename           !< File name.
   integer(I4P)  :: indent=0_I4P       !< Indent count.
   integer(I4P)  :: xml=0_I4P          !< XML Logical unit.
   integer(I4P)  :: error=0_I4P        !< Error status.
   type(xml_tag) :: tag                !< XML tags handler.
   integer(I4P)  :: procs_number=1_I4P !< Number of MPI processes.
   integer(I4P)  :: myrank=0_I4P       !< MPI ID process.
   logical       :: is_async=.false.   !< Asyncronous saving.
   type(string)  :: async_tags         !< Asyncronous tags data.
   contains
      ! public methods
      ! file methods
      procedure, pass(self) :: close_file !< Close XDMF file.
      procedure, pass(self) :: open_file  !< Open XDMF file.
      ! XDMF tag methods
      ! attribute tag
      procedure, pass(self) :: close_attribute_tag !< Close `Attribute` tag.
      procedure, pass(self) :: open_attribute_tag  !< Open `Attribute` tag.
      ! dataitem tag
      procedure, pass(self) :: dataitem_tag_attr  !< Return `DataItem` tag attributes as string.
      ! procedure, pass(self) :: dataitem_tag_str   !< Return `DataItem` tag as string.
      procedure, pass(self) :: write_dataitem_tag !< Write `DataItem` tag.
      ! domain tag
      procedure, pass(self) :: close_domain_tag !< Close `Domain` tag.
      procedure, pass(self) :: open_domain_tag  !< Open `Domain` tag.
      ! geometry tag
      procedure, pass(self) :: close_geometry_tag !< Close `Geometry` tag.
      procedure, pass(self) :: open_geometry_tag  !< Open `Geometry` tag.
      ! grid tag
      procedure, pass(self) :: close_grid_tag !< Close `Grid` tag.
      procedure, pass(self) :: open_grid_tag  !< Open `Grid` tag.
      ! async tags
      procedure, pass(self) :: write_async_tags !< Write async tags of other (not mine...) processes.
      ! header tag
      procedure, pass(self) :: write_header_tag !< Write header tag.
      ! time tag
      procedure, pass(self) :: write_time_tag !< Write `Time` tag.
      ! topology tag
      procedure, pass(self) :: write_topology_tag !< Write `Topology` tag.
      ! generic tag methods
      procedure, pass(self) :: tag_str                !< Return tag as string.
      procedure, pass(self) :: write_end_tag          !< Write `</tag_name>` end tag.
      procedure, pass(self) :: write_self_closing_tag !< Write self closing tag.
      procedure, pass(self) :: write_start_tag        !< Write start tag.
      procedure, pass(self) :: write_tag              !< Write tag.
      ! MPI methods
      procedure, pass(self) :: gather_async_tags !< Gather async tags.
endtype xdmf_file_object

contains
   ! files methods
   subroutine close_file(self)
   !< Close XDMF file.
   class(xdmf_file_object), intent(inout) :: self !< File handler.

   if (.not.self%is_async) then
      call self%write_end_tag(name='Xdmf')
      close(unit=self%xml, iostat=self%error)
   endif
   endsubroutine close_file

   subroutine open_file(self, filename)
   !< Open XDMF file.
   class(xdmf_file_object), intent(inout) :: self               !< File handler.
   character(*),            intent(in)    :: filename           !< File name.
   logical                                :: is_mpi_initialized !< MPI env status.

   ! reset file handler
   select type(self)
   type is(xdmf_file_object)
      self = xdmf_file_object()
   endselect
   self%is_async = .false. ;
   self%async_tags = ''
   self%filename = trim(adjustl(filename))
   call MPI_INITIALIZED(is_mpi_initialized, self%error)
   if (.not.is_mpi_initialized) call MPI_INIT(self%error)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, self%procs_number, self%error)
   call MPI_COMM_RANK(MPI_COMM_WORLD, self%myrank, self%error)
   if (self%myrank/=0_I4P) self%is_async = .true.
   if (.not.self%is_async) open(newunit=self%xml,           &
                                file=self%filename%chars(), &
                                form='UNFORMATTED',         &
                                access='STREAM',            &
                                action='WRITE',             &
                                status='REPLACE',           &
                                iostat=self%error)
   call self%write_header_tag
   endsubroutine open_file

   ! XDMF tag methods
   ! attribute tag
   subroutine close_attribute_tag(self)
   !< Close `Attribute` tag.
   class(xdmf_file_object), intent(inout) :: self !< File handler.

   call self%write_end_tag(name='Attribute')
   endsubroutine close_attribute_tag

   subroutine open_attribute_tag(self,attribute_name,attribute_center,attribute_type,attribute_item_type, &
                                 element_family,element_degree,element_cell)
   !< Open `Attribute` tag.
   !< @NOTE Attribute tag has the following attributes:
   !< ```
   !< ->            default value
   !< Name          (no default)
   !< Center        Node | Edge | Face | Cell | Grid | Other
   !< AttributeType Scalar | Vector | Tensor (9 values expected) | Tensor6 (summetric tensor) | Matrix (arbitrary NxM matrix)
   !< -> Only Meaningful if ItemType="FiniteElementFunction"
   !< ItemType      (no default) | FiniteElementFunction
   !< ElementFamily (no default) | CG | DG | Q | DQ | RT
   !< ElementDegree (no default) Arbitrary integer value
   !< ElementCell   (no default) | interval | triangle | tetrahedron | quadrilateral | hexahedron
   !< ```
   class(xdmf_file_object), intent(inout)        :: self                !< File handler.
   character(*),            intent(in), optional :: attribute_name      !< Attribute name.
   character(*),            intent(in), optional :: attribute_center    !< Attribute center.
   character(*),            intent(in), optional :: attribute_type      !< Attribute type.
   character(*),            intent(in), optional :: attribute_item_type !< Attribute item type.
   character(*),            intent(in), optional :: element_family      !< Attribute element family.
   character(*),            intent(in), optional :: element_degree      !< Attribute element degree.
   character(*),            intent(in), optional :: element_cell        !< Attribute element cell.
   character(:), allocatable                     :: attr                !< Attributes list.

   attr = '' ; if (present(attribute_name     )) attr =        'Name="'         //trim(adjustl(attribute_name     ))//'"'
               if (present(attribute_center   )) attr = attr//' Center="'       //trim(adjustl(attribute_center   ))//'"'
               if (present(attribute_type     )) attr = attr//' AttributeType="'//trim(adjustl(attribute_type     ))//'"'
               if (present(attribute_item_type)) attr = attr//' ItemType="'     //trim(adjustl(attribute_item_type))//'"'
               if (present(element_family     )) attr = attr//' ElementFamily="'//trim(adjustl(element_family     ))//'"'
               if (present(element_degree     )) attr = attr//' ElementDegree="'//trim(adjustl(element_degree     ))//'"'
               if (present(element_cell       )) attr = attr//' ElementCell="'  //trim(adjustl(element_cell       ))//'"'
   call self%write_start_tag(name='Attribute', attributes=attr)
   endsubroutine open_attribute_tag

   ! dataitem tag
   function dataitem_tag_attr(self,             &
                              item_type,        &
                              item_dimensions,  &
                              number_type,      &
                              number_precision, &
                              number_format,    &
                              endian,           &
                              compression,      &
                              seek) result(attr)
   !< Return `DataItem` tag attributes as string.
   !< @NOTE DataItem tag has the following attributes:
   !< ```
   !< ->          default value
   !< Name        (no default)
   !< ItemType    Uniform | Collection | tree | HyperSlab | coordinates | Function
   !< Dimensions  (no default) in KJI Order
   !< NumberType  Float | Int | UInt | Char | UChar
   !< Precision   4 | 1 | 2 (Int or UInt only) | 8
   !< Format      XML | HDF | Binary
   !< Endian      Native | Big | Little (applicable only to Binary format)
   !< Compression Raw|Zlib|BZip2 (applicable only to Binary format and depend on xdmf configuration)
   !< Seek        0 (number of bytes to skip, applicable only to Binary format with Raw compression)
   !< ```
   class(xdmf_file_object), intent(inout)        :: self             !< File handler.
   character(*),            intent(in), optional :: item_type        !< Item type.
   character(*),            intent(in), optional :: item_dimensions  !< Item dimensions, external given.
   character(*),            intent(in), optional :: number_type      !< Number type.
   character(*),            intent(in), optional :: number_precision !< Number precision.
   character(*),            intent(in), optional :: number_format    !< Number format.
   character(*),            intent(in), optional :: endian           !< Bits endian.
   character(*),            intent(in), optional :: compression      !< Data compression.
   character(*),            intent(in), optional :: seek             !< Bytes to skip.
   character(:), allocatable                     :: attr             !< Attributes list.

   attr = '' ; if (present(item_type       )) attr =        'ItemType="'    //trim(adjustl(item_type       ))//'"'
               if (present(item_dimensions )) attr = attr//' Dimensions="'  //trim(adjustl(item_dimensions ))//'"'
               if (present(number_type     )) attr = attr//' DataType="'    //trim(adjustl(number_type     ))//'"'
               if (present(number_precision)) attr = attr//' Precision="'   //trim(adjustl(number_precision))//'"'
               if (present(number_format   )) attr = attr//' Format="'      //trim(adjustl(number_format   ))//'"'
               if (present(endian          )) attr = attr//' Endian="'      //trim(adjustl(endian          ))//'"'
               if (present(compression     )) attr = attr//' Compression="' //trim(adjustl(compression     ))//'"'
               if (present(seek            )) attr = attr//' Seek="'        //trim(adjustl(seek            ))//'"'
   endfunction dataitem_tag_attr

   function dataitem_tag_str(self,             &
                             content,          &
                             item_type,        &
                             item_dimensions,  &
                             number_type,      &
                             number_precision, &
                             number_format,    &
                             endian,           &
                             compression,      &
                             seek)
   !< Return `DataItem` tag as string.
   class(xdmf_file_object), intent(inout)        :: self             !< File handler.
   character(*),            intent(in)           :: content          !< Data item content.
   character(*),            intent(in), optional :: item_type        !< Item type.
   character(*),            intent(in), optional :: item_dimensions  !< Item dimensions, external given.
   character(*),            intent(in), optional :: number_type      !< Number type.
   character(*),            intent(in), optional :: number_precision !< Number precision.
   character(*),            intent(in), optional :: number_format    !< Number format.
   character(*),            intent(in), optional :: endian           !< Bits endian.
   character(*),            intent(in), optional :: compression      !< Data compression.
   character(*),            intent(in), optional :: seek             !< Bytes to skip.
   character(:), allocatable                     :: dataitem_tag_str !< Dataitem tag as string.
   character(:), allocatable                     :: attr             !< Attributes list.

   attr = self%dataitem_tag_attr(item_type        = item_type,        &
                                 item_dimensions  = item_dimensions,  &
                                 number_type      = number_type,      &
                                 number_precision = number_precision, &
                                 number_format    = number_format,    &
                                 endian           = endian,           &
                                 compression      = compression,      &
                                 seek             = seek)
   dataitem_tag_str = self%tag_str(name='DataItem', attributes=attr, content=content)
   endfunction dataitem_tag_str

   subroutine write_dataitem_tag(self,             &
                                 content,          &
                                 item_type,        &
                                 item_dimensions,  &
                                 number_type,      &
                                 number_precision, &
                                 number_format,    &
                                 endian,           &
                                 compression,      &
                                 seek)
   !< Write `DataItem` tag.
   class(xdmf_file_object), intent(inout)        :: self             !< File handler.
   character(*),            intent(in)           :: content          !< Data item content.
   character(*),            intent(in), optional :: item_type        !< Item type.
   character(*),            intent(in), optional :: item_dimensions  !< Item dimensions, external given.
   character(*),            intent(in), optional :: number_type      !< Number type.
   character(*),            intent(in), optional :: number_precision !< Number precision.
   character(*),            intent(in), optional :: number_format    !< Number format.
   character(*),            intent(in), optional :: endian           !< Bits endian.
   character(*),            intent(in), optional :: compression      !< Data compression.
   character(*),            intent(in), optional :: seek             !< Bytes to skip.
   character(:), allocatable                     :: attr             !< Attributes list.

   attr = self%dataitem_tag_attr(item_type        = item_type,        &
                                 item_dimensions  = item_dimensions,  &
                                 number_type      = number_type,      &
                                 number_precision = number_precision, &
                                 number_format    = number_format,    &
                                 endian           = endian,           &
                                 compression      = compression,      &
                                 seek             = seek)
   call self%write_tag(name='DataItem', attributes=attr, content=content)
   endsubroutine write_dataitem_tag

   ! domain tag
   subroutine close_domain_tag(self)
   !< Close `Domain` tag.
   class(xdmf_file_object), intent(inout) :: self !< File handler.

   if (.not.self%is_async) then
      call self%write_end_tag(name='Domain')
   else
      ! in asyncronous mode do not for Domain tag, only update indent number
      self%indent = self%indent - 2
   endif
   endsubroutine close_domain_tag

   subroutine open_domain_tag(self)
   !< Open `Domain` tag.
   class(xdmf_file_object), intent(inout) :: self !< File handler.

   if (.not.self%is_async) then
      call self%write_start_tag(name='Domain')
   else
      ! in asyncronous mode do not for Domain tag, only update indent number
      self%indent = self%indent + 2
   endif
   endsubroutine open_domain_tag

   ! geometry tag
   subroutine close_geometry_tag(self)
   !< Close `Geometry` tag.
   class(xdmf_file_object), intent(inout) :: self !< File handler.

   call self%write_end_tag(name='Geometry')
   endsubroutine close_geometry_tag

   subroutine open_geometry_tag(self, geometry_type)
   !< Open `Geometry` tag.
   !< @NOTE Geometry type has the following values:
   !< ```
   !< XYZ - Interlaced locations (default)
   !< XY - Z is set to 0.0
   !< X_Y_Z - X,Y, and Z are separate arrays
   !< VXVYVZ - Three arrays, one for each axis
   !< ORIGIN_DXDYDZ - Six Values : Ox,Oy,Oz + Dx,Dy,Dz
   !< ORIGIN_DXDY - Four Values : Ox,Oy + Dx,Dy
   !< ```
   class(xdmf_file_object), intent(inout)        :: self           !< File handler.
   character(*),            intent(in), optional :: geometry_type  !< Geometry type.
   character(:), allocatable                     :: attr           !< Attributes list.

   attr = '' ; if (present(geometry_type)) attr = 'Type="'//trim(adjustl(geometry_type))//'"'
   call self%write_start_tag(name='Geometry', attributes=attr)
   endsubroutine open_geometry_tag

   ! grid tag
   subroutine close_grid_tag(self, grid_type)
   !< Close `Grid` tag.
   class(xdmf_file_object), intent(inout)        :: self       !< File handler.
   character(*),            intent(in), optional :: grid_type  !< Grid type.
   character(:), allocatable                     :: grid_type_ !< Grid type, local var.

   grid_type_ = '' ; if (present(grid_type)) grid_type_ = trim(adjustl(grid_type))
   if (grid_type_==XDMF_PARAMETERS%XDMF_GRID_TYPE_COLLECTION_ASYNC) then
      if (self%is_async) then
         ! close asyncrounous collection of async MPI process, update indent
         self%indent = self%indent - 2
      endif
      ! close asyncrounous collection of master MPI process, gather async tags before closing
      call self%gather_async_tags
      if (self%myrank==0_I4P) call self%write_async_tags(async_tags=self%async_tags)
   endif
   call self%write_end_tag(name='Grid')
   endsubroutine close_grid_tag

   subroutine open_grid_tag(self, grid_name, grid_type, grid_collection_type, grid_section)
   !< Open `Grid` tag.
   !< @NOTE Grid tag has the following attributes:
   !< ```
   !< ->              default value
   !< Name            (no default)
   !< GridType        Uniform | Collection | Tree | Subset
   !< CollectionType  Spatial | Temporal (Only Meaningful if GridType="Collection")
   !< Section         DataItem | All  (Only Meaningful if GridType="Subset")
   !< ```
   class(xdmf_file_object), intent(inout)        :: self                 !< File handler.
   character(*),            intent(in), optional :: grid_name            !< Grid name.
   character(*),            intent(in), optional :: grid_type            !< Grid type.
   character(*),            intent(in), optional :: grid_collection_type !< Grid collection type.
   character(*),            intent(in), optional :: grid_section         !< Grid section.
   character(:), allocatable                     :: grid_type_           !< Grid type, local var.
   character(:), allocatable                     :: attr                 !< Attributes list.

   grid_type_ = '' ; if (present(grid_type)) grid_type_ = trim(adjustl(grid_type))
   if (grid_type_==XDMF_PARAMETERS%XDMF_GRID_TYPE_COLLECTION_ASYNC) then
      if (self%is_async) then
         ! asyncrounous collection of async MPI process, just update indent
         self%indent = self%indent + 2
         return
      else
         ! asyncrounous collection of master MPI process, convert grid type to standard collection
         grid_type_ = XDMF_PARAMETERS%XDMF_GRID_TYPE_COLLECTION
      endif
   endif
   attr = '' ; if (present(grid_name           )) attr =        'Name="'          //trim(adjustl(grid_name           ))//'"'
               if (present(grid_type           )) attr = attr//' GridType="'      //trim(adjustl(grid_type_          ))//'"'
               if (present(grid_collection_type)) attr = attr//' CollectionType="'//trim(adjustl(grid_collection_type))//'"'
               if (present(grid_section        )) attr = attr//' Section="'       //trim(adjustl(grid_section        ))//'"'
   call self%write_start_tag(name='Grid', attributes=attr)
   endsubroutine open_grid_tag

   ! async tags
   subroutine write_async_tags(self, async_tags)
   !< Write async tags of other (not mine...) processes.
   class(xdmf_file_object), intent(inout) :: self       !< File handler.
   type(string),            intent(in)    :: async_tags !< Asyncronous tags data.

   write(unit=self%xml, iostat=self%error) async_tags%chars()
   endsubroutine write_async_tags

   ! header tag
   subroutine write_header_tag(self)
   !< Write header tag.
   class(xdmf_file_object), intent(inout) :: self   !< File handler.
   type(string)                           :: buffer !< Buffer string.

   if (.not.self%is_async) then
      buffer =         '<?xml version="1.0" encoding="utf-8"?>'//NL
      buffer = buffer//'<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="3.0">'
      write(unit=self%xml, iostat=self%error) buffer//NL
   endif
   self%indent = 2
   endsubroutine write_header_tag

   ! topolgy tag
   subroutine write_time_tag(self,time_value,time_type,time_dimensions,time_number_type,time_number_precision,time_number_format)
   !< Write `Time` tag.
   !< @NOTE Time tag has the following attributes:
   !< ```
   !< ->        default value
   !< TimeType  Single | HyperSlab | List | Range
   !< Value     (no default - Only valid for TimeType="Single")
   !< ```
   class(xdmf_file_object), intent(inout)        :: self                  !< File handler.
   character(*),            intent(in)           :: time_value            !< Time value.
   character(*),            intent(in), optional :: time_type             !< Time type.
   character(*),            intent(in), optional :: time_dimensions       !< Time (list) dimensions.
   character(*),            intent(in), optional :: time_number_type      !< Time number type.
   character(*),            intent(in), optional :: time_number_precision !< Time number precision.
   character(*),            intent(in), optional :: time_number_format    !< Time number format.
   character(:), allocatable                     :: time_type_            !< Time type, local var.

   time_type_ = XDMF_PARAMETERS%XDMF_TIME_TYPE_SINGLE ; if (present(time_type)) time_type_ = trim(adjustl(time_type))
   if (time_type_==XDMF_PARAMETERS%XDMF_TIME_TYPE_SINGLE) then
      call self%write_self_closing_tag(name='Time', attributes='Value="'//trim(adjustl(time_value))//'"')
   else
      call self%write_start_tag(name='Time', attributes='TimeType="'//time_type_//'"')
      call self%write_dataitem_tag(content          = time_value,            &
                                   item_dimensions  = time_dimensions,       &
                                   number_type      = time_number_type,      &
                                   number_precision = time_number_precision, &
                                   number_format    = time_number_format)
      call self%write_end_tag(name='Time')
   endif
   endsubroutine write_time_tag

   ! topolgy tag
   subroutine write_topology_tag(self, topology_type, topology_dimensions)
   !< Write `Topology` tag.
   !< @NOTE Topology tag has the following attributes:
   !< ```
   !< ->               default value
   !< Name             (no default)
   !< TopologyType     Polyvertex | Polyline | Polygon |
   !<                  Triangle | Quadrilateral | Tetrahedron | Pyramid| Wedge | Hexahedron |
   !<                  Edge_3 | Triangle_6 | Quadrilateral_8 | Tetrahedron_10 | Pyramid_13 |
   !<                  Wedge_15 | Hexahedron_20 |
   !<                  Mixed |
   !<                  2DSMesh | 2DRectMesh | 2DCoRectMesh |
   !<                  3DSMesh | 3DRectMesh | 3DCoRectMesh
   !< NodesPerElement  (no default) Only Important for Polyvertex, Polygon and Polyline
   !< NumberOfElements (no default)
   !<     OR
   !< Dimensions       (no default)
   !< Order            each cell type has its own default
   !< ```
   class(xdmf_file_object), intent(inout) :: self                !< File handler.
   character(*),            intent(in)    :: topology_type       !< Topology type.
   character(*),            intent(in)    :: topology_dimensions !< Topology dimensions.

   call self%write_self_closing_tag(name='Topology', &
                                    attributes='Type="'//trim(adjustl(topology_type))//'" '//&
                                               'Dimensions="'//trim(adjustl(topology_dimensions))//'"')
   endsubroutine write_topology_tag

   ! generic tag methods
   function tag_str(self, name, attributes, content)
   !< Return tag as string
   class(xdmf_file_object), intent(inout)        :: self       !< File handler.
   character(*),            intent(in)           :: name       !< Tag name.
   character(*),            intent(in), optional :: attributes !< Tag attributes.
   character(*),            intent(in), optional :: content    !< Tag content.
   character(:), allocatable                     :: tag_str    !< Tas as string.

   self%tag = xml_tag(name=name, attributes_stream=attributes, sanitize_attributes_value=.true., content=content, &
                      indent=self%indent)
   tag_str = self%tag%stringify(is_indented=.true., is_content_indented=.true.)
   endfunction tag_str

   subroutine write_end_tag(self, name)
   !< Write `</tag_name>` end tag.
   class(xdmf_file_object), intent(inout) :: self !< File handler.
   character(*),            intent(in)    :: name !< Tag name.

   self%indent = self%indent - 2
   self%tag = xml_tag(name=name, indent=self%indent)
   if (.not.self%is_async) then
      call self%tag%write(unit=self%xml, iostat=self%error, is_indented=.true., end_record=NL, only_end=.true.)
   else
      self%async_tags = self%async_tags//self%tag%stringify(is_indented=.true., only_end=.true.)//NL
   endif
   endsubroutine write_end_tag

   subroutine write_self_closing_tag(self, name, attributes)
   !< Write `<tag_name.../>` self closing tag.
   class(xdmf_file_object), intent(inout)        :: self       !< File handler.
   character(*),            intent(in)           :: name       !< Tag name.
   character(*),            intent(in), optional :: attributes !< Tag attributes.

   self%tag = xml_tag(name=name, attributes_stream=attributes, sanitize_attributes_value=.true., indent=self%indent, &
                      is_self_closing=.true.)
   if (.not.self%is_async) then
      call self%tag%write(unit=self%xml, iostat=self%error, is_indented=.true., end_record=NL)
   else
      self%async_tags = self%async_tags//self%tag%stringify(is_indented=.true.)//NL
   endif
   endsubroutine write_self_closing_tag

   subroutine write_start_tag(self, name, attributes)
   !< Write `<tag_name...>` start tag.
   class(xdmf_file_object), intent(inout)        :: self       !< File handler.
   character(*),            intent(in)           :: name       !< Tag name.
   character(*),            intent(in), optional :: attributes !< Tag attributes.

   self%tag = xml_tag(name=name, attributes_stream=attributes, sanitize_attributes_value=.true., indent=self%indent)
   if (.not.self%is_async) then
      call self%tag%write(unit=self%xml, iostat=self%error, is_indented=.true., end_record=NL, only_start=.true.)
   else
      self%async_tags = self%async_tags//self%tag%stringify(is_indented=.true., only_start=.true.)//NL
   endif
   self%indent = self%indent + 2
   endsubroutine write_start_tag

   subroutine write_tag(self, name, attributes, content)
   !< Write `<tag_name...>...</tag_name>` tag.
   class(xdmf_file_object), intent(inout)        :: self       !< File handler.
   character(*),            intent(in)           :: name       !< Tag name.
   character(*),            intent(in), optional :: attributes !< Tag attributes.
   character(*),            intent(in), optional :: content    !< Tag content.

   self%tag = xml_tag(name=name, attributes_stream=attributes, sanitize_attributes_value=.true., content=content, &
                      indent=self%indent)
   if (.not.self%is_async) then
      call self%tag%write(unit=self%xml, iostat=self%error, is_indented=.true., is_content_indented=.true., end_record=NL)
   else
      self%async_tags = self%async_tags//self%tag%stringify(is_indented=.true., is_content_indented=.true.)//NL
   endif
   endsubroutine write_tag

   ! MPI methods
   subroutine gather_async_tags(self)
   !< Gather async tags.
   !< @NOTE master process (myrank==0) gather the asyncronous tags of other MPI procs and save them
   !< into its own async_tags string.
   class(xdmf_file_object), intent(inout) :: self               !< File handler.
   integer(I4P)                           :: my_async_tags_len  !< Length of my async tags.
   integer(I4P)                           :: all_async_tags_len !< Length of all async tags, total lenght.
   integer(I4P), allocatable              :: recvcounts(:)      !< Size of chars from other processes.
   integer(I4P), allocatable              :: offset(:)          !< Offset in receive buffer.
   character(:), allocatable              :: recvbuf            !< Receive buffer.
   integer(I4P)                           :: i                  !< Counter.

   if (self%procs_number==1) return ! no multi MPI procs, nothing to gather

   call MPI_BARRIER(MPI_COMM_WORLD, self%error) ! all MPI procs must close their XDMF async tags

   ! gather all async tags lengths
   my_async_tags_len = self%async_tags%len()
   if (self%myrank == 0_I4P) then
      allocate(recvcounts(self%procs_number))
      allocate(offset(self%procs_number))
   endif
   call MPI_GATHER(my_async_tags_len, 1_I4P, MPI_INTEGER, recvcounts, 1_I4P, MPI_INTEGER, 0_I4P, MPI_COMM_WORLD, self%error)
   ! compute offset and total lenght of receive buffer
   if (self%myrank == 0_I4P) then
      offset(1) = 0_I4P
      all_async_tags_len = recvcounts(1)
      do i = 2, self%procs_number
         offset(i) = offset(i-1) + recvcounts(i-1)
         all_async_tags_len = all_async_tags_len + recvcounts(i)
      enddo
      allocate(character(len=all_async_tags_len) :: recvbuf)
   endif
   ! gather async tags
   call MPI_GATHERV(self%async_tags%chars(), my_async_tags_len, MPI_CHARACTER, &
                    recvbuf, recvcounts, offset, MPI_CHARACTER, &
                    0_I4P, MPI_COMM_WORLD, self%error)
   if (self%myrank==0_I4P) self%async_tags = recvbuf
   endsubroutine gather_async_tags
endmodule motion_xdmf_file_object
