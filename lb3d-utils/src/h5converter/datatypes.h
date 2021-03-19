//! Data types
enum DataTypes {
  UNDEFINED = -1, //!< Error value
  HDF5      =  0, //!< HDF5 float array (LB3D output)
  XDRF      =  1, //!< XDR float array (LB3D output)
  UINT8     =  2  //!< Unstructured uint8_t data (Holger's CT data)
};

