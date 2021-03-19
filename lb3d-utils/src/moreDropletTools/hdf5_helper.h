#include <stdlib.h>
#include <vector>
#include <string>
#include <iostream>
#include <map>

#ifndef VERBOSITY_LEVEL_H
#define VERBOSITY_LEVEL_H
#include "verbosity_level.h"
#endif

using namespace std;

#ifndef H5_NO_NAMESPACE
#ifndef H5_NO_STD
  using std::cout;
  using std::endl;
#endif  // H5_NO_STD
#endif

#ifndef H5_C_HEADER
#include "H5Cpp.h"
#else
#include "hdf5.h"
#endif

//! Verbosity for the HDF5 wrapper functions.
/*!
    This is set to NONE by default - one can change this by setting it in different parts of the calling code.
*/
extern VerbosityLevel HDF_VERBOSITY;

//! Parser object for an LB3D input file.
/*!
    One can use this class to load an HDF5 output file which contains an input file in its metadata.
*/
class LB3dAttributes {
  vector<string> lines;             /*!< After creating the object, contains all lines of the HDF5 metadata. */
  map<string,string> values;        /*!< After calling LB3dAttributes::parseLines, contains key-values pairs of all parameter declarations found in LB3dAttributes::lines. */

  //! Parses the values stored in LB3dAttributes::lines.
  void parseLines();
  //! Helper function to strip all leading and trailing whitespace from a string.
  /*!
      @param[in,out] str String to be stripped of leading and trailing whitespaces.
  */
  void trimSpaces(string& str);

  public:
    //! Constructor.
    /*!
        @param[in] filename Filename of the HDF5 file.
    */
    LB3dAttributes(const char* filename);

    static const int ok    =  0;    /*!< Return value for the getter functions. */
    static const int error = -1;    /*!< Return value for the getter functions. */

    //! Retrieve a string value from the input file.
    /*!
        @param[in]  key   Name of the parameter.
        @param[out] value Value of the parameter, if found.
        @return           Exit status.
    */
    int getString(char* key, string* value);

    //! Retrieve a boolean value from the input file.
    /*!
        @param[in]  key   Name of the parameter.
        @param[out] value Value of the parameter, if found.
        @return           Exit status.
    */
    int getBool(char* key, bool* value);

    //! Retrieve a double value from the input file.
    /*!
        @param[in]  key   Name of the parameter.
        @param[out] value Value of the parameter, if found.
        @return           Exit status.
    */
    int getDouble(char* key, double* value);

    //! Retrieve an integer value from the input file.
    /*!
        @param[in]  key   Name of the parameter.
        @param[out] value Value of the parameter, if found.
        @return           Exit status.
    */
    int getInteger(char* key, int* value);
};

#ifndef H5_NO_NAMESPACE
  using namespace H5;
#endif

//! Convert dimension array of hsize_t to integers.
/*!
    @param[in]  dims      Array of 3 hsize_t containing dimensions of a system (FORTRAN ORDER!).
    @param[out] nx        Size of system in x-direction.
    @param[out] ny        Size of system in y-direction.
    @param[out] nz        Size of system in z-direction.
*/
void hdfToIntDims(hsize_t dims[], int *nx, int *ny, int *nz);

//! Convert dimension array of hsize_t to integers.
/*!
    @param[in]  dims      Array of 4 hsize_t containing dimensions of a system (FORTRAN ORDER!).
    @param[out] nx        Size of system in x-direction.
    @param[out] ny        Size of system in y-direction.
    @param[out] nz        Size of system in z-direction.
    @param[out] nv        Size of system in v-direction.
*/
void hdfToIntDims(hsize_t dims[], int *nx, int *ny, int *nz, int *nv);

//! Get the dimensions of an HDF5 dataset.
/*!
    @param[in]  filename  Path to the HDF5 file.
    @param[out] pdims     Dimensions of the dataset.
    @param[in]  dsetname  Name of the dataset.
    @return               Exit status.
*/
int  hdfGetDims(char* filename, hsize_t **pdims, const char *dsetname);
//! Get the dimensions from the default dataset std::HDF_DEFAULT_DATASET of rank std::HDF_RANK.
/*!
    \overload
*/
int  hdfGetDims(char* filename, hsize_t **pdims);

//! Read all data from an HDF5 dataset of rank 3.
/*!
    @param[in]  filename  Path to the HDF5 file.
    @param[out] pdata     Array containing the dataset.
    @param[out] pdims     Dimensions of the dataset.
    @param[in]  dsetname  Name of the dataset.
    @return               Exit status.
*/
int  hdfReadAll(char* filename, float ****pdata, hsize_t **pdims, const char* dsetname);
//! Read all data from the default dataset std::HDF_DEFAULT_DATASET of rank 3.
/*!
    \overload
*/
int  hdfReadAll(char* filename, float ****pdata, hsize_t **pdims);

//! Read data from a slab of an HDF5 dataset of rank 3.
/*!
    @param[in]  filename  Path to the HDF5 file.
    @param[out] pdata     Array containing the dataset.
    @param[in]  ox        x-offset of the slab.
    @param[in]  oy        y-offset of the slab.
    @param[in]  oz        z-offset of the slab.
    @param[in]  dx        x-size of the slab.
    @param[in]  dy        y-size of the slab.
    @param[in]  dz        z-size of the slab.
    @param[in]  dsetname  Name of the dataset.
    @return               Exit status.
*/
int  hdfReadSlab(char* filename, float ****pdata, const int ox, const int oy, const int oz, const int dx, const int dy, const int dz, const char *dsetname);
//! Read data from a slab of the default dataset std::HDF_DEFAULT_DATASET of rank 3.
/*!
    \overload
*/
int  hdfReadSlab(char* filename, float ****pdata, const int ox, const int oy, const int oz, const int dx, const int dy, const int dz);

//! Read all data from an HDF5 dataset of rank 4 with offset.
/*!
    @param[in]  filename  Path to the HDF5 file.
    @param[out] pdata     Array containing the dataset.
    @param[out] pdims     Dimensions of the dataset.
    @param[in]  dsetname  Name of the dataset.
    @param[in]  voffset   Vector offset.
    @return               Exit status.
*/
int  hdfReadAllVecComp(char* filename, float **pdata, hsize_t **pdims, const int voffset, const char* dsetname);
//! Read all data from the default dataset std::HDF_DEFAULT_DATASET of rank 4 with offset.
/*!
    \overload
*/
int  hdfReadAllVecComp(char* filename, float **pdata, hsize_t **pdims, const int voffset);

//! Free a rank-3 array allocated by std::hdfReadAll or std::hdfReadSlab.
/*!
    @param[in,out] data   Array to be freed.
    @return               Exit status.
*/
int  hdfFree3DArray(float ***data);
//! Free an array allocated by std::hdfGetDims.
/*!
    @param[in,out] dims   Array to be freed.
    @return               Exit status.
*/
int  hdfFreeDims(hsize_t *dims);
//! Free a rank-4 array allocated by std::hdfReadAllVec.
/*!
    @param[in,out] data   Array to be freed.
    @return               Exit status.
*/
int  hdfFree4DArray(float ****data);

