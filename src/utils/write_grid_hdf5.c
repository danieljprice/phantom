/*
 * This subroutine performs the calls to the HDF5 library
 * to write gridded data
 *
 */
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
void write_grid_hdf5(char *filename, char *datasetname, float *data, int *nx, int *ny, int *nz, int *ierr)
   {
   hid_t     file_id, dataset_id, dataspace_id;
   herr_t    status;
   herr_t    HDF5_error = -1;

   /* Create new HDF5 file (use H5F_ACC_TRUNC: overwrite if already exists) */
   
   printf(" writing %s \n",filename);
   file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
   if (file_id == HDF5_error)
      { printf("ERROR opening %s for write \n",filename); *ierr = 1; return; }
   
   /* Create the dataspace */
   hsize_t dims[3];
   dims[0] = *nx;
   dims[1] = *ny;
   dims[2] = *nz;
   dataspace_id = H5Screate_simple(3,dims,NULL);
   if (dataspace_id == HDF5_error)
      { printf("ERROR creating data space \n"); *ierr = 2; return; }

   /* Create the dataset */
   dataset_id = H5Dcreate1(file_id, datasetname, H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT);
   if (dataset_id == HDF5_error) 
      { printf("ERROR creating data set %s \n",datasetname); *ierr = 3; return; }

   /* write dataset to file */
   status = H5Dwrite(dataset_id,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
   if (status == HDF5_error) { printf("ERROR writing to file \n"); *ierr = 4; }

   /* Close the dataset */
   status = H5Dclose(dataset_id);
   if (status == HDF5_error) { printf("ERROR closing dataset \n"); *ierr = 5; }

   /* Close the dataspace */
   status = H5Sclose(dataspace_id);
   if (status == HDF5_error) { printf("ERROR closing dataspace \n"); *ierr = 6; }
      
   /* Close the file */
   status = H5Fclose( file_id );
   if (status == HDF5_error) { printf("ERROR closing file \n"); *ierr = 7; }
   
   }
