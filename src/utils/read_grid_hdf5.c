/*
 * This subroutine reads HDF5 gridded data
 *
 */
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <hdf5.h>
void read_grid_hdf5_header(char *filename, int *nx, int *ny, int *nz, int *ncol, int *ierr)
   {
   hid_t     file_id;
   hid_t     dataset_id;
   hid_t     dataspace_id;
   herr_t    status;
   herr_t    HDF5_error = -1;

   printf(" Opening %s \n",filename);
   file_id = H5Fopen(filename,H5F_ACC_RDONLY,H5P_DEFAULT);
   if (file_id == HDF5_error)
      { printf("ERROR opening %s \n",filename); *ierr = 1; return; }

   // get number of datasets in file
   hsize_t ndatasets[1];
   H5Gget_num_objs(file_id, ndatasets);

   *ncol = ndatasets[0];

   // get names of all datasets in file
   // set function value to true (1) if it is present

   int i;
   char name[256];
   printf(" File contains %i HDF5 datasets... \n",(int) ndatasets[0]);
   for (i=0;i<ndatasets[0];i++) 
     {
       H5Gget_objname_by_idx(file_id, i, name, 256);
       printf(" Opening \"%s\" \n",name);

       dataset_id = H5Dopen1(file_id,name);
       if (dataset_id == HDF5_error) 
           { printf("ERROR opening %s data set \n",name); *ierr = 2; return; }

       dataspace_id = H5Dget_space(dataset_id);

       // get dimensional information from dataspace
       hsize_t HDFxdims[4], HDFmaxdims[4];
       int rank = H5Sget_simple_extent_dims(dataspace_id, HDFxdims, HDFmaxdims);
       if (rank > 4) { printf("RANK of dataset exceeds array bounds \n"); *ierr = 3; return; }
       if (rank != 3) { printf("RANK of dataset not equal to 3 \n"); *ierr = 4; }

       // from the dimensional info, calculate the size of the buffer.
       *nx = HDFxdims[0];
       *ny = HDFxdims[1];
       *nz = HDFxdims[2];

       status = H5Sclose(dataspace_id);
       if (status == HDF5_error) { printf("ERROR closing dataspace \n"); *ierr = 5; }

       status = H5Dclose(dataset_id);
       if (status == HDF5_error) { printf("ERROR closing dataset \n"); *ierr = 6; }
     }
   
   status = H5Fclose( file_id );
   if (status == HDF5_error) { printf("ERROR closing file \n"); *ierr = 7; }
   
   }

void read_grid_hdf5_column(char *filename, int *ireadcol, int *nx, int *ny, int *nz, double *datgrid, int *ierr)
   {
   hid_t     file_id;
   hid_t     dataset_id, SPHdataset_id;
   hid_t     dataspace_id, SPHdataspace_id;
   herr_t    status;
   herr_t    HDF5_error = -1;

   file_id = H5Fopen(filename,H5F_ACC_RDONLY,H5P_DEFAULT);
   if (file_id == HDF5_error)
      { printf("ERROR re-opening %s \n",filename); *ierr = 1; return; }
   
   // get number of datasets in file
   hsize_t ndatasets[1];
   H5Gget_num_objs(file_id, ndatasets);

   int i;
   char name[256];
   printf(" File contains %i HDF5 datasets... \n",(int) ndatasets[0]);

   H5Gget_objname_by_idx(file_id, *ireadcol-1, name, 256);
   printf(" Opening \"%s\" \n",name);

   dataset_id = H5Dopen1(file_id,name);
   if (dataset_id == HDF5_error) 
       { printf("ERROR opening %s data set \n",name); *ierr = 2; return; }

   dataspace_id = H5Dget_space(dataset_id);

   status = H5Dread(dataset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,datgrid);
   if (status == HDF5_error) { printf("ERROR reading data from %s \n",name); *ierr = 3; }

   status = H5Sclose(dataspace_id);
   if (status == HDF5_error) { printf("ERROR closing dataspace \n"); *ierr = 4; }

   status = H5Dclose(dataset_id);
   if (status == HDF5_error) { printf("ERROR closing dataset \n"); *ierr = 4; }
   
   status = H5Fclose( file_id );
   if (status == HDF5_error) { printf("ERROR closing file \n"); *ierr = 7; }

   }

void read_grid_hdf5_column_float(char *filename, int *ireadcol, int *nx, int *ny, int *nz, float *datgrid, int *ierr)
   {
   hid_t     file_id;
   hid_t     dataset_id, SPHdataset_id;
   hid_t     dataspace_id, SPHdataspace_id;
   herr_t    status;
   herr_t    HDF5_error = -1;

   file_id = H5Fopen(filename,H5F_ACC_RDONLY,H5P_DEFAULT);
   if (file_id == HDF5_error)
      { printf("ERROR re-opening %s \n",filename); *ierr = 1; return; }
   
   // get number of datasets in file
   hsize_t ndatasets[1];
   H5Gget_num_objs(file_id, ndatasets);

   int i;
   char name[256];
   printf(" File contains %i HDF5 datasets... \n",(int) ndatasets[0]);

   H5Gget_objname_by_idx(file_id, *ireadcol-1, name, 256);
   printf(" Opening \"%s\" \n",name);

   dataset_id = H5Dopen1(file_id,name);
   if (dataset_id == HDF5_error) 
       { printf("ERROR opening %s data set \n",name); *ierr = 2; return; }

   dataspace_id = H5Dget_space(dataset_id);

   status = H5Dread(dataset_id,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,datgrid);
   if (status == HDF5_error) { printf("ERROR reading data from %s \n",name); *ierr = 3; }

   status = H5Sclose(dataspace_id);
   if (status == HDF5_error) { printf("ERROR closing dataspace \n"); *ierr = 4; }

   status = H5Dclose(dataset_id);
   if (status == HDF5_error) { printf("ERROR closing dataset \n"); *ierr = 4; }
   
   status = H5Fclose( file_id );
   if (status == HDF5_error) { printf("ERROR closing file \n"); *ierr = 7; }

   }


