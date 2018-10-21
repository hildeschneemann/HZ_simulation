// initialize an HDF5 dataset
// function to write a frequency array to the dataset

#include "fisher.h"
#include <iostream>
#include <fstream>
#include "H5Cpp.h"
using namespace std;

void createHDF5(const H5std_string &FILE_NAME,const H5std_string &DSET_FREQ_NAME,
    const int &RANK_FREQ, const int &DIM0, const int &DIM1, const int &DIM2,
    const H5std_string &DSET_W_NAME, const int &RANK_W, const int &DIM_W0, const int &DIM_W1,
    const H5std_string &DSET_HI_NAME, const int &DIM_HI, const H5std_string &DSET_GEN_NAME)
{
    H5::H5File file(FILE_NAME, H5F_ACC_TRUNC);

	hsize_t dims_freq[3];
    hsize_t chunk_dims[3] = {1, DIM1, DIM2};	// chunk dimensions for freq
    hsize_t dims_w[2];
	hsize_t dims_hi[1];
    hsize_t dims_gen[1];
	dims_freq[0] = DIM0;
	dims_freq[1] = DIM1;
    dims_freq[2] = DIM2;
    dims_w[0] = DIM_W0;
	dims_w[1] = DIM_W1;
	dims_hi[0] = DIM_HI;
    dims_gen[0] = DIM0;

	H5::DataSpace dspace_freq = H5::DataSpace(RANK_FREQ, dims_freq);
    H5::DataSpace dspace_w = H5::DataSpace(RANK_W, dims_w);
	H5::DataSpace dspace_hi = H5::DataSpace(1, dims_hi);
    H5::DataSpace dspace_gen = H5::DataSpace(1, dims_gen);

    // Modify freq dataset creation property to enable chunking
    H5::DSetCreatPropList  *plist = new  H5::DSetCreatPropList;
    plist->setChunk(RANK_FREQ, chunk_dims);
    plist->setDeflate(6);

	H5::DataSet dset_freq(file.createDataSet(DSET_FREQ_NAME,
	                                 H5::PredType::NATIVE_DOUBLE, dspace_freq,
                                     *plist));
    H5::DataSet dset_w(file.createDataSet(DSET_W_NAME,
                                    H5::PredType::NATIVE_DOUBLE, dspace_w));
	H5::DataSet dset_hi(file.createDataSet(DSET_HI_NAME,
					H5::PredType::NATIVE_DOUBLE, dspace_hi));
    H5::DataSet dset_gen(file.createDataSet(DSET_GEN_NAME,
                                    H5::PredType::NATIVE_INT, dspace_gen));

    dspace_freq.close();
    delete plist;
    dspace_w.close();
	dspace_hi.close();
    dset_freq.close();
    dset_w.close();
	dset_hi.close();
    dset_gen.close();
    file.close();
}

void writeTimeStepHDF5(H5::H5File &file, H5::DataSet &dataset,
    const int &RANK, hsize_t dimsm[],
    double sdata[], int &indexGen)
{
    int i;
    H5::DataSpace dataspace;

    hsize_t* offset = new hsize_t[RANK];
	hsize_t* count = new hsize_t[RANK];
	hsize_t* stride =  new hsize_t[RANK];
	hsize_t* block = new hsize_t[RANK];

    for(i = 0; i < RANK; i++)
    {
        offset[i] = 0;
        count[i] = dimsm[i];
        stride[i] = 1;
        block[i] = 1;
    }
    offset[0] = indexGen;

    H5::DataSpace memspace(RANK, dimsm, NULL);

    dataspace = dataset.getSpace();
    dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);

    dataset.write(sdata, H5::PredType::NATIVE_DOUBLE, memspace, dataspace);

    dataspace.close();
    memspace.close();
}

void writeGenSaved(H5::H5File &file, H5::DataSet &dataset, int data[])
{
    H5::DataSpace dataspace;
    dataset.write(data, H5::PredType::NATIVE_INT);
}

void writeHISaved(H5::H5File &file, H5::DataSet &dataset, double data[])
{
    H5::DataSpace dataspace;
    dataset.write(data, H5::PredType::NATIVE_DOUBLE);
}

void writeAttributeHDF5(H5::H5File &file, H5::DataSet &dataset,
    const char *ATTR_NAME, int attr_data)
{
    hsize_t dims[1] = {1};
	H5::DataSpace attr_dataspace = H5::DataSpace(1, dims);

	H5::Attribute attribute = dataset.createAttribute(ATTR_NAME, H5::PredType::NATIVE_INT,
	                                          attr_dataspace);

	attribute.write(H5::PredType::NATIVE_INT, &attr_data);

    attr_dataspace.close();
    attribute.close();
}

void writeAttributeHDF5(H5::H5File &file, H5::DataSet &dataset,
    const char *ATTR_NAME, double attr_data)
{
    hsize_t dims[1] = {1};

	H5::DataSpace attr_dataspace = H5::DataSpace(1, dims);

	H5::Attribute attribute = dataset.createAttribute(ATTR_NAME, H5::PredType::NATIVE_DOUBLE,
	                                          attr_dataspace);

	attribute.write(H5::PredType::NATIVE_DOUBLE, &attr_data);

    attr_dataspace.close();
    attribute.close();
}
