#include <RcppEigen.h>

#include "nifti1_io.h"

#include "nifti_image.h"

using namespace Rcpp;

// Convert an S4 "nifti" object, as defined in the oro.nifti package, to a "nifti_image" struct
nifti_image * retrieveImageFromNiftiS4 (const RObject &object, const bool copyData = true)
{
    nifti_1_header header;
    header.sizeof_hdr = 348;
    
    const std::vector<short> dims = object.slot("dim_");
    for (int i=0; i<8; i++)
        header.dim[i] = dims[i];
    
    header.intent_p1 = object.slot("intent_p1");
    header.intent_p2 = object.slot("intent_p2");
    header.intent_p3 = object.slot("intent_p3");
    header.intent_code = object.slot("intent_code");
    
    header.datatype = object.slot("datatype");
    header.bitpix = object.slot("bitpix");
    
    header.slice_start = object.slot("slice_start");
    header.slice_end = object.slot("slice_end");
    header.slice_code = object.slot("slice_code");
    header.slice_duration = object.slot("slice_duration");
    
    const std::vector<float> pixdims = object.slot("pixdim");
    for (int i=0; i<8; i++)
        header.pixdim[i] = pixdims[i];
    header.xyzt_units = object.slot("xyzt_units");
    
    header.vox_offset = object.slot("vox_offset");
    
    header.scl_slope = object.slot("scl_slope");
    header.scl_inter = object.slot("scl_inter");
    header.toffset = object.slot("toffset");
    
    header.cal_max = object.slot("cal_max");
    header.cal_min = object.slot("cal_min");
    header.glmax = header.glmin = 0;
    
    strncpy(header.descrip, as<std::string>(object.slot("descrip")).c_str(), 79);
    header.descrip[79] = '\0';
    strncpy(header.aux_file, as<std::string>(object.slot("aux_file")).c_str(), 23);
    header.aux_file[23] = '\0';
    strncpy(header.intent_name, as<std::string>(object.slot("intent_name")).c_str(), 15);
    header.intent_name[15] = '\0';
    strncpy(header.magic, as<std::string>(object.slot("magic")).c_str(), 3);
    header.magic[3] = '\0';
    
    header.qform_code = object.slot("qform_code");
    header.sform_code = object.slot("sform_code");
    
    header.quatern_b = object.slot("quatern_b");
    header.quatern_c = object.slot("quatern_c");
    header.quatern_d = object.slot("quatern_d");
    header.qoffset_x = object.slot("qoffset_x");
    header.qoffset_y = object.slot("qoffset_y");
    header.qoffset_z = object.slot("qoffset_z");
    
    const std::vector<float> srow_x = object.slot("srow_x");
    const std::vector<float> srow_y = object.slot("srow_y");
    const std::vector<float> srow_z = object.slot("srow_z");
    for (int i=0; i<4; i++)
    {
        header.srow_x[i] = srow_x[i];
        header.srow_y[i] = srow_y[i];
        header.srow_z[i] = srow_z[i];
    }
    
    if (header.datatype == DT_UINT8 || header.datatype == DT_INT16 || header.datatype == DT_INT32 || header.datatype == DT_INT8 || header.datatype == DT_UINT16 || header.datatype == DT_UINT32)
        header.datatype = DT_INT32;
    else if (header.datatype == DT_FLOAT32 || header.datatype == DT_FLOAT64)
        header.datatype = DT_FLOAT64;  // This assumes that sizeof(double) == 8
    else
        throw std::runtime_error("Data type is not supported");
    
    nifti_image *image = nifti_convert_nhdr2nim(header, NULL);
    
    const SEXP data = object.slot(".Data");
    if (!copyData || Rf_length(data) == 1)
        image->data = NULL;
    else
    {
        const size_t dataSize = nifti_get_volsize(image);
        image->data = calloc(1, dataSize);
        if (header.datatype == DT_INT32)
            memcpy(image->data, INTEGER(data), dataSize);
        else
            memcpy(image->data, REAL(data), dataSize);
    }
    
    return image;
}

nifti_image * retrieveImageFromArray (const RObject &object)
{
    int dims[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
    const std::vector<int> dimVector = object.attr("dim");
    
    const int nDims = std::min(7, int(dimVector.size()));
    dims[0] = nDims;
    for (int i=0; i<nDims; i++)
        dims[i+1] = dimVector[i];
    
    short datatype = DT_UNKNOWN;
    const int sexpType = object.sexp_type();
    if (sexpType == INTSXP || sexpType == LGLSXP)
        datatype = DT_INT32;
    else if (sexpType == REALSXP)
        datatype = DT_FLOAT64;
    
    nifti_image *image = nifti_make_new_nim(dims, datatype, TRUE);
    
    const size_t dataSize = nifti_get_volsize(image);
    if (datatype == DT_INT32)
        memcpy(image->data, INTEGER(object), dataSize);
    else
        memcpy(image->data, REAL(object), dataSize);
    
    if (object.hasAttribute("pixdim"))
    {
        const std::vector<float> pixdimVector = object.attr("dim");
        const int pixdimLength = pixdimVector.size();
        for (int i=0; i<std::min(pixdimLength,nDims); i++)
            image->pixdim[i+1] = pixdimVector[i];
    }
    
    return image;
}

nifti_image * retrieveImage (const SEXP _image, const bool readData = true)
{
    nifti_image *image = NULL;
    if (Rf_isString(_image))
    {
        std::string path = as<std::string>(_image);
        image = nifti_image_read(path.c_str(), int(readData));
    }
    else
    {
        RObject imageObject(_image);
        if (imageObject.hasAttribute(".nifti_image_ptr"))
        {
            XPtr<nifti_image> imagePtr(SEXP(imageObject.attr(".nifti_image_ptr")));
            image = imagePtr;
        }
        else if (imageObject.inherits("nifti"))
            image = retrieveImageFromNiftiS4(imageObject, readData);
        else if (imageObject.hasAttribute("dim"))
            image = retrieveImageFromArray(imageObject);
        else
            throw std::runtime_error("Cannot convert object of class \"" + as<std::string>(imageObject.attr("class")) + "\" to a nifti_image");
    }
    
    return image;
}
