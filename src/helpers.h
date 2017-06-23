#ifndef _HELPERS_H_
#define _HELPERS_H_

#include "RNifti.h"

int nonunitaryDims (const RNifti::NiftiImage &image);

bool isMultichannel (const RNifti::NiftiImage &image);

RNifti::NiftiImage collapseChannels (const RNifti::NiftiImage &image);

void checkImages (const RNifti::NiftiImage &sourceImage, const RNifti::NiftiImage &targetImage);

RNifti::NiftiImage normaliseImage (const RNifti::NiftiImage &image);

RNifti::NiftiImage allocateMultiregResult (const RNifti::NiftiImage &source, const RNifti::NiftiImage &target, const bool forceDouble);

#endif
