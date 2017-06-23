#ifndef _HELPERS_H_
#define _HELPERS_H_

#include "RNifti.h"

bool isMultichannel (const RNifti::NiftiImage &image);

RNifti::NiftiImage collapseChannels (const RNifti::NiftiImage &image);

void checkImages (RNifti::NiftiImage &sourceImage, RNifti::NiftiImage &targetImage);

RNifti::NiftiImage allocateMultiregResult (const RNifti::NiftiImage &source, const RNifti::NiftiImage &target, const bool forceDouble);

#endif
