#include <RcppEigen.h>

#include "_reg_aladin.h"
#include "_reg_aladin_sym.h"

#include "config.h"
#include "aladin.h"
#include "AffineMatrix.h"
#include "DeformationField.h"

// Run the "aladin" registration algorithm
AladinResult regAladin (const NiftiImage &sourceImage, const NiftiImage &targetImage, const LinearTransformScope scope, const bool symmetric, const int nLevels, const int maxIterations, const int useBlockPercentage, const int interpolation, const NiftiImage &sourceMaskImage, const NiftiImage &targetMaskImage, const AffineMatrix &initAffine, const bool verbose, const bool estimateOnly)
{
    // Binarise the mask images
    if (!sourceMaskImage.isNull())
        reg_tools_binarise_image(sourceMaskImage);
    if (!targetMaskImage.isNull())
        reg_tools_binarise_image(targetMaskImage);
    
    // The source data type is changed for interpolation precision if necessary
    if (interpolation != 0)
        reg_tools_changeDatatype<double>(sourceImage);
    
    AladinResult result;
    
    if (nLevels == 0)
    {
        DeformationField deformationField(targetImage, initAffine);
        result.image = deformationField.resampleImage(sourceImage, interpolation);
        result.affine = initAffine;
    }
    else
    {
        reg_aladin<PRECISION_TYPE> *reg;
        if (symmetric)
            reg = new reg_aladin_sym<PRECISION_TYPE>;
        else
            reg = new reg_aladin<PRECISION_TYPE>;
    
        reg->SetMaxIterations(maxIterations);
        reg->SetNumberOfLevels(nLevels);
        reg->SetLevelsToPerform(nLevels);
        reg->SetReferenceSigma(0.0);
        reg->SetFloatingSigma(0.0);
        reg->SetAlignCentre(1);
        reg->SetPerformAffine(scope == AffineScope);
        reg->SetPerformRigid(1);
        reg->SetVerbose(int(verbose));
        reg->SetBlockStepSize(1);
        reg->SetBlockPercentage(useBlockPercentage);
        reg->SetInlierLts(50.0);
        reg->SetInterpolation(interpolation);
        reg->setPlatformCode(NR_PLATFORM_CPU);
        reg->setCaptureRangeVox(3);
        
        // Set the reference and floating images
        reg->SetInputReference(targetImage);
        reg->SetInputFloating(sourceImage);
    
        // Set the initial affine transformation
        mat44 affineMatrix = initAffine;
        reg->SetTransformationMatrix(&affineMatrix);
    
        // Set the masks if defined
        if (sourceMaskImage != NULL)
            reg->SetInputFloatingMask(sourceMaskImage);
        if (targetMaskImage != NULL)
            reg->SetInputMask(targetMaskImage);
    
        // Run the registration
        reg->Run();
    
        // Store the results
        if (!estimateOnly)
            result.image = NiftiImage(reg->GetFinalWarpedImage());
        result.affine = AffineMatrix(*reg->GetTransformationMatrix());
        result.iterations = reg->GetCompletedIterations();
    
        delete reg;
    }
    
    return result;
}
