#include <RcppEigen.h>

#include "_reg_f3d.h"
#include "_reg_f3d2.h"

#include "config.h"
#include "helpers.h"
#include "aladin.h"
#include "f3d.h"
#include "AffineMatrix.h"
#include "DeformationField.h"

using namespace RNifti;

F3dResult regF3d (const NiftiImage &sourceImage, const NiftiImage &targetImage, const int nLevels, const int maxIterations, const int interpolation, const NiftiImage &sourceMaskImage, const NiftiImage &targetMaskImage, const NiftiImage &initControlPoints, const AffineMatrix &initAffine, const int nBins, const std::vector<float> &spacing, const float bendingEnergyWeight, const float linearEnergyWeight, const float jacobianWeight, const bool symmetric, const bool verbose, const bool estimateOnly)
{
    NiftiImage source, target;
    source = normaliseImage(isMultichannel(sourceImage) ? collapseChannels(sourceImage) : sourceImage);
    target = normaliseImage(isMultichannel(targetImage) ? collapseChannels(targetImage) : targetImage);
    NiftiImage sourceMask(sourceMaskImage);
    NiftiImage targetMask(targetMaskImage);
    NiftiImage controlPoints(initControlPoints);
    
    // Binarise the mask images
    if (!sourceMask.isNull())
        reg_tools_binarise_image(sourceMask);
    if (!targetMask.isNull())
        reg_tools_binarise_image(targetMask);
    
    // Change data types for interpolation precision if necessary
    if (interpolation != 0)
    {
        reg_tools_changeDatatype<double>(source);
        if (symmetric)
            reg_tools_changeDatatype<double>(target);
    }
    
    F3dResult result;
    
    if (nLevels == 0)
    {
        if (!controlPoints.isNull())
        {
            result.forwardTransform = controlPoints;
            DeformationField deformationField(target, controlPoints);
            result.image = deformationField.resampleImage(source, interpolation);
        }
        else
        {
            DeformationField deformationField(target, initAffine);
            result.forwardTransform = deformationField.getFieldImage();
            result.image = deformationField.resampleImage(source, interpolation);
        }
    }
    else
    {
        reg_f3d<PRECISION_TYPE> *reg = NULL;

        // Create the reg_f3d object
        if (symmetric)
            reg = new reg_f3d2<PRECISION_TYPE>(target->nt, source->nt);
        else
            reg = new reg_f3d<PRECISION_TYPE>(target->nt, source->nt);
        
#ifdef _OPENMP
        const int maxThreadNumber = omp_get_max_threads();
        if (verbose)
            Rprintf("[NiftyReg F3D] Using OpenMP with %i thread(s)\n", maxThreadNumber);
#endif

        // Set the reg_f3d parameters
        reg->SetReferenceImage(target);
        reg->SetFloatingImage(source);
        
        if (verbose)
            reg->PrintOutInformation();
        else
            reg->DoNotPrintOutInformation();
        
        if (!sourceMask.isNull())
            reg->SetFloatingMask(sourceMask);
        if (!targetMask.isNull())
            reg->SetReferenceMask(targetMask);
        
        mat44 affineMatrix;
        if (!controlPoints.isNull())
            reg->SetControlPointGridImage(controlPoints);
        else
        {
            affineMatrix = initAffine;
            reg->SetAffineTransformation(&affineMatrix);
        }
        
        reg->SetBendingEnergyWeight(bendingEnergyWeight);
        reg->SetLinearEnergyWeight(linearEnergyWeight);
        reg->SetJacobianLogWeight(jacobianWeight);
        
        reg->SetMaximalIterationNumber(maxIterations);
        
        for (int i = 0; i < 3; i++)
            reg->SetSpacing(unsigned(i), PRECISION_TYPE(spacing[i]));
        
        for (int i = 0; i < target->nt; i++)
        {
            reg->UseNMISetReferenceBinNumber(i, nBins);
            reg->UseNMISetFloatingBinNumber(i, nBins);
        }
        
        reg->SetLevelNumber(nLevels);
        reg->SetLevelToPerform(nLevels);

        if (interpolation == 3)
            reg->UseCubicSplineInterpolation();
        else if (interpolation == 1)
            reg->UseLinearInterpolation();
        else
            reg->UseNeareatNeighborInterpolation();
        
        // Run the registration
        reg->Run();
        
        if (!estimateOnly)
            result.image = NiftiImage(reg->GetWarpedImage()[0]);
        result.forwardTransform = NiftiImage(reg->GetControlPointPositionImage());
        if (symmetric)
            result.reverseTransform = NiftiImage(reg->GetBackwardControlPointPositionImage());
        result.iterations = reg->GetCompletedIterations();
        
        // Erase the registration object
        delete reg;
    }
    
    return result;
}
