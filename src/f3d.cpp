#include <RcppEigen.h>

#include "_reg_f3d.h"
#include "_reg_f3d2.h"

#include "config.h"
#include "aladin.h"
#include "f3d.h"
#include "AffineMatrix.h"
#include "DeformationField.h"

F3dResult regF3d (const NiftiImage &sourceImage, const NiftiImage &targetImage, const int nLevels, const int maxIterations, const int interpolation, const NiftiImage &sourceMaskImage, const NiftiImage &targetMaskImage, const NiftiImage &initControlPoints, const AffineMatrix &initAffine, const int nBins, const std::vector<float> &spacing, const float bendingEnergyWeight, const float linearEnergyWeight, const float jacobianWeight, const bool symmetric, const bool verbose, const bool estimateOnly)
{
    // Binarise the mask images
    if (!sourceMaskImage.isNull())
        reg_tools_binarise_image(sourceMaskImage);
    if (!targetMaskImage.isNull())
        reg_tools_binarise_image(targetMaskImage);
    
    // Change data types for interpolation precision if necessary
    if (interpolation != 0)
    {
        reg_tools_changeDatatype<double>(sourceImage);
        if (symmetric)
            reg_tools_changeDatatype<double>(targetImage);
    }
    
    F3dResult result;
    
    if (nLevels == 0)
    {
        if (!initControlPoints.isNull())
        {
            result.forwardTransform = initControlPoints;
            DeformationField deformationField(targetImage, initControlPoints);
            result.image = deformationField.resampleImage(sourceImage, interpolation);
        }
        else
        {
            DeformationField deformationField(targetImage, initAffine);
            result.forwardTransform = deformationField.getFieldImage();
            result.image = deformationField.resampleImage(sourceImage, interpolation);
        }
    }
    else
    {
        reg_f3d<PRECISION_TYPE> *reg = NULL;

        // Create the reg_f3d object
        if (symmetric)
            reg = new reg_f3d2<PRECISION_TYPE>(targetImage->nt, sourceImage->nt);
        else
            reg = new reg_f3d<PRECISION_TYPE>(targetImage->nt, sourceImage->nt);
        
#ifdef _OPENMP
        const int maxThreadNumber = omp_get_max_threads();
        if (verbose)
            Rprintf("[NiftyReg F3D] Using OpenMP with %i thread(s)\n", maxThreadNumber);
#endif

        // Set the reg_f3d parameters
        reg->SetReferenceImage(targetImage);
        reg->SetFloatingImage(sourceImage);
        
        if (verbose)
            reg->PrintOutInformation();
        else
            reg->DoNotPrintOutInformation();
        
        if (!sourceMaskImage.isNull())
            reg->SetFloatingMask(sourceMaskImage);
        if (!targetMaskImage.isNull())
            reg->SetReferenceMask(targetMaskImage);
        
        mat44 affineMatrix;
        if (!initControlPoints.isNull())
            reg->SetControlPointGridImage(initControlPoints);
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
