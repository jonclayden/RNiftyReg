#include <RcppEigen.h>

#include "_reg_f3d.h"
#include "_reg_f3d_sym.h"
#include "_reg_localTransformation.h"

#include "config.h"
#include "aladin.h"
#include "f3d.h"
#include "AffineMatrix.h"
#include "DeformationField.h"

F3dResult regF3d (nifti_image *sourceImage, nifti_image *targetImage, const int nLevels, const int maxIterations, const int interpolation, nifti_image *sourceMaskImage, nifti_image *targetMaskImage, nifti_image *controlPointImage, AffineMatrix *initAffine, const int nBins, const std::vector<float> &spacing, const float bendingEnergyWeight, const float jacobianWeight, const float inverseConsistencyWeight, const bool symmetric, const bool verbose, const bool estimateOnly)
{
    if (controlPointImage == NULL && initAffine == NULL)
        initAffine = new AffineMatrix(sourceImage, targetImage);
    else if (controlPointImage != NULL)
        initAffine = NULL;
    
    // Binarise the mask images
    if (sourceMaskImage != NULL)
        reg_tools_binarise_image(sourceMaskImage);
    if (targetMaskImage != NULL)
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
        if (controlPointImage != NULL)
        {
            result.forwardControlPoints = copyCompleteImage(controlPointImage);
            DeformationField deformationField(targetImage, controlPointImage);
            result.forwardImage = deformationField.resampleImage(sourceImage, interpolation);
        }
        else
        {
            DeformationField deformationField(targetImage, *initAffine);
            result.forwardControlPoints = copyCompleteImage(deformationField.getFieldImage());
            result.forwardImage = deformationField.resampleImage(sourceImage, interpolation);
        }
        result.reverseImage = NULL;
        result.reverseControlPoints = NULL;
    }
    else
    {
        reg_f3d<PRECISION_TYPE> *reg = NULL;

        // Create the reg_f3d object
        if (symmetric)
            reg = new reg_f3d_sym<PRECISION_TYPE>(targetImage->nt, sourceImage->nt);
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
        
        if (targetMaskImage != NULL)
            reg->SetReferenceMask(targetMaskImage);
        
        if (controlPointImage != NULL)
            reg->SetControlPointGridImage(controlPointImage);
        
        mat44 affineMatrix;
        if (initAffine != NULL)
        {
            affineMatrix = *initAffine;
            reg->SetAffineTransformation(&affineMatrix);
        }
        
        reg->SetBendingEnergyWeight(bendingEnergyWeight);
        reg->SetLinearEnergyWeights(0.0, 0.0);
        reg->SetL2NormDisplacementWeight(0.0);
        reg->SetJacobianLogWeight(jacobianWeight);
        
        reg->SetMaximalIterationNumber(maxIterations);
        
        for (int i = 0; i < 3; i++)
            reg->SetSpacing((unsigned) i, (PRECISION_TYPE) spacing[i]);
        
        reg->SetLevelNumber(nLevels);
        reg->SetLevelToPerform(nLevels);

        if (interpolation == 3)
            reg->UseCubicSplineInterpolation();
        else if (interpolation == 1)
            reg->UseLinearInterpolation();
        else
            reg->UseNeareatNeighborInterpolation();
        
        // Parameters only relevant to the symmetric version of the algorithm
        if (symmetric)
        {
            if (sourceMaskImage != NULL)
                reg->SetFloatingMask(sourceMaskImage);
            
            reg->SetInverseConsistencyWeight(inverseConsistencyWeight);
        }
        
        // Run the registration
        reg->Run();
        
        result.forwardImage = NULL;
        result.reverseImage = NULL;
        
        if (!estimateOnly)
            result.forwardImage = reg->GetWarpedImage()[0];
        result.forwardControlPoints = reg->GetControlPointPositionImage();
        if (symmetric)
        {
            if (!estimateOnly)
                result.reverseImage = reg->GetWarpedImage()[1];
            result.reverseControlPoints = reg->GetBackwardControlPointPositionImage();
        }
        result.iterations = reg->GetCompletedIterations();
        
        // Erase the registration object
        delete reg;
    }
    
    return result;
}
