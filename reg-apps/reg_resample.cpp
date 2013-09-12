/**
 * @file reg_resample.cpp
 * @author Marc Modat
 * @date 18/05/2009
 *
 *  Created by Marc Modat on 18/05/2009.
 *  Copyright (c) 2009, University College London. All rights reserved.
 *  Centre for Medical Image Computing (CMIC)
 *  See the LICENSE.txt file in the nifty_reg root folder
 *
 */

#ifndef _MM_RESAMPLE_CPP
#define _MM_RESAMPLE_CPP

#include <limits>

#include "_reg_ReadWriteImage.h"
#include "_reg_resampling.h"
#include "_reg_globalTransformation.h"
#include "_reg_localTransformation.h"
#include "_reg_tools.h"
#include "reg_resample.h"

typedef struct{
    char *referenceImageName;
    char *floatingImageName;
    char *affineMatrixName;
    char *inputCPPName;
    char *inputDEFName;
    char *outputResultName;
    char *outputBlankName;
    float sourceBGValue;
    int interpolation;
    float paddingValue;
}PARAM;
typedef struct{
    bool referenceImageFlag;
    bool floatingImageFlag;
    bool affineMatrixFlag;
    bool affineFlirtFlag;
    bool inputCPPFlag;
    bool inputDEFFlag;
    bool outputResultFlag;
    bool outputBlankFlag;
    bool outputBlankXYFlag;
    bool outputBlankYZFlag;
    bool outputBlankXZFlag;
}FLAG;


void PetitUsage(char *exec)
{
    fprintf(stderr,"Usage:\t%s -target <referenceImageName> -source <floatingImageName> [OPTIONS].\n",exec);
    fprintf(stderr,"\tSee the help for more details (-h).\n");
    return;
}
void Usage(char *exec)
{
    printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    printf("Usage:\t%s -ref <filename> -flo <filename> [OPTIONS].\n",exec);
    printf("\t-ref <filename>\tFilename of the reference image (mandatory)\n");
    printf("\t-flo <filename>\tFilename of the floating image (mandatory)\n\n");
#ifdef _SVN_REV
    fprintf(stderr,"\n-v Print the subversion revision number\n");
#endif

    printf("* * OPTIONS * *\n");
    printf("\t*\tOnly one of the following tranformation is taken into account\n");
    printf("\t-aff <filename>\t\tFilename which contains an affine transformation (Affine*Reference=floating)\n");
    printf("\t-affFlirt <filename>\t\tFilename which contains a radiological flirt affine transformation\n");
    printf("\t-cpp <filename>\t\tFilename of the control point grid image (from reg_f3d)\n");
    printf("\t-def <filename>\t\tFilename of the deformation field image (from reg_transform)\n");

    printf("\t*\tThere are no limit for the required output number from the following\n");
    printf("\t-res <filename> \tFilename of the resampled image [none]\n");
    printf("\t-blank <filename> \tFilename of the resampled blank grid [none]\n");

    printf("\t*\tOthers\n");
    printf("\t-NN \t\t\tUse a Nearest Neighbor interpolation for the source resampling (cubic spline by default)\n");
    printf("\t-LIN \t\t\tUse a Linear interpolation for the source resampling (cubic spline by default)\n");
    printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    return;
}

int main(int argc, char **argv)
{
    PARAM *param = (PARAM *)calloc(1,sizeof(PARAM));
    FLAG *flag = (FLAG *)calloc(1,sizeof(FLAG));

    param->interpolation=3; // Cubic spline interpolation used by default
    param->paddingValue=0;

    /* read the input parameter */
    for(int i=1;i<argc;i++){
        if(strcmp(argv[i], "-help")==0 || strcmp(argv[i], "-Help")==0 ||
           strcmp(argv[i], "-HELP")==0 || strcmp(argv[i], "-h")==0 ||
           strcmp(argv[i], "--h")==0 || strcmp(argv[i], "--help")==0){
            Usage(argv[0]);
            return 0;
        }
        else if(strcmp(argv[i], "--xml")==0){
            printf("%s",xml_resample);
            return 0;
        }
#ifdef _SVN_REV
        if( strcmp(argv[i], "-version")==0 ||
            strcmp(argv[i], "-Version")==0 ||
            strcmp(argv[i], "-V")==0 ||
            strcmp(argv[i], "-v")==0 ||
            strcmp(argv[i], "--v")==0 ||
            strcmp(argv[i], "--version")==0){
            printf("NiftyReg revision number: %i\n",_SVN_REV);
            return 0;
        }
#endif
        else if((strcmp(argv[i],"-ref")==0) || (strcmp(argv[i],"-target")==0) ||
                (strcmp(argv[i],"--ref")==0)){
            param->referenceImageName=argv[++i];
            flag->referenceImageFlag=1;
        }
        else if((strcmp(argv[i],"-flo")==0) || (strcmp(argv[i],"-source")==0) ||
                (strcmp(argv[i],"--flo")==0)){
            param->floatingImageName=argv[++i];
            flag->floatingImageFlag=1;
        }
        else if(strcmp(argv[i], "-aff") == 0 ||
                (strcmp(argv[i],"--aff")==0)){
            param->affineMatrixName=argv[++i];
            flag->affineMatrixFlag=1;
        }
        else if(strcmp(argv[i], "-affFlirt") == 0 ||
                (strcmp(argv[i],"--affFlirt")==0)){
            param->affineMatrixName=argv[++i];
            flag->affineMatrixFlag=1;
            flag->affineFlirtFlag=1;
        }
        else if((strcmp(argv[i],"-res")==0) || (strcmp(argv[i],"-result")==0) ||
                (strcmp(argv[i],"--res")==0)){
                param->outputResultName=argv[++i];
                flag->outputResultFlag=1;
        }
        else if(strcmp(argv[i], "-cpp") == 0 ||
                (strcmp(argv[i],"--cpp")==0)){
            param->inputCPPName=argv[++i];
            flag->inputCPPFlag=1;
        }
        else if(strcmp(argv[i], "-def") == 0 ||
                (strcmp(argv[i],"--def")==0)){
            param->inputDEFName=argv[++i];
            flag->inputDEFFlag=1;
        }
        else if(strcmp(argv[i], "-inter") == 0 ||
                (strcmp(argv[i],"--inter")==0)){
            param->interpolation=atoi(argv[++i]);
        }
        else if(strcmp(argv[i], "-pad") == 0 ||
                (strcmp(argv[i],"--pad")==0)){
            param->paddingValue=(float)atof(argv[++i]);
        }
        else if(strcmp(argv[i], "-NN") == 0){
            param->interpolation=0;
        }
        else if(strcmp(argv[i], "-LIN") == 0 ||
                (strcmp(argv[i],"-TRI")==0)){
            param->interpolation=1;
        }
        else if(strcmp(argv[i], "-CUB") == 0 ||
                (strcmp(argv[i],"-SPL")==0)){
            param->interpolation=3;
        }
        else if(strcmp(argv[i], "-blank") == 0 ||
                (strcmp(argv[i],"--blank")==0)){
            param->outputBlankName=argv[++i];
            flag->outputBlankFlag=1;
        }
        else if(strcmp(argv[i], "-blankXY") == 0 ||
                (strcmp(argv[i],"--blankXY")==0)){
            param->outputBlankName=argv[++i];
            flag->outputBlankXYFlag=1;
        }
        else if(strcmp(argv[i], "-blankYZ") == 0 ||
                (strcmp(argv[i],"--blankYZ")==0)){
            param->outputBlankName=argv[++i];
            flag->outputBlankYZFlag=1;
        }
        else if(strcmp(argv[i], "-blankXZ") == 0 ||
                (strcmp(argv[i],"--blankXZ")==0)){
            param->outputBlankName=argv[++i];
            flag->outputBlankXZFlag=1;
        }
        else{
            fprintf(stderr,"Err:\tParameter %s unknown.\n",argv[i]);
            PetitUsage(argv[0]);
            return 1;
        }
    }

    if(!flag->referenceImageFlag || !flag->floatingImageFlag){
    fprintf(stderr,"[NiftyReg ERROR] The reference and the floating image have both to be defined.\n");
            PetitUsage(argv[0]);
            return 1;
    }
	
    /* Check the number of input images */
    if(((unsigned int)flag->affineMatrixFlag +
        (unsigned int)flag->inputCPPFlag +
        (unsigned int)flag->inputDEFFlag) > 1){
        fprintf(stderr,"[NiftyReg ERROR] Only one input transformation has to be assigned.\n");
        PetitUsage(argv[0]);
        return 1;
    }

	/* Read the target image */
    nifti_image *referenceImage = reg_io_ReadImageHeader(param->referenceImageName);
    if(referenceImage == NULL){
        fprintf(stderr,"[NiftyReg ERROR] Error when reading the target image: %s\n",
                param->referenceImageName);
        return 1;
    }
    reg_checkAndCorrectDimension(referenceImage);
	
    /* Read the source image */
    nifti_image *floatingImage = reg_io_ReadImageFile(param->floatingImageName);
    if(floatingImage == NULL){
        fprintf(stderr,"[NiftyReg ERROR] Error when reading the source image: %s\n",
                param->floatingImageName);
        return 1;
    }
    reg_checkAndCorrectDimension(floatingImage);

    // Tell the CLI that the process has started
    startProgress("reg_resample");

    // Set up progress indicators
    float iProgressStep=1, nProgressSteps;

    /* *********************************** */
    /* DISPLAY THE RESAMPLING PARAMETERS */
    /* *********************************** */
    printf("\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
    printf("Command line:\n");
    for(int i=0;i<argc;i++) printf(" %s", argv[i]);
    printf("\n\n");
    printf("Parameters\n");
    printf("Target image name: %s\n",referenceImage->fname);
    printf("\t%ix%ix%i voxels, %i volumes\n",referenceImage->nx,referenceImage->ny,referenceImage->nz,referenceImage->nt);
    printf("\t%gx%gx%g mm\n",referenceImage->dx,referenceImage->dy,referenceImage->dz);
    printf("Source image name: %s\n",floatingImage->fname);
    printf("\t%ix%ix%i voxels, %i volumes\n",floatingImage->nx,floatingImage->ny,floatingImage->nz,floatingImage->nt);
    printf("\t%gx%gx%g mm\n",floatingImage->dx,floatingImage->dy,floatingImage->dz);
    printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n\n");

    /* *********************** */
    /* READ THE TRANSFORMATION */
    /* *********************** */
    nifti_image *controlPointImage = NULL;
    nifti_image *deformationFieldImage = NULL;
    int currentDatatype=NIFTI_TYPE_FLOAT32;
    int currentNbyper=sizeof(float);
    mat44 *affineTransformationMatrix = (mat44 *)calloc(1,sizeof(mat44));
    if(flag->inputCPPFlag){
#ifndef NDEBUG
        printf("[NiftyReg DEBUG] Name of the control point image: %s\n", param->inputCPPName);
#endif
        controlPointImage = reg_io_ReadImageFile(param->inputCPPName);
        if(controlPointImage == NULL){
            fprintf(stderr,"[NiftyReg ERROR] Error when reading the control point image: %s\n",param->inputCPPName);
            return 1;
        }
        reg_checkAndCorrectDimension(controlPointImage);
        currentDatatype=controlPointImage->datatype;
        currentNbyper=controlPointImage->nbyper;
    }
    else if(flag->inputDEFFlag){
#ifndef NDEBUG
        printf("[NiftyReg DEBUG] Name of the deformation field image: %s\n", param->inputDEFName);
#endif
        deformationFieldImage = reg_io_ReadImageFile(param->inputDEFName);
        if(deformationFieldImage == NULL){
            fprintf(stderr,"[NiftyReg ERROR] Error when reading the deformation field image: %s\n",param->inputDEFName);
            return 1;
        }
        reg_checkAndCorrectDimension(deformationFieldImage);
        currentDatatype=deformationFieldImage->datatype;
        currentNbyper=deformationFieldImage->nbyper;
    }
    else if(flag->affineMatrixFlag){
#ifndef NDEBUG
        printf("[NiftyReg DEBUG] Name of affine transformation: %s\n", param->affineMatrixName);
#endif
        // Check first if the specified affine file exist
        if(FILE *aff=fopen(param->affineMatrixName, "r")){
            fclose(aff);
        }
        else{
            fprintf(stderr,"The specified input affine file (%s) can not be read\n",param->affineMatrixName);
            return 1;
        }
        reg_tool_ReadAffineFile(affineTransformationMatrix,
                                referenceImage,
                                floatingImage,
                                param->affineMatrixName,
                                flag->affineFlirtFlag);
    }
    else{
        // identity transformation is considered
        affineTransformationMatrix->m[0][0]=1.0;
        affineTransformationMatrix->m[1][1]=1.0;
        affineTransformationMatrix->m[2][2]=1.0;
        affineTransformationMatrix->m[3][3]=1.0;
    }

    // Update progress via CLI
    progressXML(1, "Transform loaded...");

    // Allocate and compute the deformation field if necessary
    if(!flag->inputDEFFlag){
#ifndef NDEBUG
        printf("[NiftyReg DEBUG] Allocation of the deformation field\n");
#endif
        // Allocate
        deformationFieldImage = nifti_copy_nim_info(referenceImage);
        deformationFieldImage->dim[0]=deformationFieldImage->ndim=5;
        deformationFieldImage->dim[1]=deformationFieldImage->nx=referenceImage->nx;
        deformationFieldImage->dim[2]=deformationFieldImage->ny=referenceImage->ny;
        deformationFieldImage->dim[3]=deformationFieldImage->nz=referenceImage->nz;
        deformationFieldImage->dim[4]=deformationFieldImage->nt=1;deformationFieldImage->pixdim[4]=deformationFieldImage->dt=1.0;
        if(referenceImage->nz>1) deformationFieldImage->dim[5]=deformationFieldImage->nu=3;
        else deformationFieldImage->dim[5]=deformationFieldImage->nu=2;
        deformationFieldImage->pixdim[5]=deformationFieldImage->du=1.0;
        deformationFieldImage->dim[6]=deformationFieldImage->nv=1;deformationFieldImage->pixdim[6]=deformationFieldImage->dv=1.0;
        deformationFieldImage->dim[7]=deformationFieldImage->nw=1;deformationFieldImage->pixdim[7]=deformationFieldImage->dw=1.0;
        deformationFieldImage->nvox =(size_t)deformationFieldImage->nx*(size_t)deformationFieldImage->ny*(size_t)deformationFieldImage->nz*
                (size_t)deformationFieldImage->nt*(size_t)deformationFieldImage->nu;
        deformationFieldImage->datatype = currentDatatype;
        deformationFieldImage->nbyper = currentNbyper;
        deformationFieldImage->data = (void *)calloc(deformationFieldImage->nvox, deformationFieldImage->nbyper);
        //Computation
        if(flag->inputCPPFlag){
#ifndef NDEBUG
            printf("[NiftyReg DEBUG] Computation of the deformation field from the CPP image\n");
#endif
            if(controlPointImage->intent_p1==SPLINE_VEL_GRID){
                reg_spline_getDeformationFieldFromVelocityGrid(controlPointImage,
                                                               deformationFieldImage,
                                                               false // the number of step is not automatically updated
                                                               );
            }
            else{
                reg_tools_multiplyValueToImage(deformationFieldImage,deformationFieldImage,0.f);
                reg_getDeformationFromDisplacement(deformationFieldImage);
                reg_spline_getDeformationField(controlPointImage,
                                               deformationFieldImage,
                                               NULL, // mask
                                               true, //composition
                                               true // bspline
                                               );
            }
        }
        else{
#ifndef NDEBUG
            printf("[NiftyReg DEBUG] Computation of the deformation field from the affine transformation\n");
#endif
            reg_affine_positionField(affineTransformationMatrix,
                                     referenceImage,
                                     deformationFieldImage);
        }
    }

    // Update progress via CLI
    progressXML(2, "Deformation field ready...");

    /* ************************* */
    /* RESAMPLE THE SOURCE IMAGE */
    /* ************************* */
    if(flag->outputResultFlag){
        switch(param->interpolation){
        case 0:
            param->interpolation=0;
            break;
        case 1:
            param->interpolation=1;
            break;
        default:
            param->interpolation=3;
            break;
        }
        nifti_image *resultImage = nifti_copy_nim_info(referenceImage);
        resultImage->dim[0]=resultImage->ndim=floatingImage->dim[0];
        resultImage->dim[4]=resultImage->nt=floatingImage->dim[4];
        resultImage->cal_min=floatingImage->cal_min;
        resultImage->cal_max=floatingImage->cal_max;
        resultImage->scl_slope=floatingImage->scl_slope;
        resultImage->scl_inter=floatingImage->scl_inter;
        resultImage->datatype = floatingImage->datatype;
        resultImage->nbyper = floatingImage->nbyper;
        resultImage->nvox = (size_t)resultImage->dim[1] * (size_t)resultImage->dim[2] *
                (size_t)resultImage->dim[3] * (size_t)resultImage->dim[4];
        resultImage->data = (void *)calloc(resultImage->nvox, resultImage->nbyper);

        if(floatingImage->dim[4]==6 || floatingImage->dim[4]==7){
#ifndef _NDEBUG
            printf("[NiftyReg DEBUG] DTI-based resampling\n");
#endif
            // Compute first the Jacobian matrices
            mat33 *jacobian = (mat33 *)malloc(deformationFieldImage->nx *
                                              deformationFieldImage->ny *
                                              deformationFieldImage->nz *
                                              sizeof(mat33));
            reg_defField_getJacobianMatrix(deformationFieldImage,
                                           jacobian);
            // resample the DTI image
            bool timepoints[7]; for(int i=0;i<7;++i) timepoints[i]=true;
            if(floatingImage->dim[4]==7) timepoints[0]=false;
            reg_resampleImage(floatingImage,
                              resultImage,
                              deformationFieldImage,
                              NULL,
                              param->interpolation,
                              std::numeric_limits<float>::quiet_NaN(),
                              timepoints,
                              jacobian
                              );
        }
        else{
            reg_resampleImage(floatingImage,
                              resultImage,
                              deformationFieldImage,
                              NULL,
                              param->interpolation,
                              param->paddingValue);
        }
        memset(resultImage->descrip, 0, 80);
        strcpy (resultImage->descrip,"Warped image using NiftyReg (reg_resample)");
        reg_io_WriteImageFile(resultImage,param->outputResultName);
        printf("[NiftyReg] Resampled image has been saved: %s\n", param->outputResultName);
        nifti_image_free(resultImage);
    }

    /* *********************** */
    /* RESAMPLE A REGULAR GRID */
    /* *********************** */
    if(flag->outputBlankFlag ||
       flag->outputBlankXYFlag ||
       flag->outputBlankYZFlag ||
       flag->outputBlankXZFlag ){
        nifti_image *gridImage = nifti_copy_nim_info(floatingImage);
        gridImage->cal_min=0;
        gridImage->cal_max=255;
        gridImage->datatype = NIFTI_TYPE_UINT8;
        gridImage->nbyper = sizeof(unsigned char);
        gridImage->data = (void *)calloc(gridImage->nvox, gridImage->nbyper);
        unsigned char *gridImageValuePtr = static_cast<unsigned char *>(gridImage->data);
        for(int z=0; z<gridImage->nz;z++){
            for(int y=0; y<gridImage->ny;y++){
                for(int x=0; x<gridImage->nx;x++){
                    if(referenceImage->nz>1){
                        if(flag->outputBlankXYFlag){
                            if( x/10==(float)x/10.0 || y/10==(float)y/10.0)
                                *gridImageValuePtr = 255;
                        }
                        else if(flag->outputBlankYZFlag){
                            if( y/10==(float)y/10.0 || z/10==(float)z/10.0)
                                *gridImageValuePtr = 255;
                        }
                        else if(flag->outputBlankXZFlag){
                            if( x/10==(float)x/10.0 || z/10==(float)z/10.0)
                                *gridImageValuePtr = 255;
                        }
                        else{
                            if( x/10==(float)x/10.0 || y/10==(float)y/10.0 || z/10==(float)z/10.0)
                                *gridImageValuePtr = 255;
                        }
                    }
                    else{
                        if( x/10==(float)x/10.0 || x==referenceImage->nx-1 || y/10==(float)y/10.0 || y==referenceImage->ny-1)
                            *gridImageValuePtr = 255;
                    }
                    gridImageValuePtr++;
                }
            }
        }

        nifti_image *resultImage = nifti_copy_nim_info(referenceImage);
        resultImage->dim[0]=resultImage->ndim=3;
        resultImage->dim[4]=resultImage->nt=1;
        resultImage->cal_min=floatingImage->cal_min;
        resultImage->cal_max=floatingImage->cal_max;
        resultImage->scl_slope=floatingImage->scl_slope;
        resultImage->scl_inter=floatingImage->scl_inter;
        resultImage->datatype =NIFTI_TYPE_UINT8;
        resultImage->nbyper = sizeof(unsigned char);
        resultImage->data = (void *)calloc(resultImage->nvox, resultImage->nbyper);
        reg_resampleImage(gridImage,
                          resultImage,
                          deformationFieldImage,
                          NULL,
                          1, // linear interpolation
                          0);
        memset(resultImage->descrip, 0, 80);
        strcpy (resultImage->descrip,"Warped regular grid using NiftyReg (reg_resample)");
        reg_io_WriteImageFile(resultImage,param->outputBlankName);
        nifti_image_free(resultImage);
        nifti_image_free(gridImage);
        printf("[NiftyReg] Resampled grid has been saved: %s\n", param->outputBlankName);
    }

    // Tell the CLI that we finished
    closeProgress("reg_resample", "Normal exit");

    nifti_image_free(referenceImage);
    nifti_image_free(floatingImage);
    nifti_image_free(controlPointImage);
    nifti_image_free(deformationFieldImage);
    free(affineTransformationMatrix);

    free(flag);
    free(param);
    return 0;
}

#endif
