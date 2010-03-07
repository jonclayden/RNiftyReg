/*
 *  _reg_bspline.h
 *  
 *
 *  Created by Marc Modat on 25/03/2009.
 *  Copyright (c) 2009, University College London. All rights reserved.
 *  Centre for Medical Image Computing (CMIC)
 *  See the LICENSE.txt file in the nifty_reg root folder
 *
 */

#ifndef _REG_BSPLINE_H
#define _REG_BSPLINE_H

#include "nifti1_io.h"
#include "_reg_affineTransformation.h"

#if _USE_SSE
	#include <emmintrin.h>
#endif

extern "C++" template<class PrecisionTYPE>
void reg_bspline(	nifti_image *splineControlPoint,
			nifti_image *targetImage,
			nifti_image *deformationField,
            int *mask,
			int type
 			);

extern "C++" template<class PrecisionTYPE>
PrecisionTYPE reg_bspline_bendingEnergy(	nifti_image *splineControlPoint,
					nifti_image *targetImage,
					int type
					);
extern "C++" template<class PrecisionTYPE>
void reg_bspline_bendingEnergyGradient(	nifti_image *splineControlPoint,
									  nifti_image *targetImage,
									  nifti_image *gradientImage,
									  float weight
									  );
extern "C++" template<class PrecisionTYPE>
void reg_bspline_jacobianDeterminantGradient(	nifti_image *splineControlPoint,
												nifti_image *targetImage,
												nifti_image *gradientImage,
												float weight,
												bool approx
												);

extern "C++" template<class PrecisionTYPE>
PrecisionTYPE reg_bspline_jacobian(	nifti_image *splineControlPoint,
								  nifti_image *targetImage,
								  int type
								  );

extern "C++"
void reg_bspline_GetJacobianMap(	nifti_image *splineControlPoint,
							   nifti_image *jacobianImage
							   );

extern "C++"
void reg_bspline_GetJacobianMatrix(	nifti_image *splineControlPoint,
									nifti_image *jacobianImage
							   );

extern "C++"
void reg_voxelCentric2NodeCentric(	nifti_image *nodeImage,
									nifti_image *voxelImage,
									float weight
								  );

extern "C++"
void reg_bspline_refineControlPointGrid(	nifti_image *targetImage,
					nifti_image *splineControlPoint);

extern "C++"
int reg_bspline_initialiseControlPointGridWithAffine(	mat44 *affineTransformation,
							nifti_image *controlPointImage
						      );




extern "C++"
int reg_square_cpp(nifti_image *positionGridImage,
				   nifti_image *decomposedGridImage);

extern "C++" template<class PrecisionTYPE>
int reg_getDisplacementFromPosition(nifti_image *controlPointImage);

extern "C++" template<class PrecisionTYPE>
int reg_getPositionFromDisplacement(nifti_image *controlPointImage);

extern "C++"
int reg_spline_Interpolant2Interpolator(nifti_image *inputImage,
										nifti_image *outputImage);
extern "C++"
int reg_bspline_GetJacobianMapFromVelocityField(nifti_image* velocityFieldImage,
												nifti_image* jacobianImage);

#endif
