/** @file reg_aladin.h
 * @date 20/06/2012
 * @author Marc Modat
 * @brief Header file that contains the string to be returned
 * for the slicer extension
 */

char xml_aladin[] =
"<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"
"<executable>\n"
"  <category>Registration</category>\n"
"  <title>RegAladin</title>\n"
"  <description><![CDATA[Module/executable for global registration (rigid and/or affine) based on a block-matching approach and a trimmed least squared optimisation.]]></description>\n"
"  <version>0.0.1</version>\n"
"  <documentation-url> TODO</documentation-url>\n"
"  <license>BSD</license>\n"
"  <contributor>Marc Modat, Pankaj Daga, David Cash (UCL)</contributor>\n"
"  <parameters advanced=\"false\">\n"
"    <label>Input images. Reference and floating images are mandatory</label>\n"
"    <description>Input images to perform the registration</description>\n"
"    <file fileExtensions=\".nii,.nii.gz,.nrrd,.png\">\n"
"      <name>referenceImageName</name>\n"
"      <description>Reference image filename (also called Target of Fixed)</description>\n"
"      <label>Reference image</label>\n"
"      <default>required</default>\n"
"      <longflag>ref</longflag>\n"
"      <channel>input</channel>\n"
"    </file>\n"
"    <file fileExtensions=\".nii,.nii.gz,.nrrd,.png\">\n"
"      <name>referenceMaskImageName</name>\n"
"      <description>Reference mask image filename</description>\n"
"      <label>Ref. mask</label>\n"
"      <default></default>\n"
"      <longflag>rmask</longflag>\n"
"      <channel>input</channel>\n"
"    </file>\n"
"    <float>\n"
"      <name>smoothReferenceWidth</name>\n"
"      <description>Standard deviation in mm (voxel if negative) of the Gaussian kernel used to smooth the reference image</description>\n"
"      <label>Ref .Smooth</label>\n"
"      <default>0</default>\n"
"      <longflag>smooR</longflag>\n"
"      <channel>input</channel>\n"
"    </float>\n"
"    <file fileExtensions=\".nii,.nii.gz,.nrrd,.png\">\n"
"      <name>floatingImageName</name>\n"
"      <description>Floating image filename (also called Source of moving)</description>\n"
"      <label>Floating image</label>\n"
"      <default>required</default>\n"
"      <longflag>flo</longflag>\n"
"      <channel>input</channel>\n"
"    </file>\n"
#ifdef _BUILD_NR_DEV
"    <file fileExtensions=\".nii,.nii.gz,.nrrd,.png\">\n"
"      <name>floatingMaskImageName</name>\n"
"      <description>Floating mask image filename</description>\n"
"      <label>Flo. mask</label>\n"
"      <default></default>\n"
"      <longflag>fmask</longflag>\n"
"      <channel>input</channel>\n"
"    </file>\n"
#endif // _BUILD_NR_DEV
"    <float>\n"
"      <name>smoothFloatingWidth</name>\n"
"      <description>Standard deviation in mm (voxel if negative) of the Gaussian kernel used to smooth the Floating image</description>\n"
"      <label>Flo. smooth</label>\n"
"      <default>0</default>\n"
"      <longflag>smooF</longflag>\n"
"      <channel>input</channel>\n"
"    </float>\n"
"  </parameters>\n"
"  <parameters advanced=\"false\">\n"
"    <label>Input affine parametrisation</label>\n"
"    <description>Optional input affine transformation</description>\n"
"    <file fileExtensions=\".txt\">\n"
"      <name>inputAffineName</name>\n"
"      <description>Affine registration matrix stored as a text file</description>\n"
"      <label>Input affine trans. from NiftyReg</label>\n"
"      <default></default>\n"
"      <longflag>inaff</longflag>\n"
"      <channel>input</channel>\n"
"    </file>\n"
"    <file fileExtensions=\".txt\">\n"
"      <name>inputFlirtAffineName</name>\n"
"      <description>Affine registration matrix from flirt (FSL) stored as a text file</description>\n"
"      <label>Input affine trans. from Flirt(FSL)</label>\n"
"      <default></default>\n"
"      <longflag>affFlirt</longflag>\n"
"      <channel>input</channel>\n"
"    </file>\n"
"  </parameters>\n"
"  <parameters advanced=\"false\">\n"
"    <label>Registration output</label>\n"
"    <description>Final affine trnansformation and warped image</description>\n"
"    <string>\n"
"      <name>outputAffineFileName</name>\n"
"      <description>Affine registration matrix output, saved as a text file</description>\n"
"      <label>Output affine filename</label>\n"
"      <default>outputAffineResult.txt</default>\n"
"      <longflag>aff</longflag>\n"
"      <channel>output</channel>\n"
"    </string>\n"
"    <string>\n"
"      <name>outputWarpedImageName</name>\n"
"      <description>Warped floating image</description>\n"
"      <label>Output warped image</label>\n"
"      <default>outputAffineResult.nii</default>\n"
"      <longflag>res</longflag>\n"
"      <channel>output</channel>\n"
"    </string>\n"
"  </parameters>\n"
"  <parameters advanced=\"true\">\n"
"    <label>Various optimisation parameters</label>\n"
"    <description>Various optimisation parameters such as the size of the pyramid or the number of level to use in the pyramidal approach.</description>\n"
"    <integer>\n"
"      <name>levelPyramidNumber</name>\n"
"      <description>Number of level to use to generate the pyramids for the coarse-to-fine approach</description>\n"
"      <label>Level number</label>\n"
"      <default>3</default>\n"
"      <longflag>ln</longflag>\n"
"      <channel>input</channel>\n"
"      <constraints>\n"
"        <minimum>0</minimum>\n"
"        <maximum>10</maximum>\n"
"      </constraints>\n"
"    </integer>\n"
"    <integer>\n"
"      <name>levelToPerformNumber</name>\n"
"      <description>Number of level to use to run the registration once the pyramids have been created</description>\n"
"      <label>Level to perform</label>\n"
"      <default>3</default>\n"
"      <longflag>lp</longflag>\n"
"      <channel>input</channel>\n"
"      <constraints>\n"
"        <minimum>0</minimum>\n"
"        <maximum>10</maximum>\n"
"      </constraints>\n"
"    </integer>\n"
"    <integer>\n"
"      <name>iterationNumber</name>\n"
"      <description>Maximal number of iteration of the trimmed least square approach to perform per total</description>\n"
"      <label>Iteration number</label>\n"
"      <default>5</default>\n"
"      <longflag>maxit</longflag>\n"
"      <channel>input</channel>\n"
"      <constraints>\n"
"        <minimum>1</minimum>\n"
"        <maximum>100</maximum>\n"
"      </constraints>\n"
"    </integer>\n"
"    <float>\n"
"      <name>blockPercentage</name>\n"
"      <description>Percentage of block to use in the optimisation scheme</description>\n"
"      <label>Percentage block</label>\n"
"      <default>50</default>\n"
"      <longflag>v</longflag>\n"
"      <channel>input</channel>\n"
"      <constraints>\n"
"        <minimum>1</minimum>\n"
"        <maximum>100</maximum>\n"
"      </constraints>\n"
"    </float>\n"
"    <float>\n"
"      <name>inlierPercentage</name>\n"
"      <description>Percentage of block to consider as inlier in the optimisation scheme</description>\n"
"      <label>Percentage inlier</label>\n"
"      <default>50</default>\n"
"      <longflag>i</longflag>\n"
"      <channel>input</channel>\n"
"      <constraints>\n"
"        <minimum>1</minimum>\n"
"        <maximum>100</maximum>\n"
"      </constraints>\n"
"    </float>\n"
#ifdef _BUILD_NR_DEV
"    <boolean>\n"
"      <name>useSym</name>\n"
"      <description>Performs a symmetric registration where both, forward and backward transformations are optimised</description>\n"
"      <label>Use symmetry</label>\n"
"      <default>false</default>\n"
"      <longflag>sym</longflag>\n"
"      <channel>input</channel>\n"
"    </boolean>\n"
#endif // _BUILD_NR_DEV
"    <boolean>\n"
"      <name>rigidOnly</name>\n"
"      <description>Performs only a rigid registration, rigid then affine by default</description>\n"
"      <label>Rigid only</label>\n"
"      <default>false</default>\n"
"      <longflag>rigOnly</longflag>\n"
"      <channel>input</channel>\n"
"    </boolean>\n"
"    <boolean>\n"
"      <name>affineOnly</name>\n"
"      <description>Performs only an affine registration, rigid then affine by default</description>\n"
"      <label>Affine only</label>\n"
"      <default>false</default>\n"
"      <longflag>affDirect</longflag>\n"
"      <channel>input</channel>\n"
"    </boolean>\n"
"    <boolean>\n"
"      <name>useHeaderOrigin</name>\n"
"      <description>Use the header origin to initialise the transformation. Image centres are used by default</description>\n"
"      <label>Use header</label>\n"
"      <default>false</default>\n"
"      <longflag>nac</longflag>\n"
"      <channel>input</channel>\n"
"    </boolean>\n"
"    <integer-enumeration>\n"
"      <name>interpolation</name>\n"
"      <description>Interpolation order to use internally to warp the floating image</description>\n"
"      <label>Interpolation order</label>\n"
"      <default>1</default>\n"
"      <longflag>interp</longflag>\n"
"      <channel>input</channel>\n"
"      <element>0</element>\n"
"      <element>1</element>\n"
"      <element>3</element>\n"
"    </integer-enumeration>\n"
"  </parameters>\n"
"</executable>\n"
;
