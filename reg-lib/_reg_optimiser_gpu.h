#ifndef _REG_OPTIMISER_GPU_H
#define _REG_OPTIMISER_GPU_H

#include "_reg_optimiser.h"
#include "_reg_common_gpu.h"
#include "_reg_tools_gpu.h"


/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
/** @class Global optimisation class for GPU
 * @brief Standard gradient acent optimisation
 */
class reg_optimiser_gpu : public reg_optimiser<float>
{
protected:
  float4 *currentDOF_gpu; // pointer to the cpp nifti image
  float4 *bestDOF_gpu;
  float4 *gradient_gpu;

public:
  reg_optimiser_gpu();
  ~reg_optimiser_gpu();

  virtual float* GetBestDOF(){return reinterpret_cast<float *>(this->bestDOF_gpu);}
  virtual float* GetCurrentDOF(){return reinterpret_cast<float *>(this->currentDOF_gpu);}
  virtual float* GetGradient(){return reinterpret_cast<float *>(this->gradient_gpu);}

  virtual void RestoreBestDOF();
  virtual void StoreCurrentDOF();

  virtual void Initialise(size_t nvox,
                          int dim,
                          bool optX,
                          bool optY,
                          bool optZ,
                          size_t maxit,
                          size_t start,
                          InterfaceOptimiser *o,
                          float *cppData,
                          float *gradData=NULL,
                          size_t a=0,
                          float *b=NULL,
                          float *c=NULL);
  virtual void Optimise(float maxLength,
                        float smallLength,
                        float &startLength);
  virtual void Perturbation(float length);
  virtual void NormaliseGradient();
};
/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
/** @class Conjugate gradient optimisation class for GPU
 * @brief
 */
class reg_conjugateGradient_gpu : public reg_optimiser_gpu
{
  protected:
    float4 *array1;
    float4 *array2;
    bool firstcall;
    void UpdateGradientValues(); /// @brief Update the gradient array

  public:
    reg_conjugateGradient_gpu();
    ~reg_conjugateGradient_gpu();

    virtual void Initialise(size_t nvox,
                            int dim,
                            bool optX,
                            bool optY,
                            bool optZ,
                            size_t maxit,
                            size_t start,
                            InterfaceOptimiser *o,
                            float *cppData,
                            float *gradData=NULL,
                            size_t a=0,
                            float *b=NULL,
                            float *c=NULL);
    virtual void Optimise(float maxLength,
                          float smallLength,
                          float &startLength);
    virtual void Perturbation(float length);
};
/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
/* \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ */
/** @brief
 */
extern "C++"
void reg_initialiseConjugateGradient_gpu(float4 **nodeNMIGradientArray_d,
                                         float4 **conjugateG_d,
                                         float4 **conjugateH_d,
                                         int nodeNumber);

/** @brief
 */
extern "C++"
void reg_GetConjugateGradient_gpu(float4 **nodeNMIGradientArray_d,
                                  float4 **conjugateG_d,
                                  float4 **conjugateH_d,
                                  int nodeNumber);

/** @brief
 */
extern "C++"
float reg_getMaximalLength_gpu(float4 **nodeNMIGradientArray_d,
                               int nodeNumber);

/** @brief
 */
extern "C++"
void reg_updateControlPointPosition_gpu(nifti_image *controlPointImage,
                                        float4 **controlPointImageArray_d,
                                        float4 **bestControlPointPosition_d,
                                        float4 **nodeNMIGradientArray_d,
                                        float currentLength);

#endif // _REG_OPTIMISER_GPU_H