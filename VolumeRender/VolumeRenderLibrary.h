#ifndef Volume_Render_Library
#define Volume_Render_Library

#include <stdint.h>

#define BLOCK_SIZE_X 8 
#define BLOCK_SIZE_Y 8
#define BLOCK_SIZE_Z 16

enum {
  VRL_OK = 0,
  VRL_INVALID_VALUE ,
  VRL_NO_MEMORY,
  VRL_NO_EFFECT,
  VRL_ERROR
};

enum {
  VRL_IMAGE_NULL = 0,
  VRL_IMAGE_SET,
  VRL_IMAGE_ALLOC
};

enum VRL_IMAGE_QUALITY {
  VRL_IMG_QUA_LOW = 0,
  VRL_IMG_QUA_MEDIUM,
  VRL_IMG_QUA_BEST,
  VRL_IMG_QUA_BEST_EXPERIMENTAL,
  VRL_IMG_QUA_REALISTIC
};

enum ANTIALIASING_VALUE{
  VRL_ANTIALIASING_X1  = 1,
  VRL_ANTIALIASING_X4  = 4,
  VRL_ANTIALIASING_X16 = 16
};

typedef struct _vrlColor {
  float r;
  float g;
  float b;
  float a;

  float ambient;
  float diffuse;
  float specular;
  //only for material
  float emission;
  float shininess;
} vrlColor;

typedef struct _LUTPoint {
    float density;
    vrlColor sColor;
    void *next;
} LUTPoint;

class VRL {
private:
  int16_t *volume;
  uint32_t *xOffset;
  uint32_t *yOffset;
  uint32_t *zOffset;

  uint32_t xlen;
  uint32_t ylen;
  uint32_t zlen;
  int32_t minDensity;
  int32_t maxDensity;
  float volumeBoxMaxLen;

  LUTPoint *headLUTPoint;

  float *interpR;
  float *interpG;
  float *interpB;
  float *interpA;
  float *interpAmbient;
  float *interpDiffuse;
  float *interpSpecular;
  float *interpEmission;
  float *interpShininess;

  float lightDirection[4];
  vrlColor lightColor;

  unsigned char *privateRenderimage;
  unsigned char *image;
  uint32_t imageWidth;
  uint32_t imageHeight;
  int8_t imageStatus;

  uint32_t xClipMin;
  uint32_t xClipMax;
  uint32_t yClipMin;
  uint32_t yClipMax;
  uint32_t zClipMin;
  uint32_t zClipMax;

  float cameraPosition[3];
  float cameraTarget[3];
  float cameraUp[3];
  float zNear;

  float cameraTranslateMatrix[16];
  float cameraRotateMatrix[16];
  float cameraPerspectiveMatrix[16];

  float volumeTranslateMatrix[16];
  float volumeRotateMatrix[16];

  float renderStepSize;
  ANTIALIASING_VALUE antialiasingValue;
  uint32_t imageScale;

  void cameraPerspective(float angle, uint32_t width, uint32_t height, float zNear, float zFar);
  void cameraTranslate(float x, float y, float z);
  void cameraRotate(float *cameraTarget, float *cameraUp);
  void calcCameraMatrixs(uint32_t width, uint32_t height);

  void vecCross(float *a, float *b, float *out);
  void vecNormalize(float *in);

  void multMatrix4x4(float *lMatrix, float *rMatrix, float *out);
  void matrix4MultVector4(float *matrix, float *vector, float *outVector);

  void deleteVolume();
  void deleteInterpolate();

  // functions for render
  float* renderGetMVPMatrix();
  float* renderCreateBoxForVolume();
  void renderPipeLine(float *mvpMatrix, float *pointsTriangle, float *depthBuffer);
  int32_t renderFragmentShader(float *rayPosition, unsigned char *pixelColor);

  float getDensityFromVolume(float x, float y, float z);
public:
  VRL();
  ~VRL();

/****************************
 * Filling volume functions *
 ****************************/
  // xlen: slice width, ylen: slice height, zlen: number of slices
  int setVolume(int16_t *volume, uint32_t xlen, uint32_t ylen, uint32_t zlen, int32_t minDensity, int32_t maxDensity);

/***********************
 * Setup LUT functions *
 ***********************/
  int clearLUT();
  int addLUTPoint(float density, vrlColor *sColor);
  int removeLUTPoint(float density);
  int interpolateLUT();

/***************************
 *  Setup light functions  *
 * (only one light source) *
 ***************************/
  // Adds or replace the directional light source to the scene
  int setLight(float *direction, vrlColor *sColor);

/*************************************
 * Result image allocation functions *
 *     (only RGBA(8888) format)      *
 *************************************/
  // Sets the buffer for the result output image
  int setImage(unsigned char *image, uint32_t width, uint32_t height);
  // Allocates a buffer for the result output image
  int allocImage(uint32_t width, uint32_t height);
  unsigned char *getImage();

/*************             
 * Rendering *
 *************/
  // Sets the rendering clipping box
  int setClippingBox(uint32_t xClipMin, uint32_t xClipMax, uint32_t yClipMin, uint32_t yClipMax, uint32_t zClipMin, uint32_t zClipMax);
  // Renders the scene
  int render();

  int setCameraConfigure(float *position, float *target, float *up, float zNear);

/***********************************
 * Affine transformation functions *
 ***********************************/
  int setRotateVolume(float angle, float *axis);
  int addRotateVolume(float angle, float *axis);
  int resetRotationVolume();
  int setTranslateVolume(float x, float y, float z);
  int resetTranslateVolume();

/********************************
 * Quality/Perfomance functions *
 ********************************/
  // Enables / Disables Antialiasing
  int enableAntialiasing(ANTIALIASING_VALUE value);
  // Set size step for render
  int setImageQuality(VRL_IMAGE_QUALITY imgQua);
};
#endif
