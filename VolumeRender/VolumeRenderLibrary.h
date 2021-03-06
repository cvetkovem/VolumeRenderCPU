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
  float *volume;
  uint32_t *xOffset;
  uint32_t *yOffset;
  uint32_t *zOffset;

  uint32_t xlen;
  uint32_t ylen;
  uint32_t zlen;
  float minDensity;
  float maxDensity;
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

  int32_t enableDrawBox;

  // M: model matrix
  float modelMatrix[16];
  // V: view matrix 
  float viewMatrix[16];
  // P: perspective matrix
  // float cameraPerspectiveMatrix[16];
  // VP: view perspective matrix 
  float viewPerspectiveMatrix[16];
  // MVP: model view perspective matrix
  float mvpMatrix[16];
  // MV: model view matrix (for fragment shader)
  float modelViewMatrix[16];
  
  // newCameraPosition (for fragment shader)
  float newCameraPosition[3];

  // const for all ray
  float rayStepSize;
  float boxXdiv2;
  float boxYdiv2;
  float boxZdiv2;

  void cameraPerspective(float angle, uint32_t width, uint32_t height, float zNear, float zFar);
  void cameraTranslate(float x, float y, float z);
  void cameraRotate(float *cameraTarget, float *cameraUp);
  void calcCameraMatrixs(uint32_t width, uint32_t height);

  void deleteVolume();
  void deleteInterpolate();

  // functions for render
  void renderCalculateMVPMatrix();
  float* renderCreateBoxForVolume();
  void renderDrawBox();
  void brezenhem(int x0, int y0, int x1, int y1, unsigned char *rgba);
  void renderPipeLine(float *pointsTriangle, float *depthBuffer);
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
  void setEnableDrawBox(int32_t en);
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
