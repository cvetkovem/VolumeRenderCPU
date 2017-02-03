#include "VolumeRenderLibrary.h"
#include "utils\VectorMath.h"

#include <stdio.h>

#define EPSILON 0.00001f

VRL::VRL() {
  this->volume = nullptr;
  this->xOffset = nullptr;
  this->yOffset = nullptr;
  this->zOffset = nullptr;
  this->volumeBoxMaxLen = 1.0f;

  this->headLUTPoint = nullptr;

  this->interpR = nullptr;
  this->interpG = nullptr;
  this->interpB = nullptr;
  this->interpA = nullptr;
  this->interpAmbient = nullptr;
  this->interpDiffuse = nullptr;
  this->interpSpecular = nullptr;
  this->interpEmission = nullptr;
  this->interpShininess = nullptr;

  this->privateRenderimage = nullptr;
  this->image = nullptr;
  this->imageStatus = VRL_IMAGE_NULL;

  this->renderStepSize = 1.0f;
  this->antialiasingValue = VRL_ANTIALIASING_X1;
  this->imageScale = 1;

  this->enableDrawBox = 0;
}

VRL::~VRL() {
  deleteVolume();
  clearLUT();

  if (this->image && this->imageStatus == VRL_IMAGE_ALLOC) {
    delete[] (this->image);
  }
}

void VRL::deleteVolume() {
  if (this->volume) {
    delete[](this->volume);
    this->volume = nullptr;
  }

  if (this->xOffset) {
    delete[](this->xOffset);
    this->xOffset = nullptr;
  }

  if (this->yOffset) {
    delete[](this->yOffset);
    this->yOffset = nullptr;
  }

  if (this->zOffset) {
    delete[](this->zOffset);
    this->zOffset = nullptr;
  }
}

void VRL::deleteInterpolate() {
  if (this->interpR) {
    delete[](this->interpR);
  }
  this->interpR = nullptr;
  this->interpG = nullptr;
  this->interpB = nullptr;
  this->interpA = nullptr;
  this->interpAmbient = nullptr;
  this->interpDiffuse = nullptr;
  this->interpSpecular = nullptr;
  this->interpEmission = nullptr;
  this->interpShininess = nullptr;
}

int VRL::setVolume(int16_t *_volume, uint32_t _xlen, uint32_t _ylen, uint32_t _zlen, int32_t _minDensity, int32_t _maxDensity) {
  uint32_t numberBlock;
  uint32_t index;
  uint32_t xlenNew = _xlen;
  uint32_t ylenNew = _ylen;
  uint32_t zlenNew = _zlen;

  this->xlen = _xlen;
  this->ylen = _ylen;
  this->zlen = _zlen;
  this->minDensity = _minDensity;
  this->maxDensity = _maxDensity;

  // delete old data
  deleteVolume();
  
  if (_xlen % BLOCK_SIZE_X) {
    xlenNew = ((_xlen / BLOCK_SIZE_X) + 1) * BLOCK_SIZE_X;
  }
  if (_ylen % BLOCK_SIZE_Y) {
    ylenNew = ((_ylen / BLOCK_SIZE_Y) + 1) * BLOCK_SIZE_Y;
  }
  if (_zlen % BLOCK_SIZE_Z) {
    zlenNew = ((_zlen / BLOCK_SIZE_Z) + 1) * BLOCK_SIZE_Z;
  }

  this->volume = new int16_t[xlenNew * ylenNew * zlenNew];
  this->xOffset = new uint32_t[xlenNew];
  this->yOffset = new uint32_t[ylenNew];
  this->zOffset = new uint32_t[zlenNew];

  // calc xOffset
  numberBlock = 0;
  index = 0;
  for (uint32_t i = 0; i < _xlen; i++) {
    xOffset[i] = numberBlock * (BLOCK_SIZE_X * BLOCK_SIZE_Y * BLOCK_SIZE_Z) + index;
    index++;
    if (index == BLOCK_SIZE_X) {
      index = 0;
      numberBlock++;
    }
  }
  // calc yOffset
  numberBlock = 0;
  index = 0;
  for (uint32_t i = 0; i < _ylen; i++) {
    yOffset[i] = numberBlock * (BLOCK_SIZE_Y * BLOCK_SIZE_Z * xlenNew) + index * BLOCK_SIZE_X;
    index++;
    if (index == BLOCK_SIZE_Y) {
      index = 0;
      numberBlock++;
    }
  }
  // calc zOffset
  numberBlock = 0;
  index = 0;
  for (uint32_t i = 0; i < _zlen; i++) {
    zOffset[i] = numberBlock * (BLOCK_SIZE_Z * xlenNew * ylenNew) + index * (BLOCK_SIZE_X * BLOCK_SIZE_Y);
    index++;
    if (index == BLOCK_SIZE_Z) {
      index = 0;
      numberBlock++;
    }
  }

  int16_t tmp;
  for (uint32_t k = 0; k < _zlen; k++) {
    for (uint32_t j = 0; j < _ylen; j++) {
      for (uint32_t i = 0; i < _xlen; i++) {
        tmp = _volume[k*(_xlen*_ylen) + j*_xlen + i];
        
        tmp = (tmp < _minDensity) ? (_minDensity) : (tmp);
        tmp = (tmp > _maxDensity) ? (_maxDensity) : (tmp);

        this->volume[xOffset[i] + yOffset[j] + zOffset[k]] = tmp;
      }
    }
  }

  return VRL_OK;
}

int VRL::clearLUT() {
  LUTPoint *next = this->headLUTPoint;
  
  while (this->headLUTPoint) {
      next = (LUTPoint *)(this->headLUTPoint->next);
      delete (this->headLUTPoint);
      this->headLUTPoint = next;
  }

  deleteInterpolate();

  return VRL_OK;
}

int VRL::addLUTPoint(float density, vrlColor *sColor) {
  LUTPoint *current = this->headLUTPoint;
  LUTPoint *previous = nullptr;

  LUTPoint *newLUTPoint = new LUTPoint();
  newLUTPoint->density = density;
  newLUTPoint->sColor = *sColor;

  removeLUTPoint(density);

  if (!current) {
      this->headLUTPoint = newLUTPoint;
      newLUTPoint->next = nullptr;

      return VRL_OK;
  }

  while (current) {
      if (current->density > density) {
          if (previous) {
            previous->next = (void *)newLUTPoint;
            newLUTPoint->next = (void *)current;
          } else {
            newLUTPoint->next = (void *)(this->headLUTPoint);
            this->headLUTPoint = newLUTPoint;
          }

          return VRL_OK;
      }
      previous = current;
      current = (LUTPoint *)(current->next);
  }

  previous->next = (void *)newLUTPoint;
  newLUTPoint->next = nullptr;

  return VRL_OK;
}
  
int VRL::removeLUTPoint(float density) {
  LUTPoint *current = this->headLUTPoint;
  LUTPoint *previous = nullptr;

  while (current) {
      if (fabs(current->density - density) < EPSILON) {
          if (previous) {
              previous->next = current->next;
          } else {
              this->headLUTPoint = (LUTPoint *)(current->next);
          }
          delete current;
          return VRL_OK;
      }
      previous = current;
      current = (LUTPoint *)(current->next);
  }

  return VRL_NO_EFFECT;
}

int VRL::interpolateLUT() {
  deleteInterpolate();

  int32_t sizeInterp = this->maxDensity - this->minDensity;

  this->interpR         = new float[sizeInterp * 9];
  this->interpG         = &(this->interpR[sizeInterp * 1]);
  this->interpB         = &(this->interpR[sizeInterp * 2]);
  this->interpA         = &(this->interpR[sizeInterp * 3]);
  this->interpAmbient   = &(this->interpR[sizeInterp * 4]);
  this->interpDiffuse   = &(this->interpR[sizeInterp * 5]);
  this->interpSpecular  = &(this->interpR[sizeInterp * 6]);
  this->interpEmission  = &(this->interpR[sizeInterp * 7]);
  this->interpShininess = &(this->interpR[sizeInterp * 8]);

  LUTPoint *current = this->headLUTPoint;
  LUTPoint *previous = nullptr;
  uint32_t index = 0;
  while (current) {
      this->interpR[(int)(current->density - this->minDensity)] = current->sColor.r;
      this->interpG[(int)(current->density - this->minDensity)] = current->sColor.g;
      this->interpB[(int)(current->density - this->minDensity)] = current->sColor.b;
      this->interpA[(int)(current->density - this->minDensity)] = current->sColor.a;
      this->interpAmbient[(int)(current->density   - this->minDensity)] = current->sColor.ambient;
      this->interpDiffuse[(int)(current->density   - this->minDensity)] = current->sColor.diffuse;
      this->interpSpecular[(int)(current->density  - this->minDensity)] = current->sColor.specular;
      this->interpEmission[(int)(current->density  - this->minDensity)] = current->sColor.emission;
      this->interpShininess[(int)(current->density - this->minDensity)] = current->sColor.shininess;

      index++;
      previous = current;
      current = (LUTPoint *)current->next;
  }

  if (index < 2) {
    deleteInterpolate();
    return VRL_ERROR;
  }

  // left extrapolation
  for (int i = 0; i < (int)(this->headLUTPoint->density - this->minDensity); i++) {
    this->interpR[i] = this->headLUTPoint->sColor.r;
    this->interpG[i] = this->headLUTPoint->sColor.g;
    this->interpB[i] = this->headLUTPoint->sColor.b;
    this->interpA[i] = this->headLUTPoint->sColor.a;
    this->interpAmbient[i]   = this->headLUTPoint->sColor.ambient;
    this->interpDiffuse[i]   = this->headLUTPoint->sColor.diffuse;
    this->interpSpecular[i]  = this->headLUTPoint->sColor.specular;
    this->interpEmission[i]  = this->headLUTPoint->sColor.emission;
    this->interpShininess[i] = this->headLUTPoint->sColor.shininess;
  }

  // right extrapolation
  for (int i = (int)(previous->density - this->minDensity); i < sizeInterp; i++) {
    this->interpR[i] = previous->sColor.r;
    this->interpG[i] = previous->sColor.g;
    this->interpB[i] = previous->sColor.b;
    this->interpA[i] = previous->sColor.a;
    this->interpAmbient[i]   = previous->sColor.ambient;
    this->interpDiffuse[i]   = previous->sColor.diffuse;
    this->interpSpecular[i]  = previous->sColor.specular;
    this->interpEmission[i]  = previous->sColor.emission;
    this->interpShininess[i] = previous->sColor.shininess;
  }

  // linear interpolation
  current = this->headLUTPoint;
  uint32_t start;
  uint32_t stop;
  while (current != previous) {
    start = (uint32_t)( current->density - this->minDensity );
    stop  = (uint32_t)( ((LUTPoint *)current->next)->density - this->minDensity );
    for (uint32_t i = start; i < stop; i++) {
      this->interpR[i] = current->sColor.r + (i - start) * (previous->sColor.r - current->sColor.r) / (stop - start);
      this->interpG[i] = current->sColor.g + (i - start) * (previous->sColor.g - current->sColor.g) / (stop - start);
      this->interpB[i] = current->sColor.b + (i - start) * (previous->sColor.b - current->sColor.b) / (stop - start);
      this->interpA[i] = current->sColor.a + (i - start) * (previous->sColor.a - current->sColor.a) / (stop - start);
      this->interpAmbient[i]   = current->sColor.ambient   + (i - start) * (previous->sColor.ambient   - current->sColor.ambient)   / (stop - start);
      this->interpDiffuse[i]   = current->sColor.diffuse   + (i - start) * (previous->sColor.diffuse   - current->sColor.diffuse)   / (stop - start);
      this->interpSpecular[i]  = current->sColor.specular  + (i - start) * (previous->sColor.specular  - current->sColor.specular)  / (stop - start);
      this->interpEmission[i]  = current->sColor.emission  + (i - start) * (previous->sColor.emission  - current->sColor.emission)  / (stop - start);
      this->interpShininess[i] = current->sColor.shininess + (i - start) * (previous->sColor.shininess - current->sColor.shininess) / (stop - start);
    }
    current = (LUTPoint *)current->next;
  }

  return VRL_OK;
}

int VRL::setLight(float *direction, vrlColor *sColor) {
    this->lightDirection[0] = direction[0];
    this->lightDirection[1] = direction[1];
    this->lightDirection[2] = direction[2];
    this->lightDirection[3] = 0.0f;
    this->lightColor = *sColor;

    return VRL_OK;
}

int VRL::setImage(unsigned char *_image, uint32_t width, uint32_t height) {
  if (this->image && this->imageStatus == VRL_IMAGE_ALLOC) {
    delete[] (this->image);
  }

  this->image = _image;
  this->imageWidth = width;
  this->imageHeight = height;

  this->imageStatus = VRL_IMAGE_SET;
  return VRL_OK;
}

int VRL::allocImage(uint32_t width, uint32_t height) {
  if (this->image && this->imageStatus == VRL_IMAGE_ALLOC) {
    delete[] (this->image);
  }

  this->image = new unsigned char[4 * width * height];
  this->imageWidth = width;
  this->imageHeight = height;

  this->imageStatus = VRL_IMAGE_ALLOC;
  
  return VRL_OK;
}

unsigned char* VRL::getImage() {
  return (this->image);
}

void VRL::calcCameraMatrixs(uint32_t width, uint32_t height) {
  cameraTranslate(this->cameraPosition[0], this->cameraPosition[1], this->cameraPosition[2]);
  cameraRotate(this->cameraTarget, this->cameraUp);
  cameraPerspective(60.0f, width, height, zNear, 10.0f);
}

int VRL::setCameraConfigure(float *position, float *target, float *up, float zNear) {
  copyVec3f(this->cameraPosition, position);
  copyVec3f(this->cameraTarget, target);
  copyVec3f(this->cameraUp, up);

  this->zNear = zNear;

  return VRL_OK;
}

void VRL::cameraTranslate(float x, float y, float z) {
  setVec4f(&(this->cameraTranslateMatrix[0]),  1.0f, 0.0f, 0.0f, -x);
  setVec4f(&(this->cameraTranslateMatrix[4]),  0.0f, 1.0f, 0.0f, -y);
  setVec4f(&(this->cameraTranslateMatrix[8]),  0.0f, 0.0f, 1.0f, -z);
  setVec4f(&(this->cameraTranslateMatrix[12]), 0.0f, 0.0f, 0.0f, 1.0f);
}

void VRL::cameraRotate(float *_cameraTarget, float *_cameraUp) {
  float N[3] = { _cameraTarget[0], _cameraTarget[1], _cameraTarget[2] };
  normalizeVec3f(N);

  float U[3] = { _cameraUp[0], _cameraUp[1], _cameraUp[2] };
  normalizeVec3f(U);

  crossVec3f(U, _cameraTarget, U);
  
  float V[3];
  crossVec3f(N, U, V);

  setVec4f(&(this->cameraRotateMatrix[0]),  U[0], U[1], U[2], 0.0f);
  setVec4f(&(this->cameraRotateMatrix[4]),  V[0], V[1], V[2], 0.0f);
  setVec4f(&(this->cameraRotateMatrix[8]),  N[0], N[1], N[2], 0.0f);
  setVec4f(&(this->cameraRotateMatrix[12]), 0.0f, 0.0f, 0.0f, 1.0f);
}

void VRL::cameraPerspective(float angle, uint32_t width, uint32_t height, float _zNear, float zFar) {
  const double ar = (double)(width) / (double)(height);
  const float zRange = _zNear - zFar;
  const double tanHalfFOV = tan((double)ToRadian(angle / 2.0f));

  setVec4f(&(this->cameraPerspectiveMatrix[0]),  (float)(1.0L / (tanHalfFOV * ar)), 0.0f, 0.0f, 0.0f);
  setVec4f(&(this->cameraPerspectiveMatrix[4]),  0.0f, (float)(1.0L / tanHalfFOV), 0.0f, 0.0f);
  setVec4f(&(this->cameraPerspectiveMatrix[8]),  0.0f, 0.0f, (-_zNear - zFar) / zRange, 2.0f * zFar * _zNear / zRange);
  setVec4f(&(this->cameraPerspectiveMatrix[12]), 0.0f, 0.0f, 1.0f, 0.0f);
}

int VRL::setRotateVolume(float angle, float *axis) {
  double rotationAngle = (double)ToRadian(angle);
  double qx = (double)(axis[0]) * sin(rotationAngle / 2.0L);
  double qy = (double)(axis[1]) * sin(rotationAngle / 2.0L);
  double qz = (double)(axis[2]) * sin(rotationAngle / 2.0L);
  double qw = cos(rotationAngle / 2.0L);

  this->volumeRotateMatrix[0]  = (float)(1.0L - 2.0L * qy*qy - 2.0L * qz * qz);
  this->volumeRotateMatrix[1]  = (float)(2.0L * qx*qy - 2.0L * qz*qw);
  this->volumeRotateMatrix[2]  = (float)(2.0L * qx*qz + 2.0L * qy*qw);
  this->volumeRotateMatrix[3]  = 0.0f;

  this->volumeRotateMatrix[4]  = (float)(2.0L * qx*qy + 2.0L * qz*qw);
  this->volumeRotateMatrix[5]  = (float)(1.0L - 2.0L * qx*qx - 2.0L * qz*qz);
  this->volumeRotateMatrix[6]  = (float)(2.0L * qy*qz - 2.0L * qx*qw);
  this->volumeRotateMatrix[7]  = 0.0f;

  this->volumeRotateMatrix[8]  = (float)(2.0L * qx*qz - 2.0L * qy*qw);
  this->volumeRotateMatrix[9]  = (float)(2.0L * qy*qz + 2.0L * qx*qw);
  this->volumeRotateMatrix[10] = (float)(1.0L - 2.0L * qx*qx - 2.0L * qy*qy);
  this->volumeRotateMatrix[11] = 0.0f;

  this->volumeRotateMatrix[12] = 0.0f;
  this->volumeRotateMatrix[13] = 0.0f;
  this->volumeRotateMatrix[14] = 0.0f;
  this->volumeRotateMatrix[15] = 1.0f;

  return VRL_OK;
}

int VRL::addRotateVolume(float angle, float *axis) {
    float oldRotateVolumeMatrix[16];

    for (int32_t i = 0; i < 16; i++) {
      oldRotateVolumeMatrix[i] = this->volumeRotateMatrix[i];
    }

    setRotateVolume(angle, axis);

    multMatrix4fx4f(this->volumeRotateMatrix, oldRotateVolumeMatrix, this->volumeRotateMatrix);

    return VRL_OK;
}

int VRL::resetRotationVolume() {
  setVec4f(&(this->volumeRotateMatrix[0]),  1.0f, 0.0f, 0.0f, 0.0f);
  setVec4f(&(this->volumeRotateMatrix[4]),  0.0f, 1.0f, 0.0f, 0.0f);
  setVec4f(&(this->volumeRotateMatrix[8]),  0.0f, 0.0f, 1.0f, 0.0f);
  setVec4f(&(this->volumeRotateMatrix[12]), 0.0f, 0.0f, 0.0f, 1.0f);

  return VRL_OK;
}

int VRL::setTranslateVolume(float x, float y, float z) {
  setVec4f(&(this->volumeTranslateMatrix[0]),  1.0f, 0.0f, 0.0f, x);
  setVec4f(&(this->volumeTranslateMatrix[4]),  0.0f, 1.0f, 0.0f, y);
  setVec4f(&(this->volumeTranslateMatrix[8]),  0.0f, 0.0f, 1.0f, z);
  setVec4f(&(this->volumeTranslateMatrix[12]), 0.0f, 0.0f, 0.0f, 1.0f);

  return VRL_OK;
}

int VRL::resetTranslateVolume() {
  int ret = setTranslateVolume(0.0f, 0.0f, 0.0f);
  return ret;
}

int VRL::setImageQuality(VRL_IMAGE_QUALITY imgQua) {
  switch (imgQua) {
  case VRL_IMG_QUA_LOW:
    this->renderStepSize = 4.0f;
    break;

  case VRL_IMG_QUA_MEDIUM:
    this->renderStepSize = 2.0f;
    break;
  
  case VRL_IMG_QUA_BEST:
    this->renderStepSize = 1.0f;
    break;

  case VRL_IMG_QUA_BEST_EXPERIMENTAL:
    this->renderStepSize = 0.5f;
    break;

  case VRL_IMG_QUA_REALISTIC:
    this->renderStepSize = 0.01f;
  break;

  default:
    this->renderStepSize = 1.0f;
    break;
  }

  return VRL_OK;
}

int VRL::enableAntialiasing(ANTIALIASING_VALUE value) {
  switch (value) {
  case VRL_ANTIALIASING_X1:
    this->antialiasingValue = VRL_ANTIALIASING_X1;
    this->imageScale = 1;
    break;

  case VRL_ANTIALIASING_X4:
    this->antialiasingValue = VRL_ANTIALIASING_X4;
    this->imageScale = 2;
    break;

  case VRL_ANTIALIASING_X16:
    this->antialiasingValue = VRL_ANTIALIASING_X16;
    this->imageScale = 4;
    break;

  default:
    this->antialiasingValue = VRL_ANTIALIASING_X1;
    this->imageScale = 1;
    break;
  }

  return VRL_OK;
}

/*******************
 * RENDER funtions *
 *******************/

float* VRL::renderGetMVPMatrix() {
  calcCameraMatrixs(this->imageWidth * this->imageScale, this->imageHeight * this->imageScale);

  // M: model matrix
  float modelMatrix[16];
  multMatrix4fx4f(this->volumeTranslateMatrix, this->volumeRotateMatrix, modelMatrix);

  // V: view matrix 
  float viewMatrix[16];
  multMatrix4fx4f(this->cameraRotateMatrix, this->cameraTranslateMatrix, viewMatrix);

  // P: projection matrix
  // this->cameraPerspectiveMatrix;

  // MVP matrix
  float *mvpMatrix = new float[16];
  multMatrix4fx4f(this->cameraPerspectiveMatrix, viewMatrix, mvpMatrix);
  multMatrix4fx4f(mvpMatrix, modelMatrix, mvpMatrix);

  return mvpMatrix;
}

float* VRL::renderCreateBoxForVolume(){
  float *pointsTrianglesForBox = new float[108]; // (3 points) * (3 value for point(xyz)) * (2 triangles for one face) * (6 faces)

  // normalize the length of edges
  float maxLen = (float)((this->xlen > this->ylen) ? (this->xlen) : (this->ylen));
  maxLen = (float)((this->zlen > maxLen) ? (this->zlen) : (maxLen));
  float boxX = (float)(this->xlen) / maxLen;
  float boxY = (float)(this->ylen) / maxLen;
  float boxZ = (float)(this->zlen) / maxLen;
  
  float boxXdiv2 = boxX / 2.0f;
  float boxYdiv2 = boxY / 2.0f;
  float boxZdiv2 = boxZ / 2.0f;
  
  // copy for use in fragment shader
  this->volumeBoxMaxLen = maxLen;

  // bottom triangles
  // first
  pointsTrianglesForBox[0]  = -boxXdiv2; pointsTrianglesForBox[1]  = -boxYdiv2; pointsTrianglesForBox[2]  = -boxZdiv2; 
  pointsTrianglesForBox[3]  =  boxXdiv2; pointsTrianglesForBox[4]  = -boxYdiv2; pointsTrianglesForBox[5]  = -boxZdiv2; 
  pointsTrianglesForBox[6]  =  boxXdiv2; pointsTrianglesForBox[7]  = -boxYdiv2; pointsTrianglesForBox[8]  =  boxZdiv2;
  // second
  pointsTrianglesForBox[9]  = -boxXdiv2; pointsTrianglesForBox[10] = -boxYdiv2; pointsTrianglesForBox[11] = -boxZdiv2;
  pointsTrianglesForBox[12] = -boxXdiv2; pointsTrianglesForBox[13] = -boxYdiv2; pointsTrianglesForBox[14] =  boxZdiv2; 
  pointsTrianglesForBox[15] =  boxXdiv2; pointsTrianglesForBox[16] = -boxYdiv2; pointsTrianglesForBox[17] =  boxZdiv2;

  // top triangles
  // first
  pointsTrianglesForBox[18] = -boxXdiv2; pointsTrianglesForBox[19] =  boxYdiv2; pointsTrianglesForBox[20] = -boxZdiv2; 
  pointsTrianglesForBox[21] =  boxXdiv2; pointsTrianglesForBox[22] =  boxYdiv2; pointsTrianglesForBox[23] = -boxZdiv2; 
  pointsTrianglesForBox[24] =  boxXdiv2; pointsTrianglesForBox[25] =  boxYdiv2; pointsTrianglesForBox[26] =  boxZdiv2; 
  // second
  pointsTrianglesForBox[27] = -boxXdiv2; pointsTrianglesForBox[28] =  boxYdiv2; pointsTrianglesForBox[29] = -boxZdiv2;
  pointsTrianglesForBox[30] = -boxXdiv2; pointsTrianglesForBox[31] =  boxYdiv2; pointsTrianglesForBox[32] =  boxZdiv2; 
  pointsTrianglesForBox[33] =  boxXdiv2; pointsTrianglesForBox[34] =  boxYdiv2; pointsTrianglesForBox[35] =  boxZdiv2;

  // front triangles
  // first
  pointsTrianglesForBox[36] = -boxXdiv2; pointsTrianglesForBox[37] = -boxYdiv2; pointsTrianglesForBox[38] = -boxZdiv2;
  pointsTrianglesForBox[39] = -boxXdiv2; pointsTrianglesForBox[40] =  boxYdiv2; pointsTrianglesForBox[41] = -boxZdiv2;
  pointsTrianglesForBox[42] =  boxXdiv2; pointsTrianglesForBox[43] =  boxYdiv2; pointsTrianglesForBox[44] = -boxZdiv2;
  // second
  pointsTrianglesForBox[45] = -boxXdiv2; pointsTrianglesForBox[46] = -boxYdiv2; pointsTrianglesForBox[47] = -boxZdiv2;
  pointsTrianglesForBox[48] =  boxXdiv2; pointsTrianglesForBox[49] = -boxYdiv2; pointsTrianglesForBox[50] = -boxZdiv2;
  pointsTrianglesForBox[51] =  boxXdiv2; pointsTrianglesForBox[52] =  boxYdiv2; pointsTrianglesForBox[53] = -boxZdiv2;

  // back triangles
  // first
  pointsTrianglesForBox[54] = -boxXdiv2; pointsTrianglesForBox[55] = -boxYdiv2; pointsTrianglesForBox[56] =  boxZdiv2;
  pointsTrianglesForBox[57] = -boxXdiv2; pointsTrianglesForBox[58] =  boxYdiv2; pointsTrianglesForBox[59] =  boxZdiv2;
  pointsTrianglesForBox[60] =  boxXdiv2; pointsTrianglesForBox[61] =  boxYdiv2; pointsTrianglesForBox[62] =  boxZdiv2;
  // second
  pointsTrianglesForBox[63] = -boxXdiv2; pointsTrianglesForBox[64] = -boxYdiv2; pointsTrianglesForBox[65] =  boxZdiv2;
  pointsTrianglesForBox[66] =  boxXdiv2; pointsTrianglesForBox[67] = -boxYdiv2; pointsTrianglesForBox[68] =  boxZdiv2;
  pointsTrianglesForBox[69] =  boxXdiv2; pointsTrianglesForBox[70] =  boxYdiv2; pointsTrianglesForBox[71] =  boxZdiv2;

  // left triangles
  // first
  pointsTrianglesForBox[72] = -boxXdiv2; pointsTrianglesForBox[73] = -boxYdiv2; pointsTrianglesForBox[74] = -boxZdiv2;
  pointsTrianglesForBox[75] = -boxXdiv2; pointsTrianglesForBox[76] =  boxYdiv2; pointsTrianglesForBox[77] = -boxZdiv2;
  pointsTrianglesForBox[78] = -boxXdiv2; pointsTrianglesForBox[79] =  boxYdiv2; pointsTrianglesForBox[80] =  boxZdiv2;
  // second
  pointsTrianglesForBox[81] = -boxXdiv2; pointsTrianglesForBox[82] = -boxYdiv2; pointsTrianglesForBox[83] = -boxZdiv2;
  pointsTrianglesForBox[84] = -boxXdiv2; pointsTrianglesForBox[85] = -boxYdiv2; pointsTrianglesForBox[86] =  boxZdiv2;
  pointsTrianglesForBox[87] = -boxXdiv2; pointsTrianglesForBox[88] =  boxYdiv2; pointsTrianglesForBox[89] =  boxZdiv2;

  // right triangles
  // first
  pointsTrianglesForBox[90] =  boxXdiv2; pointsTrianglesForBox[91] = -boxYdiv2; pointsTrianglesForBox[92] = -boxZdiv2;
  pointsTrianglesForBox[93] =  boxXdiv2; pointsTrianglesForBox[94] =  boxYdiv2; pointsTrianglesForBox[95] = -boxZdiv2;
  pointsTrianglesForBox[96] =  boxXdiv2; pointsTrianglesForBox[97] =  boxYdiv2; pointsTrianglesForBox[98] =  boxZdiv2;
  // second
  pointsTrianglesForBox[99] =  boxXdiv2; pointsTrianglesForBox[100]= -boxYdiv2; pointsTrianglesForBox[101]= -boxZdiv2;
  pointsTrianglesForBox[102]=  boxXdiv2; pointsTrianglesForBox[103]= -boxYdiv2; pointsTrianglesForBox[104]=  boxZdiv2;
  pointsTrianglesForBox[105]=  boxXdiv2; pointsTrianglesForBox[106]=  boxYdiv2; pointsTrianglesForBox[107]=  boxZdiv2;

  return pointsTrianglesForBox;
}

float sign(float *p1, float *p2, float *p3) {
    return (p1[0] - p3[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p3[1]);
}

int32_t pointInTriangle(float *pt, float *v1, float *v2, float *v3) {
    int32_t b1, b2, b3;

    b1 = sign(pt, v1, v2) < 0.0f;
    b2 = sign(pt, v2, v3) < 0.0f;
    b3 = sign(pt, v3, v1) < 0.0f;

    return ((b1 == b2) && (b2 == b3));
}

void getBarycentricCoordinate(float *barycentricCoordinate, float *pt, float *v1, float *v2, float *v3) {
  barycentricCoordinate[0] = ((v2[1] - v3[1]) * (pt[0] - v3[0]) + (v3[0] - v2[0]) * (pt[1] - v3[1])) / ((v2[1] - v3[1]) * (v1[0] - v3[0]) + (v3[0] - v2[0]) * (v1[1] - v3[1]));
  barycentricCoordinate[1] = ((v3[1] - v1[1]) * (pt[0] - v3[0]) + (v1[0] - v3[0]) * (pt[1] - v3[1])) / ((v2[1] - v3[1]) * (v1[0] - v3[0]) + (v3[0] - v2[0]) * (v1[1] - v3[1]));
  barycentricCoordinate[2] = 1.0f - barycentricCoordinate[0] - barycentricCoordinate[1];

  barycentricCoordinate[3] = v1[3];
  barycentricCoordinate[4] = v2[3];
  barycentricCoordinate[5] = v3[3];
}

float barycentricInterp(float *barycentricCoordinate, float value1, float value2, float value3) {
  float w1, w2, w3;
  w1 = barycentricCoordinate[3];
  w2 = barycentricCoordinate[4];
  w3 = barycentricCoordinate[5];
  return ((value1 * barycentricCoordinate[0]) / w1  + (value2 * barycentricCoordinate[1]) / w2 + (value3 * barycentricCoordinate[2]) / w3) /
         (barycentricCoordinate[0] / w1 + barycentricCoordinate[1] / w2 + barycentricCoordinate[2] / w3);
}

float cubicInterpolate(float *array4, float x) {
  float x2 = x * x;
  float x3 = x2 * x;

  return array4[1] + (-0.5f * array4[0] + 0.5f * array4[2]) * x
    + (array4[0] - 2.5f * array4[1] + 2.0f * array4[2] - 0.5f * array4[3]) * x2
    + (-0.5f * array4[0] + 1.5f * array4[1] - 1.5f * array4[2] + 0.5f * array4[3]) * x3;
}

float threeCubicInterpolation(float *array64, float *pointPosition) {
  float array4[4];
  float array4Final[4];

  array4[0] = cubicInterpolate(&(array64[0]), pointPosition[0]);
  array4[1] = cubicInterpolate(&(array64[4]), pointPosition[0]);
  array4[2] = cubicInterpolate(&(array64[8]), pointPosition[0]);
  array4[3] = cubicInterpolate(&(array64[12]), pointPosition[0]);

  array4Final[0] = cubicInterpolate(array4, pointPosition[1]);

  array4[0] = cubicInterpolate(&(array64[16]), pointPosition[0]);
  array4[1] = cubicInterpolate(&(array64[20]), pointPosition[0]);
  array4[2] = cubicInterpolate(&(array64[24]), pointPosition[0]);
  array4[3] = cubicInterpolate(&(array64[28]), pointPosition[0]);

  array4Final[1] = cubicInterpolate(array4, pointPosition[1]);

  array4[0] = cubicInterpolate(&(array64[32]), pointPosition[0]);
  array4[1] = cubicInterpolate(&(array64[36]), pointPosition[0]);
  array4[2] = cubicInterpolate(&(array64[40]), pointPosition[0]);
  array4[3] = cubicInterpolate(&(array64[44]), pointPosition[0]);

  array4Final[2] = cubicInterpolate(array4, pointPosition[1]);

  array4[0] = cubicInterpolate(&(array64[48]), pointPosition[0]);
  array4[1] = cubicInterpolate(&(array64[52]), pointPosition[0]);
  array4[2] = cubicInterpolate(&(array64[56]), pointPosition[0]);
  array4[3] = cubicInterpolate(&(array64[60]), pointPosition[0]);

  array4Final[3] = cubicInterpolate(array4, pointPosition[1]);

  return cubicInterpolate(array4Final, pointPosition[2]);
}

float VRL::getDensityFromVolume(float x, float y, float z) {
  int32_t xInt = (int32_t)(x);
  int32_t yInt = (int32_t)(y);
  int32_t zInt = (int32_t)(z);

  float array64[64];
  float localPointPosition[3];
  float density;
  
  if ((xInt - 1) >= 0 && ((uint32_t)xInt - 1) + 3 < this->xlen &&
      (yInt - 1) >= 0 && ((uint32_t)yInt - 1) + 3 < this->ylen &&
      (zInt - 1) >= 0 && ((uint32_t)zInt - 1) + 3 < this->zlen) {
    int32_t index = 0;
    for (int32_t startZposition = zInt - 1; startZposition < (zInt - 1) + 4; startZposition++) {
        for (int32_t startYposition = yInt - 1; startYposition < (yInt - 1) + 4; startYposition++) {
            for (int32_t startXposition = xInt - 1; startXposition < (xInt - 1) + 4; startXposition++) {
                array64[index] = this->volume[this->xOffset[startXposition] + this->yOffset[startYposition] + this->zOffset[startZposition]];
                index++;
            }
        }
    }

    localPointPosition[0] = modff(x, &density);
    localPointPosition[1] = modff(y, &density);
    localPointPosition[2] = modff(z, &density);

    density = threeCubicInterpolation(array64, localPointPosition);
    
    density = (density < this->minDensity) ? (this->minDensity) : (density);
    density = (density > this->maxDensity - 1.0f) ? (this->maxDensity - 1.0f) : (density);

    return density;
  }

  return (float)this->minDensity;
}

int32_t VRL::renderFragmentShader(float *rayPosition, unsigned char *pixelColor) {
  float rayStepSize = this->renderStepSize / this->volumeBoxMaxLen;
  float rayDirection[3];
  float newRayPosition[3];
  int32_t inVolume;
  int32_t needDrawPixel;
  float density;

  // M: model matrix
  float modelMatrix[16];
  multMatrix4fx4f(this->volumeTranslateMatrix, this->volumeRotateMatrix, modelMatrix);
  // V: view matrix 
  float viewMatrix[16];
  multMatrix4fx4f(this->cameraRotateMatrix, this->cameraTranslateMatrix, viewMatrix);

  float MV[16];
  multMatrix4fx4f(viewMatrix, modelMatrix, MV);

  float vecPoint[4];
  float newCamera[3];
  vecPoint[0] = 0.0f - MV[3];
  vecPoint[1] = 0.0f - MV[7];
  vecPoint[2] = 0.0f - MV[11];
  vecPoint[3] = 1.0f - MV[15];

  newCamera[0] = vecPoint[0] * MV[0] + vecPoint[1] * MV[4] + vecPoint[2] * MV[8] + vecPoint[3] * MV[12];
  newCamera[1] = vecPoint[0] * MV[1] + vecPoint[1] * MV[5] + vecPoint[2] * MV[9] + vecPoint[3] * MV[13];
  newCamera[2] = vecPoint[0] * MV[2] + vecPoint[1] * MV[6] + vecPoint[2] * MV[10] + vecPoint[3] * MV[14];

  //c4 point_transformed_with_inverse = vec4(vec3((point - matrix[3]) * matrix), 1.0);

  // calculate ray direction for current point
  rayDirection[0] = rayPosition[0] - newCamera[0];
  rayDirection[1] = rayPosition[1] - newCamera[1];
  rayDirection[2] = rayPosition[2] - newCamera[2];

  normalizeVec3f(rayDirection);

  rayDirection[0] = rayStepSize * rayDirection[0];
  rayDirection[1] = rayStepSize * rayDirection[1];
  rayDirection[2] = rayStepSize * rayDirection[2];

  float boxX = (float)(this->xlen) /  this->volumeBoxMaxLen;
  float boxY = (float)(this->ylen) /  this->volumeBoxMaxLen;
  float boxZ = (float)(this->zlen) /  this->volumeBoxMaxLen;
  float boxXdiv2 = boxX / 2.0f;
  float boxYdiv2 = boxY / 2.0f;
  float boxZdiv2 = boxZ / 2.0f;

  newRayPosition[0] = rayPosition[0] + rayDirection[0];
  newRayPosition[1] = rayPosition[1] + rayDirection[1];
  newRayPosition[2] = rayPosition[2] + rayDirection[2];

  inVolume = 1;
  needDrawPixel = 0;
  int32_t index;
  float rgb[3];
  float diffuseFactor;
  float lightDirectionMinus[3];
  lightDirectionMinus[0] = -this->lightDirection[0];
  lightDirectionMinus[1] = -this->lightDirection[1];
  lightDirectionMinus[2] = -this->lightDirection[2];
  float currentNormal[3];
  float densityNeighbors[6]; // +-x; +-y; +-z
  float specularFactor;

  while (inVolume) {
    // inside step
    newRayPosition[0] = newRayPosition[0] + rayDirection[0];
    newRayPosition[1] = newRayPosition[1] + rayDirection[1];
    newRayPosition[2] = newRayPosition[2] + rayDirection[2];

    // we are still inside?
    if (newRayPosition[0] > -boxXdiv2 && newRayPosition[0] < boxXdiv2 &&
        newRayPosition[1] > -boxYdiv2 && newRayPosition[1] < boxYdiv2 &&
        newRayPosition[2] > -boxZdiv2 && newRayPosition[2] < boxZdiv2) {
      // get density
      density = getDensityFromVolume((newRayPosition[0] + boxXdiv2) * this->volumeBoxMaxLen,  // x
                                     (newRayPosition[1] + boxYdiv2) * this->volumeBoxMaxLen,  // y
                                     (newRayPosition[2] + boxZdiv2) * this->volumeBoxMaxLen); // z


// TEST START
      if (density > 200.0f) {
        // Phong 
        index = (int32_t)(density - this->minDensity);
        // Ambient
        rgb[0] = this->interpR[index] * this->interpAmbient[index];
        rgb[1] = this->interpG[index] * this->interpAmbient[index];
        rgb[2] = this->interpB[index] * this->interpAmbient[index];
        // Diffuse
        // +x -x +y -y +z -z
        densityNeighbors[0] = getDensityFromVolume((newRayPosition[0] + boxXdiv2 + rayStepSize) * this->volumeBoxMaxLen, (newRayPosition[1] + boxYdiv2) * this->volumeBoxMaxLen, (newRayPosition[2] + boxZdiv2) * this->volumeBoxMaxLen);
        densityNeighbors[1] = getDensityFromVolume((newRayPosition[0] + boxXdiv2 - rayStepSize) * this->volumeBoxMaxLen, (newRayPosition[1] + boxYdiv2) * this->volumeBoxMaxLen, (newRayPosition[2] + boxZdiv2) * this->volumeBoxMaxLen);
        densityNeighbors[2] = getDensityFromVolume((newRayPosition[0] + boxXdiv2) * this->volumeBoxMaxLen, (newRayPosition[1] + boxYdiv2 + rayStepSize) * this->volumeBoxMaxLen, (newRayPosition[2] + boxZdiv2) * this->volumeBoxMaxLen);
        densityNeighbors[3] = getDensityFromVolume((newRayPosition[0] + boxXdiv2) * this->volumeBoxMaxLen, (newRayPosition[1] + boxYdiv2 - rayStepSize) * this->volumeBoxMaxLen, (newRayPosition[2] + boxZdiv2) * this->volumeBoxMaxLen);
        densityNeighbors[4] = getDensityFromVolume((newRayPosition[0] + boxXdiv2) * this->volumeBoxMaxLen, (newRayPosition[1] + boxYdiv2) * this->volumeBoxMaxLen, (newRayPosition[2] + boxZdiv2 + rayStepSize) * this->volumeBoxMaxLen);
        densityNeighbors[5] = getDensityFromVolume((newRayPosition[0] + boxXdiv2) * this->volumeBoxMaxLen, (newRayPosition[1] + boxYdiv2) * this->volumeBoxMaxLen, (newRayPosition[2] + boxZdiv2 - rayStepSize) * this->volumeBoxMaxLen);
        currentNormal[0] = -(densityNeighbors[0] - densityNeighbors[1]);
        currentNormal[1] = -(densityNeighbors[2] - densityNeighbors[3]);
        currentNormal[2] = -(densityNeighbors[4] - densityNeighbors[5]);
        normalizeVec3f(currentNormal);
        matrix4fMultVec4f(modelMatrix, currentNormal, currentNormal);
        //if(dot(normalVec,ray)>0) normalVec*=-1;
        dotVec3f(currentNormal, lightDirectionMinus, &diffuseFactor);
        diffuseFactor = (diffuseFactor < 0.0f)?(0.0f):(diffuseFactor);
        rgb[0] += this->interpR[index] * this->interpDiffuse[index] * diffuseFactor;
        rgb[1] += this->interpG[index] * this->interpDiffuse[index] * diffuseFactor;
        rgb[2] += this->interpB[index] * this->interpDiffuse[index] * diffuseFactor;
        // Specular
        //specularFactor = ;
        //rgb[0] += this->interpR[index] * this->interpSpecular[index] * specularFactor;
        //rgb[1] += this->interpG[index] * this->interpSpecular[index] * specularFactor;
        //rgb[2] += this->interpB[index] * this->interpSpecular[index] * specularFactor;

        rgb[0] = (rgb[0] < 1.0f)?(rgb[0]):(1.0f);
        rgb[1] = (rgb[1] < 1.0f)?(rgb[1]):(1.0f);
        rgb[2] = (rgb[2] < 1.0f)?(rgb[2]):(1.0f);
        pixelColor[0] = (unsigned char)(rgb[0] * 255.0f);
        pixelColor[1] = (unsigned char)(rgb[1] * 255.0f);
        pixelColor[2] = (unsigned char)(rgb[2] * 255.0f);
        pixelColor[3] = (unsigned char)(this->interpA[(int)(density - this->minDensity)] * 255.0f);
        //this->interpAmbient[(int)(density - this->minDensity)];
        //this->interpDiffuse[(int)(density - this->minDensity)];
        //this->interpSpecular[(int)(density - this->minDensity)];
        //this->interpEmission[(int)(density - this->minDensity)];
        //this->interpShininess[(int)(density - this->minDensity)];
        needDrawPixel = 1;
        break;
      }
// TEST END

    // TODO:

    } else {
      break;
    }
  }

  if (needDrawPixel) {
    //printf("%1.6f %1.6f %1.6f\n", rayPosition[0], rayPosition[1], rayPosition[2]);
    //pixelColor[0] = 255;
    //pixelColor[1] = 255;
    //pixelColor[2] = 255;
    //pixelColor[3] = 255;

    return 1;
  }

  return 0;
}

void VRL::renderPipeLine(float *mvpMatrix, float *pointsTriangle, float *depthBuffer) {
  // Screen coordinates
  float xImg;
  float yImg;
  // Vertices of the triangle
  float inPoint[4];
  float outPoints[12]; // [0]:x, [1]:y, [3]:depth

  // Calculate the position on the screen and the depth of each point of the triangle
  for (uint32_t i = 0; i < 3; i++) {
    inPoint[0] = pointsTriangle[i * 3 + 0];
    inPoint[1] = pointsTriangle[i * 3 + 1];
    inPoint[2] = pointsTriangle[i * 3 + 2];
    inPoint[3] = 1.0f;

    matrix4fMultVec4f(mvpMatrix, inPoint, &(outPoints[i * 4]));

    xImg = outPoints[i * 4 + 0] / outPoints[i * 4 + 3];
    yImg = outPoints[i * 4 + 1] / outPoints[i * 4 + 3];

    xImg = (xImg + 1.0f) * (float)(this->imageWidth * this->imageScale) / 2.0f;
    yImg = (float)(this->imageHeight * this->imageScale) - (yImg + 1.0f) * (float)(this->imageHeight * this->imageScale) / 2.0f;

    xImg = roundf(xImg);
    yImg = roundf(yImg);

    outPoints[i * 4 + 0] = xImg;
    outPoints[i * 4 + 1] = yImg;
  }

  // Rasterizer
  // find box for triangle
  int32_t xMin = (int32_t)outPoints[0];
  int32_t xMax = (int32_t)outPoints[0];
  int32_t yMin = (int32_t)outPoints[1];
  int32_t yMax = (int32_t)outPoints[1];
  
  for (uint32_t i = 1; i < 3; i++) {
    xMin = (xMin > outPoints[i * 4]) ? ((int32_t)outPoints[i * 4]):(xMin);
    xMax = (xMax < outPoints[i * 4]) ? ((int32_t)outPoints[i * 4]):(xMax);
    
    yMin = (yMin > outPoints[i * 4 + 1]) ? ((int32_t)outPoints[i * 4 + 1]):(yMin);
    yMax = (yMax < outPoints[i * 4 + 1]) ? ((int32_t)outPoints[i * 4 + 1]):(yMax);
  }

  unsigned char pixelColor[4];
  int32_t pointInTriangleBool;
  int32_t fragmentShaderReturnBool;
  float interpDepth;
  float testPoint[2];

  // M: model matrix
  //float modelMatrix[16];
  //multMatrix4x4(this->volumeTranslateMatrix, this->volumeRotateMatrix, modelMatrix);

  // Barycentric coordinate
  float barycentricCoordinate[6]; // l1, l2, l3, w1, w2, w3
  float worldPosition[3];

  for (int32_t j = yMin; j < yMax + 1; j++) {
    for (int32_t i = xMin; i < xMax + 1; i++) {
      // if a point on the screen
      if (i >= 0 && (uint32_t)i < (this->imageWidth * this->imageScale) && j >= 0 && (uint32_t)j < (this->imageHeight * this->imageScale)) {
        testPoint[0] = (float)i;
        testPoint[1] = (float)j;
        pointInTriangleBool = pointInTriangle(testPoint, &(outPoints[0]), &(outPoints[4]), &(outPoints[8]));
        // if the point belongs to the triangle
        if (pointInTriangleBool) {
          getBarycentricCoordinate(barycentricCoordinate, testPoint, &(outPoints[0]), &(outPoints[4]), &(outPoints[8]));
          // interpDepth
          interpDepth = barycentricInterp(barycentricCoordinate, outPoints[3], outPoints[7], outPoints[11]);
          // depth buffer test
          if (interpDepth < depthBuffer[((uint32_t)j) * (this->imageWidth * this->imageScale) + (uint32_t)i]) {
            depthBuffer[((uint32_t)j) * (this->imageWidth * this->imageScale) + (uint32_t)i] = interpDepth;
            // calculate world position for current point
            worldPosition[0] = barycentricInterp(barycentricCoordinate, pointsTriangle[0 * 3 + 0], pointsTriangle[1 * 3 + 0], pointsTriangle[2 * 3 + 0]); // x
            worldPosition[1] = barycentricInterp(barycentricCoordinate, pointsTriangle[0 * 3 + 1], pointsTriangle[1 * 3 + 1], pointsTriangle[2 * 3 + 1]); // y
            worldPosition[2] = barycentricInterp(barycentricCoordinate, pointsTriangle[0 * 3 + 2], pointsTriangle[1 * 3 + 2], pointsTriangle[2 * 3 + 2]); // z
            // Fragment Shader
            fragmentShaderReturnBool = renderFragmentShader(worldPosition, pixelColor);
            // Set color to image
            if (fragmentShaderReturnBool) {
              this->privateRenderimage[((uint32_t)j * (this->imageWidth * this->imageScale) + (uint32_t)i) * 4 + 0] = pixelColor[0];
              this->privateRenderimage[((uint32_t)j * (this->imageWidth * this->imageScale) + (uint32_t)i) * 4 + 1] = pixelColor[1];
              this->privateRenderimage[((uint32_t)j * (this->imageWidth * this->imageScale) + (uint32_t)i) * 4 + 2] = pixelColor[2];
              this->privateRenderimage[((uint32_t)j * (this->imageWidth * this->imageScale) + (uint32_t)i) * 4 + 3] = pixelColor[3];
            }
          }
        }
      }
    }
  }

  return;
}

int VRL::render() {
  uint32_t imageSize = 4 * (this->imageWidth) * (this->imageHeight);
  uint32_t privateImageSize = imageSize * (this->imageScale) * (this->imageScale);
  
  // Init inside buffer for render
  this->privateRenderimage = new unsigned char[privateImageSize];

  for (uint32_t i = 0; i < privateImageSize; i+=4) {
    privateRenderimage[i + 0] = 0; // r
    privateRenderimage[i + 1] = 0; // g
    privateRenderimage[i + 2] = 0; // b
    privateRenderimage[i + 3] = 0; // a
  }

  // MVP matrix
  float *mvpMatrix = renderGetMVPMatrix();

  // Create box for volume
  float *pointsTrianglesForBox = renderCreateBoxForVolume();

  // Create Z-buffer
  float *depthBuffer = new float[privateImageSize / 4];

  for (uint32_t i = 0; i < privateImageSize / 4; i++) {
    depthBuffer[i] = 1000.0f;
  }

  // for all 12 triangles
  for (uint32_t i = 0; i < 108; i+=9) {
    //printf("%u Triangle\n", i / 9);
    renderPipeLine(mvpMatrix, &(pointsTrianglesForBox[i]), depthBuffer);
  }

  delete[]mvpMatrix;
  delete[]pointsTrianglesForBox;
  delete[]depthBuffer;

  if (this->enableDrawBox) {
    renderDrawBox();
  }

#if 0
//TEST START
  uint32_t cubePointNumber[8] = {18, 21, 3, 0, 30, 33, 15, 12};
  float vector[4];
  float outPoints[4 * 8];
  float x;
  float y;
  for (uint32_t i = 0; i < 8; i++) {
    vector[0] = pointsTrianglesForBox[cubePointNumber[i] + 0];
    vector[1] = pointsTrianglesForBox[cubePointNumber[i] + 1];
    vector[2] = pointsTrianglesForBox[cubePointNumber[i] + 2];
    vector[3] = 1.0f;

    matrix4MultVector4(mvpMatrix, vector, &(outPoints[i * 4]));

    // set point
    x = outPoints[i * 4 + 0] / outPoints[i * 4 + 3];
    y = outPoints[i * 4 + 1] / outPoints[i * 4 + 3];
    
    x = (x + 1.0f) * (float)(this->imageWidth * this->imageScale) / 2.0f;
    y = (float)(this->imageHeight * this->imageScale) - (y + 1.0f) * (float)(this->imageHeight * this->imageScale) / 2.0f;
    
    x = roundf(x);
    y = roundf(y);
    
    if ((uint32_t)x < (this->imageWidth * this->imageScale) && (uint32_t)y < (this->imageHeight * this->imageScale)) {
      this->privateRenderimage[((uint32_t)y * (this->imageWidth * this->imageScale) + (uint32_t)x) * 4 + 0] = 255;
      this->privateRenderimage[((uint32_t)y * (this->imageWidth * this->imageScale) + (uint32_t)x) * 4 + 1] = 255;
      this->privateRenderimage[((uint32_t)y * (this->imageWidth * this->imageScale) + (uint32_t)x) * 4 + 2] = 255;
      this->privateRenderimage[((uint32_t)y * (this->imageWidth * this->imageScale) + (uint32_t)x) * 4 + 3] = 255;
    }
  }
//TEST STOP 
#endif

  // output image (RGBA8888 4 byte per pixel)
  float tmpColor;
  for (uint32_t y = 0; y < this->imageHeight; y++) {
    for (uint32_t x = 0; x < this->imageWidth; x++) {
      for (uint32_t i = 0; i < 4; i++) {
        tmpColor = 0.0f;
        for (uint32_t row = 0; row < this->imageScale; row++) {
          for (uint32_t col = 0; col < this->imageScale; col++) {
            tmpColor += (float)(this->privateRenderimage[((((this->imageScale * y) + row) * (this->imageWidth)) + x) * this->imageScale * 4 + 4 * col + i]);
          }
        }
        tmpColor /= (float)(this->antialiasingValue); // 1, 4 or 16
        tmpColor = (tmpColor < 255.0f) ? (tmpColor) : (255.0f);

        this->image[((y * (this->imageWidth)) + x) * 4 + i] = (unsigned char)tmpColor;
      }
    }
  }

  delete[](this->privateRenderimage);

  return VRL_OK;
}

void VRL::renderDrawBox() {
  uint32_t cubePointNumber[8] = { 18, 21, 3, 0, 30, 33, 15, 12 };
  float points2D[8 * 2];
  float tmpVec[4];

  // MVP matrix
  float *mvpMatrix = renderGetMVPMatrix();

  // Create box for volume
  float *pointsTrianglesForBox = renderCreateBoxForVolume();

  float x;
  float y;
  for (uint32_t i = 0; i < 8; i++) {
    tmpVec[0] = pointsTrianglesForBox[cubePointNumber[i] + 0];
    tmpVec[1] = pointsTrianglesForBox[cubePointNumber[i] + 1];
    tmpVec[2] = pointsTrianglesForBox[cubePointNumber[i] + 2];
    tmpVec[3] = 1.0f;

    matrix4fMultVec4f(mvpMatrix, tmpVec, tmpVec);

    // set point
    x = tmpVec[0] / tmpVec[3];
    y = tmpVec[1] / tmpVec[3];

    x = (x + 1.0f) * (float)(this->imageWidth * this->imageScale) / 2.0f;
    y = (float)(this->imageHeight * this->imageScale) - (y + 1.0f) * (float)(this->imageHeight * this->imageScale) / 2.0f;

    x = roundf(x);
    y = roundf(y);

    points2D[i*2] = x;
    points2D[i*2 + 1] = y;
  }

  unsigned char color[4] = {255, 255, 255, 255};

  brezenhem(points2D[0], points2D[1], points2D[2], points2D[3], color);
  brezenhem(points2D[2], points2D[3], points2D[4], points2D[5], color);
  brezenhem(points2D[4], points2D[5], points2D[6], points2D[7], color);
  brezenhem(points2D[6], points2D[7], points2D[0], points2D[1], color);

  brezenhem(points2D[8], points2D[9], points2D[10], points2D[11], color);
  brezenhem(points2D[10], points2D[11], points2D[12], points2D[13], color);
  brezenhem(points2D[12], points2D[13], points2D[14], points2D[15], color);
  brezenhem(points2D[14], points2D[15], points2D[8], points2D[9], color);

  brezenhem(points2D[0], points2D[1], points2D[8], points2D[9], color);
  brezenhem(points2D[2], points2D[3], points2D[10], points2D[11], color);
  brezenhem(points2D[4], points2D[5], points2D[12], points2D[13], color);
  brezenhem(points2D[6], points2D[7], points2D[14], points2D[15], color);
}

void VRL::brezenhem(int x0, int y0, int x1, int y1, unsigned char *rgba) { 
  int A, B, sign;
  A = y1 - y0;
  B = x0 - x1;
  
  if (abs(A) > abs(B))
    sign = 1;
  else 
    sign = -1;
    
  int signa, signb;

  if (A < 0)
    signa = -1;
  else
    signa = 1;
    
  if (B < 0)
    signb = -1;  
  else 
    signb = 1;

  int f = 0;
  this->privateRenderimage[(y0 * this->imageWidth + x0) * 4 + 0] = rgba[0];
  this->privateRenderimage[(y0 * this->imageWidth + x0) * 4 + 1] = rgba[1];
  this->privateRenderimage[(y0 * this->imageWidth + x0) * 4 + 2] = rgba[2];
  this->privateRenderimage[(y0 * this->imageWidth + x0) * 4 + 3] = rgba[3];

  int x = x0, y = y0;
  if (sign == -1) {
    do {
      f += A*signa;
      if (f > 0) {
        f -= B*signb;
        y += signa;
      }
      x -= signb;
      this->privateRenderimage[(y * this->imageWidth + x) * 4 + 0] = rgba[0];
      this->privateRenderimage[(y * this->imageWidth + x) * 4 + 1] = rgba[1];
      this->privateRenderimage[(y * this->imageWidth + x) * 4 + 2] = rgba[2];
      this->privateRenderimage[(y * this->imageWidth + x) * 4 + 3] = rgba[3];
    } while (x != x1 || y != y1);
  } else {
    do {
      f += B*signb;
      if (f > 0) {
        f -= A*signa;
        x -= signb;
      }
      y += signa;
      this->privateRenderimage[(y * this->imageWidth + x) * 4 + 0] = rgba[0];
      this->privateRenderimage[(y * this->imageWidth + x) * 4 + 1] = rgba[1];
      this->privateRenderimage[(y * this->imageWidth + x) * 4 + 2] = rgba[2];
      this->privateRenderimage[(y * this->imageWidth + x) * 4 + 3] = rgba[3];
    } while (x != x1 || y != y1);
  }
}

void VRL::setEnableDrawBox(int32_t en) {
  this->enableDrawBox = en;
}
