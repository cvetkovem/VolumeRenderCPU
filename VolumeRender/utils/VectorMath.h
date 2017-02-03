#ifndef VECTOR_MATH_H
#define	VECTOR_MATH_H

#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define ToRadian(x) ((x) * M_PI / 180.0f)
#define ToDegree(x) ((x) * 180.0f / M_PI)

inline void copyVec3f(float *out, float *in) {
  out[0] = in[0];
  out[1] = in[1];
  out[2] = in[2];
}

inline void setVec3f(float *out, float x, float y, float z) {
  out[0] = x;
  out[1] = y;
  out[2] = z;
}

inline void addVec3f(float *a, float *b, float *out) {
  out[0] = a[0] + b[0];
  out[1] = a[1] + b[1];
  out[2] = a[2] + b[2];
}

inline void subVec3f(float *a, float *b, float *out) {
  out[0] = a[0] - b[0];
  out[1] = a[1] - b[1];
  out[2] = a[2] - b[2];
}

inline void multToConstVec3f(float *a, float c) {
  a[0] *= c;
  a[1] *= c;
  a[2] *= c;
}

inline void crossVec3f(float *a, float *b, float *out) {
  const float tmp0 = a[1] * b[2] - a[2] * b[1];
  const float tmp1 = a[2] * b[0] - a[0] * b[2];
  const float tmp2 = a[0] * b[1] - a[1] * b[0];

  out[0] = tmp0;
  out[1] = tmp1;
  out[2] = tmp2;
}

inline void normalizeVec3f(float *a) {
  const float length = sqrtf(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);

  a[0] /= length;
  a[1] /= length;
  a[2] /= length;
}

inline void dotVec3f(float *a, float *b, float *out) {
  *out = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

inline void setVec4f(float *out, float x, float y, float z, float w) {
  out[0] = x;
  out[1] = y;
  out[2] = z;
  out[3] = w;
}

inline void multMatrix4fx4f(float *lMatrix, float *rMatrix, float *out) {
  float tmp[16];

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      tmp[(i * 4) + j] =
        lMatrix[(i * 4) + 0] * rMatrix[(0 * 4) + j] +
        lMatrix[(i * 4) + 1] * rMatrix[(1 * 4) + j] +
        lMatrix[(i * 4) + 2] * rMatrix[(2 * 4) + j] +
        lMatrix[(i * 4) + 3] * rMatrix[(3 * 4) + j];
    }
  }

  for (int i = 0; i < 16; i++) {
    out[i] = tmp[i];
  }
}

inline void matrix4fMultVec4f(float *matrix, float *vector, float *outVector) {
  float tmp[4];

  for (int i = 0; i < 4; i++) {
    tmp[i] =
      matrix[i * 4 + 0] * vector[0] +
      matrix[i * 4 + 1] * vector[1] +
      matrix[i * 4 + 2] * vector[2] +
      matrix[i * 4 + 3] * vector[3];
  }

  outVector[0] = tmp[0];
  outVector[1] = tmp[1];
  outVector[2] = tmp[2];
  outVector[3] = tmp[3];
}

#endif	/* VECTOR_MATH_H */
