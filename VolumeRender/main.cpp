#include <stdio.h>
#include <ctime>

#include "VolumeRenderLibrary.h"

typedef struct _VolumeLoad {
  uint32_t volumeWidth;
  uint32_t volumeHeight;
  uint32_t volumeNumber;

  int16_t *volume;
} VolumeLoad;

int loadVolume(VolumeLoad *data);
void setLUTPoints(VRL *vrl);
void saveImageToFile(unsigned char *image, uint32_t imgWidth, uint32_t imgHeight);
void saveAsBMP(char *filename, unsigned char *image, uint32_t imgWidth, uint32_t imgHeight);

int main() {
  int error = 0;
  VolumeLoad sVolumeData;

  error = loadVolume(&sVolumeData);
  if (error == -1) {
      printf("ERROR: Could not load file with dicom data!\n");
      return 0;
  }

  printf("x:%u y:%u z:%u\n", sVolumeData.volumeWidth, sVolumeData.volumeHeight, sVolumeData.volumeNumber);

  VRL *vrl = new VRL();
  vrl->setVolume(sVolumeData.volume, sVolumeData.volumeWidth, sVolumeData.volumeHeight, sVolumeData.volumeNumber, -2000, 2000);
  delete[] sVolumeData.volume;

  printf("Volume copy complited!\n");

  // Set LUT tabel
  setLUTPoints(vrl);

  // Set volume object in space
  vrl->resetRotationVolume();
  vrl->resetTranslateVolume();
  vrl->setTranslateVolume(0.0f, 0.0f, 0.0f);
  float angle = -90.0f;
  float axis[3] = {0.0f, 0.0f, 1.0f};
  vrl->setRotateVolume(angle, axis);

  angle = 155.0f + 180.0f;
  axis[0] = 0.0f; axis[1] = 1.0f; axis[2] = 0.0f;
  vrl->addRotateVolume(angle, axis);

  angle = -35.0f;
  axis[0] = 1.0f; axis[1] = 0.0f; axis[2] = 0.0f;
  vrl->addRotateVolume(angle, axis);

  // Set light
  float direction[4] = {0.0f, -1.0f, 0.0f, 0.0f};
  vrlColor lightColor;
  lightColor.r = 0.98f;
  lightColor.g = 0.98f;
  lightColor.b = 0.98f;
  lightColor.a = 1.0f;
  lightColor.ambient   = 0.0f;
  lightColor.diffuse   = 1.0f;
  lightColor.specular  = 1.0f;
  lightColor.emission  = 0.0f;
  lightColor.shininess = 0.0f;
  vrl->setLight(direction, &lightColor);

  // Set camera
  float cameraPosition[3] = { 0.0f,  0.75f, -1.5f };
  float cameraTarget[3]   = { 0.0f, -0.5f,   1.0f };
  float cameraUp[3]       = { 0.0f,  1.0f,   0.0f };
  float zNear = 0.1f;
  vrl->setCameraConfigure(cameraPosition, cameraTarget, cameraUp, zNear);

  // Set quality / perfomance
  vrl->enableAntialiasing(VRL_ANTIALIASING_X4);
  vrl->setImageQuality(VRL_IMG_QUA_REALISTIC);

  // Set output image
  uint32_t imgWidth = 512;
  uint32_t imgHeight = 512;
  vrl->allocImage(imgWidth, imgHeight);

  // Set volume box draw
  vrl->setEnableDrawBox(0);

  // Render
  printf("Render start.\n");
  unsigned int startTime = clock();
  vrl->render();
  printf("Render stop.\nTime: %5.2f sec\n", (float)((clock() - startTime) / 1000.0f));
  
  unsigned char *image = vrl->getImage();

  //saveImageToFile(image, imgWidth, imgHeight);
  saveAsBMP("img.bmp", image, imgWidth, imgHeight);

  vrl->~VRL();

  return 0;
}

int loadVolume(VolumeLoad *data) {
  FILE *fp;

  fopen_s(&fp, "dicomMatrix.bin", "rb");

  if (fp) {
      fread(&(data->volumeWidth), sizeof(int), 1, fp);
      fread(&(data->volumeHeight), sizeof(int), 1, fp);
      fread(&(data->volumeNumber), sizeof(int), 1, fp);

      data->volume = new int16_t[data->volumeWidth * data->volumeHeight * data->volumeNumber];
      fread(data->volume, sizeof(int16_t), data->volumeWidth * data->volumeHeight * data->volumeNumber, fp);

      fclose(fp);

      // test goto HU
      for (uint32_t i = 0; i < data->volumeWidth * data->volumeHeight * data->volumeNumber; i++) {
        data->volume[i] = data->volume[i] - 1024;
      }

      return 0;
  }
  return -1;
}

void setLUTPoints(VRL *vrl) {
/*
  name="Tissue"
  <color density = "25.684" opacity = "0"     r = "0"   g = "0"   b = "0"   ambient = "0.25" diffuse = "0.8" specular = "0.8" / >
  <color density = "58.001" opacity = "0.195" r = "255" g = "0"   b = "0"   ambient = "0.25" diffuse = "0.8" specular = "0.8" / >
  <color density = "85.83"  opacity = "0.51"  r = "255" g = "255" b = "0"   ambient = "0.25" diffuse = "0.8" specular = "0.8" / >
  <color density = "100"    opacity = "1"     r = "255" g = "255" b = "255" ambient = "0.25" diffuse = "0.8" specular = "0.8" / >
  <shiny value="15"/>
*/

  vrl->clearLUT();

  vrlColor color;
  color.r = 0.0f;
  color.g = 0.0f;
  color.b = 0.0f;
  color.a = 0.0f;
  color.ambient   = 0.25f;
  color.diffuse   = 0.8f;
  color.specular  = 0.8f;
  color.emission  = 0.0f;
  color.shininess = 15.0f; // 15

  //w = 400 c = 50; min = -150 max = 250;
  vrl->addLUTPoint(/*-758.0f*/ -47.264f - 800.0f, &color); // density = "25.684" %

  color.r = 0.78f;
  color.a = 0.195f;
  vrl->addLUTPoint(/*-630.0f*/ 82.004f - 800.0f, &color);  // density = "58.001" %

  color.g = .68f;
  color.a = 0.51f;
  vrl->addLUTPoint(/*-470.0f*/ 193.32f - 800.0f, &color);  // density = "85.83"  %

  color.b = 0.6f;
  color.a = 1.0f;
  vrl->addLUTPoint(/*-462.0f*/ 250.0f - 800.0f, &color);  // density = "100.0"  %

  vrl->interpolateLUT();
}

void saveImageToFile(unsigned char *image, uint32_t imgWidth, uint32_t imgHeight) {
  FILE *fp;

  fopen_s(&fp, "SAVE_IMG.TXT", "w");

  if (fp) {
    for (uint32_t i = 0; i < imgWidth * imgHeight * 4; i+=4) {
      fprintf(fp, "%u %u %u %u\n", image[i + 0], image[i + 1], image[i + 2], image[i + 3]);
    }

    fclose(fp);
  }
}

void saveAsBMP(char *filename, unsigned char *image, uint32_t imgWidth, uint32_t imgHeight) {
  unsigned int headers[13];
  FILE *outfile;
  int extrabytes;
  int paddedsize;
  uint32_t x; int y; int n;
  int red, green, blue;

  extrabytes = 4 - ((imgWidth * 3) % 4);            // How many bytes of padding to add to each
                                                    // horizontal line - the size of which must
                                                    // be a multiple of 4 bytes.
  if (extrabytes == 4)
    extrabytes = 0;

  paddedsize = ((imgWidth * 3) + extrabytes) * imgHeight;

  // Note that the "BM" identifier in bytes 0 and 1 is NOT included in these "headers".                   
  headers[0]  = paddedsize + 54;      // bfSize (whole file size)
  headers[1]  = 0;                    // bfReserved (both)
  headers[2]  = 54;                   // bfOffbits
  headers[3]  = 40;                   // biSize
  headers[4]  = imgWidth;  // biWidth
  headers[5]  = imgHeight; // biHeight

  // Would have biPlanes and biBitCount in position 6, but they're shorts.
  // It's easier to write them out separately (see below) than pretend
  // they're a single int, especially with endian issues

  headers[7]  = 0;                    // biCompression
  headers[8]  = paddedsize;           // biSizeImage
  headers[9]  = 0;                    // biXPelsPerMeter
  headers[10] = 0;                    // biYPelsPerMeter
  headers[11] = 0;                    // biClrUsed
  headers[12] = 0;                    // biClrImportant

  fopen_s(&outfile, filename, "wb");

  fprintf(outfile, "BM");

  for (n = 0; n <= 5; n++) {
   fprintf(outfile, "%c", headers[n] & 0x000000FF);
   fprintf(outfile, "%c", (headers[n] & 0x0000FF00) >> 8);
   fprintf(outfile, "%c", (headers[n] & 0x00FF0000) >> 16);
   fprintf(outfile, "%c", (headers[n] & (unsigned int) 0xFF000000) >> 24);
  }

  // These next 4 characters are for the biPlanes and biBitCount fields.

  fprintf(outfile, "%c", 1);
  fprintf(outfile, "%c", 0);
  fprintf(outfile, "%c", 24);
  fprintf(outfile, "%c", 0);

  for (n = 7; n <= 12; n++) {
   fprintf(outfile, "%c", headers[n] & 0x000000FF);
   fprintf(outfile, "%c", (headers[n] & 0x0000FF00) >> 8);
   fprintf(outfile, "%c", (headers[n] & 0x00FF0000) >> 16);
   fprintf(outfile, "%c", (headers[n] & (unsigned int) 0xFF000000) >> 24);
  }

  // Headers done, now write the data
  for (y = imgHeight - 1; y >= 0; y--) {     // BMP image format is written from bottom to top
    for (x = 0; x <= imgWidth - 1; x++) {
      red   = image[(y * imgWidth + x) * 4 + 0];
      green = image[(y * imgWidth + x) * 4 + 1];
      blue  = image[(y * imgWidth + x) * 4 + 2];
      
      if (red > 255) red = 255; if (red < 0) red = 0;
      if (green > 255) green = 255; if (green < 0) green = 0;
      if (blue > 255) blue = 255; if (blue < 0) blue = 0;
      
      // Also, it's written in (b,g,r) format

      fprintf(outfile, "%c", blue);
      fprintf(outfile, "%c", green);
      fprintf(outfile, "%c", red);
    }
  
    if (extrabytes) {      // See above - BMP lines must be of lengths divisible by 4.
      for (n = 1; n <= extrabytes; n++) {
        fprintf(outfile, "%c", 0);
      }
    }
  }

  fclose(outfile);
}
