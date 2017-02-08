function createDicomMatrix(folder)
  [ImagesMatrix, PixelSpacing, SpacingBetweenSlices, SliceThickness, ~] = load_dicom_series(folder);

  % get dicom_width, dicom_height
  dicomWidth = length(ImagesMatrix(:,1,1));
  dicomHeight = length(ImagesMatrix(1,:,1));
  numberOfImages = length(ImagesMatrix(1,1,:));
  
  xStep_mm = PixelSpacing(1);
  yStep_mm = PixelSpacing(2);
  
  if(SpacingBetweenSlices == 0)
    zStep_mm = SliceThickness;
  else
    zStep_mm = SpacingBetweenSlices;
  end
  
  % I suppose x is always equal to y
  dicomBigMatrix = interp3(1:dicomWidth, 1:dicomHeight, 1:numberOfImages, ImagesMatrix, 1:dicomWidth, 1:dicomHeight, [1:1/(zStep_mm / xStep_mm):numberOfImages]', 'cubic');
  
  newNumberOfImages = length(dicomBigMatrix(1,1,:));
  
  outFile = fopen('dicomMatrix.bin', 'wb');

  fwrite(outFile, dicomWidth, 'uint32');
  fwrite(outFile, dicomHeight, 'uint32');
  fwrite(outFile, newNumberOfImages, 'uint32');
  
  for z=1:newNumberOfImages
    disp(z);
    for y=1:dicomHeight
        fwrite(outFile, dicomBigMatrix(:,y,z), 'int16');
    end
  end
  
  fclose(outFile);
end

%LOAD
function [ImagesMatrix, PixelSpacing, SpacingBetweenSlices, SliceThickness, SeriesDescription] = load_dicom_series(path)
  %get files (dicom folder)
  files = dir(path);
  files = files(3:end); % delete '.' and '..'
  
  %get number of files
  number_of_files = length(files(:));

  %get dicom info
  info = dicominfo(sprintf('%s/%s', path, files(1).name));  
  
  %get .SeriesDescription
  try
    SeriesDescription = info.SeriesDescription;
  catch
    SeriesDescription = 'NONE';
  end
  
  %get dicom_width
  dicom_width = info.Width;
  %get dicom_height
  dicom_height = info.Height;

  %get .PixelSpacing
  PixelSpacing = info.PixelSpacing;
  
  %get .SpacingBetweenSlices
  try
    SpacingBetweenSlices = info.SpacingBetweenSlices;
  catch
    SpacingBetweenSlices = 0;
  end
  
  %get .SliceThickness
  try
    SliceThickness = info.SliceThickness;
  catch
    SliceThickness = 0;
  end

  %set matrix images
  ImagesMatrix = zeros(dicom_width, dicom_height, number_of_files);
  for i = 1:number_of_files
    dicom_filename = sprintf('%s/%s', path, files(i).name);
    info = dicominfo(dicom_filename);
    ImagesMatrix(:, :, info.InstanceNumber) = dicomread(sprintf('%s/%s', path, files(i).name));
  end;
end
