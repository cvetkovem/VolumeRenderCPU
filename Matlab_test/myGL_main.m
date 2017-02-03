function img_out = myGL_main()
  window_width = 512;
  window_height = 512;
  
  img_out = zeros(window_height, window_width);
  
  %%%
  %triangle_points = [128, 300;
  %                   256, 100;
  %                   384, 350];
  %img_out = myGL_rasterization_triangle(img_out, triangle_points, window_width, window_height);
  %return;
  %%%
  
  camera_position_vec3 = [0, 1.0, -2];
  camera_target_vec3 = [0.0, -0.5, 1.0];
  camera_up_vec3 = [0, 1, 0];

  camera_translate_mat4 = myGL_camera_translate_mat4(camera_position_vec3(1), camera_position_vec3(2), camera_position_vec3(3));
  camera_rotate_mat4 = myGL_camera_rotate_mat4(camera_target_vec3, camera_up_vec3);
  
  perspective_matrix_mat4 = myGL_get_perspective_mat4(60, window_width, window_height, 0.1, 1000);
  
  model_matrix_mat4 = [1, 0, 0, 0;
                       0, 1, 0, 0;
                       0, 0, 1, 0;
                       0, 0, 0, 1];
  
  figure(1);
  lHandle = imshow(img_out);

for j = 1:1
                   
  img_out = zeros(window_height, window_width);                 
  
  model_matrix_mat4 = myGL_rotate_mat4(25, [0, 1, 0]);
  %model_matrix_mat4 = myGL_rotate_mat4(6*j, [0, 1, 0]);
  %model_matrix_mat4 = myGL_rotate_mat4(-25, [1, 0, 0]) * model_matrix_mat4;
  
  %camera_translate_mat4 = eye(4,4);
  %camera_rotate_mat4 = eye(4,4);
  
  transformation_mat4 = myGL_get_transformation_mat4(perspective_matrix_mat4, camera_rotate_mat4, camera_translate_mat4, model_matrix_mat4);
  
  %test BEGIN
  %disp(camera_translate_mat4);
  %disp(camera_rotate_mat4);
  %test END
 
  %      ^ Y 
  %      |
  %  *1  |  *2  
  % -----|-----> X
  %  *4  |  *3

  point_front1_vec4 = [-0.5,  0.5, -0.5, 1]';
  point_front2_vec4 = [ 0.5,  0.5, -0.5, 1]';
  point_front3_vec4 = [ 0.5, -0.5, -0.5, 1]';
  point_front4_vec4 = [-0.5, -0.5, -0.5, 1]';
  
  point_back1_vec4 = [-0.5,  0.5, 0.5, 1]';
  point_back2_vec4 = [ 0.5,  0.5, 0.5, 1]';
  point_back3_vec4 = [ 0.5, -0.5, 0.5, 1]';
  point_back4_vec4 = [-0.5, -0.5, 0.5, 1]';
    
  p1_front_out = transformation_mat4 * point_front1_vec4;
  p2_front_out = transformation_mat4 * point_front2_vec4;
  p3_front_out = transformation_mat4 * point_front3_vec4;
  p4_front_out = transformation_mat4 * point_front4_vec4;
  
  p1_back_out = transformation_mat4 * point_back1_vec4;
  p2_back_out = transformation_mat4 * point_back2_vec4;
  p3_back_out = transformation_mat4 * point_back3_vec4;
  p4_back_out = transformation_mat4 * point_back4_vec4;
  
  out = [p1_front_out'; p2_front_out'; p3_front_out'; p4_front_out';
         p1_back_out'; p2_back_out'; p3_back_out'; p4_back_out'];
  
  out_point = [];
  for i=1:8
    x = out(i,1) / out(i,4);
    y = out(i,2) / out(i,4);
    
    x = (x + 1) * window_width / 2;
    y = window_height - (y + 1) * window_height / 2;
    
    x = round(x);
    y = round(y);
    
    disp([x y out(i,4)]);
    
    out_point = [out_point; [x y]];
    
    img_out = set_pixel(img_out, x, y, window_width, window_height);
  end  
  
  %disp(out_point);
  
  img_out = myGL_brezenhem(img_out, out_point(1,1), out_point(1,2), out_point(2,1), out_point(2,2), window_width, window_height);
  img_out = myGL_brezenhem(img_out, out_point(2,1), out_point(2,2), out_point(3,1), out_point(3,2), window_width, window_height);
  img_out = myGL_brezenhem(img_out, out_point(3,1), out_point(3,2), out_point(4,1), out_point(4,2), window_width, window_height);
  img_out = myGL_brezenhem(img_out, out_point(4,1), out_point(4,2), out_point(1,1), out_point(1,2), window_width, window_height);
  
  img_out = myGL_brezenhem(img_out, out_point(5,1), out_point(5,2), out_point(6,1), out_point(6,2), window_width, window_height);
  img_out = myGL_brezenhem(img_out, out_point(6,1), out_point(6,2), out_point(7,1), out_point(7,2), window_width, window_height);
  img_out = myGL_brezenhem(img_out, out_point(7,1), out_point(7,2), out_point(8,1), out_point(8,2), window_width, window_height);
  img_out = myGL_brezenhem(img_out, out_point(8,1), out_point(8,2), out_point(5,1), out_point(5,2), window_width, window_height); 
  
  img_out = myGL_brezenhem(img_out, out_point(1,1), out_point(1,2), out_point(5,1), out_point(5,2), window_width, window_height);
  img_out = myGL_brezenhem(img_out, out_point(2,1), out_point(2,2), out_point(6,1), out_point(6,2), window_width, window_height);
  img_out = myGL_brezenhem(img_out, out_point(3,1), out_point(3,2), out_point(7,1), out_point(7,2), window_width, window_height);
  img_out = myGL_brezenhem(img_out, out_point(4,1), out_point(4,2), out_point(8,1), out_point(8,2), window_width, window_height); 
  
  triangle_points = [out_point(1,1), out_point(1,2);
                     out_point(2,1), out_point(2,2);
                     out_point(3,1), out_point(3,2)];
  %img_out = myGL_rasterization_triangle(img_out, triangle_points, window_width, window_height);
  
  set(lHandle, 'CData', img_out);
  
  %filename = strcat('C:\Users\DJK\Desktop\cube\cub', num2str(j), '.png');
  %imwrite(img_out, filename)
  
  pause(0.1);
end

end

function img = set_pixel(img, x, y, window_width, window_height)
 if(x > 0 && x < window_width + 1 && y > 0 && y < window_height + 1)
   img(y, x) = 1;
 end
end

% glm::mat4 myMatrix = glm::translate(10.0f, 0.0f, 0.0f);
% use: %glm::vec4 transformedVector = myMatrix * myVector;
function translate_mat4 = myGL_translate_mat4(x, y, z)
  translate_mat4 = [1 0 0 x;
                    0 1 0 y;
                    0 0 1 z;
                    0 0 0 1];
end

% glm::mat4 myScalingMatrix = glm::scale(2.0f, 2.0f ,2.0f);
% use: glm::vec4 transformedVector = myScalingMatrix * myVector;
function scale_mat4 = myGL_scale_mat4(x, y ,z)
  scale_mat4 = [x 0 0 0;
                0 y 0 0;
                0 0 z 0;
                0 0 0 1];
end

% glm::rotate( angle_in_degrees, myRotationAxis );
function rotate_mat4 = myGL_rotate_mat4(angle_in_degrees, myRotationAxis_vec3)
  % RotationAngle is in radians
  RotationAngle = angle_in_degrees * (pi / 180);
  qx = myRotationAxis_vec3(1) * sin(RotationAngle / 2);
  qy = myRotationAxis_vec3(2) * sin(RotationAngle / 2);
  qz = myRotationAxis_vec3(3) * sin(RotationAngle / 2);
  qw = cos(RotationAngle / 2);

  rotate_mat4 = [1 - 2*qy*qy - 2*qz*qz, 2*qx*qy - 2*qz*qw,     2*qx*qz + 2*qy*qw,     0;
                 2*qx*qy + 2*qz*qw,     1 - 2*qx*qx - 2*qz*qz, 2*qy*qz - 2*qx*qw,     0;
                 2*qx*qz - 2*qy*qw,     2*qy*qz + 2*qx*qw,     1 - 2*qx*qx - 2*qy*qy, 0;
                 0,                     0,                     0,                     1];
end

% glm::mat4 myModelMatrix = myTranslationMatrix * myRotationMatrix * myScaleMatrix;
% glm::vec4 myTransformedVector = myModelMatrix * myOriginalVector;
function model_matrix_mat4 = myGL_get_model_matrix_mat4(translate_mat4, rotate_mat4, scale_mat4)
  model_matrix_mat4 = translate_mat4 * rotate_mat4 * scale_mat4;
end

% glm::mat4 Projection = glm::perspective(45.0f, 1.0f, 0.1f, 100.0f);
function perspective_matrix_mat4 = myGL_get_perspective_mat4(angle_in_degrees, window_width, window_height, zNear, zFar)
  ar = window_width / window_height;
  zRange = zNear - zFar;
  angle_in_radian = angle_in_degrees * (pi / 180);
  tanHalfFOV = tan(angle_in_radian / 2);
 
  perspective_matrix_mat4 = [1 / (tanHalfFOV * ar), 0,              0,                        0;
                             0,                     1 / tanHalfFOV, 0,                        0;
                             0,                     0,              (-zNear - zFar) / zRange, 2 * zFar * zNear / zRange;
                             0,                     0,              1,                        0];
end

% cross: C = cross(A, B)
% normalize: vec3_normalize = V/sqrt(sum(V.*V))

function camera_translate_mat4 = myGL_camera_translate_mat4(x, y, z)
 camera_translate_mat4 = [1 0 0 -x;
                          0 1 0 -y;
                          0 0 1 -z;
                          0 0 0  1];
end

function camera_rotate_mat4 = myGL_camera_rotate_mat4(m_camera_target_vec3, m_camera_up_vec3)
 N = m_camera_target_vec3;
 N = N/sqrt(sum(N.*N));
 
 U = m_camera_up_vec3;
 U = U/sqrt(sum(U.*U));
 U = cross(U, m_camera_target_vec3);
 
 V = cross(N, U);
 
 camera_rotate_mat4 = [U(1), U(2), U(3), 0;
                       V(1), V(2), V(3), 0;
                       N(1), N(2), N(3), 0;
                       0,    0,    0,    1];
end

function transformation_mat4 = myGL_get_transformation_mat4(perspective_matrix_mat4, camera_rotate_mat4, camera_translate_mat4, model_matrix_mat4)
  transformation_mat4 = perspective_matrix_mat4 * camera_rotate_mat4 * camera_translate_mat4 * model_matrix_mat4;
end               

function img = myGL_brezenhem(img, x0, y0, x1, y1, window_width, window_height)
  A = y1 - y0;
  B = x0 - x1;
  
  if (abs(A) > abs(B)) 
    sign = 1;
  else
    sign = -1;
  end
  
  if (A < 0)
    signa = -1;
  else
    signa = 1;
  end
  
  if (B < 0)
    signb = -1;
  else
    signb = 1;
  end
  
  f = 0;
  %img(y0, x0) = 1; % set point
  img = set_pixel(img, x0, y0, window_width, window_height);
  
  x = x0;
  y = y0;
  
  if (sign == -1)
    while(1)
      f = f + A*signa;
      if (f > 0)
        f = f - B*signb;
        y = y + signa;
      end
      x = x - signb;
      %img(y, x) = 1; % set point
      img = set_pixel(img, x, y, window_width, window_height);
      if (x ~= x1 || y ~= y1)
          continue;
      else
          break;
      end
    end
  else
    while(1)
      f = f + B*signb;
      if (f > 0)
        f = f - A*signa;
        x = x - signb;
      end
      y = y + signa;
      %img(y, x) = 1; % set point
      img = set_pixel(img, x, y, window_width, window_height);
     if (x ~= x1 || y ~= y1)
         continue;
     else
         break;
     end
    end
  end
end

%% Rasterization
function out =  myGL_rast_sign(p1x, p1y, p2x, p2y, p3x, p3y)
  out =  (p1x - p3x) * (p2y - p3y) - (p2x - p3x) * (p1y - p3y);
end

function out = myGL_rast_point_in_triangle(ptx, pty, v1x, v1y, v2x, v2y, v3x, v3y)
  b1 = myGL_rast_sign(ptx, pty, v1x, v1y, v2x, v2y) < 0;
  b2 = myGL_rast_sign(ptx, pty, v2x, v2y, v3x, v3y) < 0;
  b3 = myGL_rast_sign(ptx, pty, v3x, v3y, v1x, v1y) < 0;

  out =  ((b1 == b2) && (b2 == b3));
end

function img = myGL_rasterization_triangle(img, triangle_points, window_width, window_height)
  %find box for triangle
  min_x = triangle_points(1,1);
  max_x = triangle_points(1,1);
  min_y = triangle_points(1,2);
  max_y = triangle_points(1,2);
  
  for i=2:3
    if(min_x > triangle_points(i,1))
      min_x = triangle_points(i,1);
    end
    
    if(max_x < triangle_points(i,1))
      max_x = triangle_points(i,1);
    end
    
    if(min_y > triangle_points(i,2))
      min_y = triangle_points(i,2);
    end
    
    if(max_y < triangle_points(i,2))
      max_y = triangle_points(i,2);
    end
  end

  for i=min_y:max_y
    for j=min_x:max_x
      if(i > 0 && i < window_height + 1 && j > 0 && j < window_width + 1)
          if(myGL_rast_point_in_triangle(j, i, triangle_points(1,1), triangle_points(1,2), triangle_points(2,1), triangle_points(2,2), triangle_points(3,1), triangle_points(3,2)))
            img(i, j) = 1;
          end
      end
    end
  end
end 
