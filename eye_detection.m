% Please edit this function only, and submit this Matlab file in a zip file
% along with your PDF report

function [left_x, right_x, left_y, right_y] = eye_detection(img)
% INPUT: RGB image
% OUTPUT: x and y coordinates of left and right eye.
% Please rewrite this function, and submit this file in Moodle (in a zip file with the report). 

left_x = 100;
right_x = 300;
left_y = 50;
right_y = 50;

X= img;
X2 = X;
[y1,x1,z1] = size(X);
res=false;
if(y1>1300 || x1>1300)
    display('resized');
    res = true;
    X = imresize(X,.5);
end

forcrop=X;
X=imcomplement(X);
sup = X;

imagesc(X)
%figure,imshow(rgb2gray(X));
X=rgb2hsv(X);
X=X(:,:,2); 
%figure,imshow(X);
% Threshold image - adaptive threshold
BW = imbinarize(X, 'adaptive', 'Sensitivity', 0.900000, 'ForegroundPolarity', 'bright');

BW = imcomplement(BW);
radius = 4;
decomposition = 0;
se = strel('disk', radius, decomposition);
BW = imopen(BW, se);
BW2= im2bw(BW);
[B,L] = bwboundaries(BW);
%figure, imshow(label2rgb(L, @jet, [.5 .5 .5]))



stats = regionprops(L,'Centroid',...
    'MajorAxisLength','MinorAxisLength','Area','BoundingBox','Orientation')


pos= [];
circles = []; 
rad= [];
met=[];
bbnum=[];
for k = 1:size(B)
    imarea = y1*x1;
   if (stats(k).Orientation > -80) && (stats(k).Area > 100)  && (stats(k).Area < (imarea*.40))
      thisBB = stats(k).BoundingBox;
        I2 = imcrop(sup,stats(k).BoundingBox);
I2 = rgb2gray(I2);
rgbImage = cat(3, I2, I2, I2);
I = rgb2hsv(rgbImage);

% Define thresholds for channel 1 based on histogram settings
channel1Min = 0.000;
channel1Max = 1.000;

% Define thresholds for channel 2 based on histogram settings
channel2Min = 0.000;
channel2Max = 1.000;

% Define thresholds for channel 3 based on histogram settings
channel3Min = 0.256;
channel3Max = 1.000;

% Create mask based on chosen histogram thresholds
BW = (I(:,:,1) >= channel1Min ) & (I(:,:,1) <= channel1Max) & ...
    (I(:,:,2) >= channel2Min ) & (I(:,:,2) <= channel2Max) & ...
    (I(:,:,3) >= channel3Min ) & (I(:,:,3) <= channel3Max);

% Initialize output masked image based on input image.
maskedRGBImage = rgbImage;

% Set background pixels where BW is false to zero.
maskedRGBImage(repmat(~BW,[1 1 3])) = 0;
rectangle('Position', [thisBB(1),thisBB(2),thisBB(3),thisBB(4)],...
      'EdgeColor','r','LineWidth',2 );

  
        [centers,radii,metric]=imfindcircles(maskedRGBImage,[10 30],'Method','TwoStage','ObjectPolarity','bright','Sensitivity',.87,'EdgeThreshold',.03);
        [centers2,radii2,metric2]=imfindcircles(maskedRGBImage,[30 80],'Method','TwoStage','ObjectPolarity','bright','Sensitivity',.87,'EdgeThreshold',.03);
        
        if (~isempty(centers2))
        centers= cat(1,centers,centers2);
        radii=cat(1,radii,radii2);
        metric=cat(1,metric,metric2);
        end
        
    
        [m,n]=size(centers);
        if m~=0 && m<6
            count = m;
            circles = cat(1,circles,centers);
            rad = cat(1,rad,radii);
            met =cat(1,met,metric);
            %while count>0
            pos = cat(1,pos,[thisBB(1),thisBB(2),thisBB(3),thisBB(4)]);
            %count = count-1;
            %end
         end
         
            
        
        end
   end   
   
  


    %find top circles)
    center=[];
    radius=[];
    [m n] =size(pos);
    figure,imshow(sup);
    for j=1:m
        if(res == false)
        main= sup;
       
        crop = imcrop(main,[pos(j,1),pos(j,2),pos(j,3),pos(j,4)]);
         [centers,radii,metric]=imfindcircles(crop,[5 20],'Method','TwoStage','ObjectPolarity','bright', 'Sensitivity',.96,'EdgeThreshold',.02 );
        viscircles([pos(j,1) + circles(j,1),(pos(j,2)+ circles(j,2))],rad(j));
        
        [k na] = size(centers);
        
         if(~isempty(centers))
            center = cat(1,center,[pos(j,1) + centers(1,1),(pos(j,2)+ centers(1,2))]);
            radius= cat(1,radius,radii(1,1));
            rectangle('Position', [pos(j,1),pos(j,2),pos(j,3),pos(j,4)],...
      'EdgeColor','r','LineWidth',2 );   
         end
        
        else
            main= sup;
            rectangle('Position', [pos(j,1),pos(j,2),pos(j,3),pos(j,4)],...
      'EdgeColor','r','LineWidth',2);
       
        crop = imcrop(main,[pos(j,1),pos(j,2),pos(j,3),pos(j,4)]);
         [centers,radii,metric]=imfindcircles(crop,[5 20],'Method','TwoStage','ObjectPolarity','bright', 'Sensitivity',.96,'EdgeThreshold',.02 );
        viscircles([pos(j,1) + circles(j,1),(pos(j,2)+ circles(j,2))],rad(j));
        
        [k na] = size(centers);
        
         if(~isempty(centers))
            center = cat(1,center,[pos(j,1)/.5 + centers(1,1),(pos(j,2)/.5+ centers(1,2))]);
            radius= cat(1,radius,radii(1,1)/.5);
            rectangle('Position', [pos(j,1),pos(j,2),pos(j,3),pos(j,4)],...
      'EdgeColor','r','LineWidth',2 );   

         end
        end
    end
    figure, imshow(X2)
    viscircles(center,radius);
    
    %figure, imtool(sup);
    
    [m n] = size(center);
    
    max_corr= 0;
    max_sim=0;
    scalar=8;
    if(res == true)
    sup=X2;
    sup= imcomplement(sup);
    scalar=16;
    end
    
    for i= 1:m

    temp = imcrop(sup,[center(i,1)-(radius(i)*2),center(i,2)-(radius(i)*2),(radius(i)*scalar),radius(i)*scalar]);
    
    temp= rgb2gray(temp);
    temp=flip(temp,2);
        for o=i+1:m-1
         
            cr= imcrop(sup,[center(o,1)-(radius(o)*2),center(o,2)-(radius(o)*2),(radius(i)*scalar),radius(i)*scalar]);
            cr= imresize(cr,size(temp));
            cr= rgb2gray(cr);
            
            [optimizer, metric] = imregconfig('multimodal')

            optimizer.InitialRadius = 0.009;
            optimizer.Epsilon = 1.5e-4;
            optimizer.GrowthFactor = 1.01;
            optimizer.MaximumIterations = 300;
            tform = imregtform(cr,temp,'similarity',optimizer,metric);
            cr = imwarp(cr,tform,'OutputView',imref2d(size(temp)));
            C1 = normxcorr2((temp),(cr));
            subplot(1,2,1); imshow(temp); title('Reference Image');
            subplot(1,2,2); imshow(cr);   title('compare Image');
            if(max(C1(:))>max_corr)
                max_corr = max(C1(:));
                left_x=center(i,1);
                left_y=center(i,2);
                
                right_x=center(o,1);
                right_y=center(o,2);
                
            end
            
        end
        
    
    end
    


end