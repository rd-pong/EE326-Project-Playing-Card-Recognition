I = im2double(rgb2gray(imread('training/Qspades.jpg')));
BW = im2bw(I, 0.65);
figure;
imshowpair(I, BW, 'montage');
title('Synthetic Image & Binary Image');

%% Remove small size and only keep white card
% label the complement of the Binary Image
imLabel = bwlabel(1-BW);
% Find the card shape by finding the largest connected region
stats = regionprops(imLabel,'centroid', 'Area');
[b,index]=sort([stats.Area],'descend');
if length(stats)<2
    BW2=imLabel;
else
    BW2=ismember(imLabel,index(1:2));
end
% Reverse complement
BW2 = 1-BW2;
imshowpair(BW, BW2, 'montage');
title('Binary Image & Edited Binary Image');

%% Find lines
findingLine = edge(BW2, 'canny');
[H,T,R] = hough(findingLine);
P  = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));
% change parameter if different data set 
lines = houghlines(findingLine,T,R,P,'FillGap',2000,'MinLength',400);
figure, imshow(I), hold on
max_len = 0;
for k = 1:length(lines)
   xy = [lines(k).point1; lines(k).point2];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

   % Plot beginnings and ends of lines
   plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
   plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');

   % Determine the endpoints of the longest line segment
   len = norm(lines(k).point1 - lines(k).point2);
   if ( len > max_len)
      max_len = len;
      xy_long = xy;
   end
end

%% Find Line Equation

%% Find Centroid
s = regionprops(max, I, {'Centroid','WeightedCentroid'});
imshow(I)
title('Weighted (red) and Unweighted (blue) Centroids'); 
hold on
numObj = numel(s);
for k = 1 : numObj
    plot(s(k).WeightedCentroid(1), s(k).WeightedCentroid(2), 'r*');
    plot(s(k).Centroid(1), s(k).Centroid(2), 'o');
end
