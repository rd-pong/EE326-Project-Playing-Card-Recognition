clear;
clc;
I = im2double(rgb2gray(imread('training/4spades.jpg')));
% test.jpg/8hearts.jpg(3edges)/10hearts.jpg(con-clockwise)/Ahearts.jpg(clockwise)
detectedCards = [];
rankTemplates = {...
    imread('templates/1.jpg') 1; ...
    imread('templates/2.jpg') 2; ...
    imread('templates/3.jpg') 3; ...
    imread('templates/4.jpg') 4; ...
    imread('templates/5.jpg') 5; ...
    imread('templates/6.jpg') 6; ...
    imread('templates/7.jpg') 7; ...
    imread('templates/8.jpg') 8; ...
    imread('templates/9.jpg') 9; ...
    imread('templates/10.jpg') 10; ...
    imread('templates/J.jpg') 11; ...
    imread('templates/Q.jpg') 12; ...
    imread('templates/K.jpg') 13; ...
    };

suitTemplates = { ...
    imread('templates/club.jpg') 1; ...
    imread('templates/heart.jpg') 2; ...
    imread('templates/spade.jpg') 3; ...
    imread('templates/diamond.jpg') 4; ...
    };

BW = imbinarize(I, 0.66);
BW = ~bwareaopen(~BW, 100);
BW = bwareaopen(BW, 1000);
% figure;
% imshowpair(I, BW, 'montage');
% title('Synthetic Image & Binary Image');

%% Remove small size and only keep white card
% label the complement of the Binary Image
imLabeled = bwlabel(BW); % Judge carefully if it's BW or 1-BW
% Find the card shape by finding the largest connected region
stats = regionprops(imLabeled,'Centroid', 'Area', 'Image', 'FilledImage', 'BoundingBox');
[b,index]=sort([stats.Area],'descend');
if length(stats)<1
    BW2=imLabeled;
else
    BW2=ismember(imLabeled,index(1:1));
end
% Reverse complement
% BW2 = 1-BW2;

% figure;
% imshowpair(BW, BW2, 'montage');
% title('Binary Image & Edited Binary Image');

%% Find lines
findingLine = edge(BW2, 'canny');
[H,T,R] = hough(findingLine);
P  = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));
% change parameter if different data set
% lines = houghlines(findingLine,T,R,P,'FillGap',3000,'MinLength',200);
lines = houghlines(findingLine,T,R,P,'FillGap',3000);

%plotEdges(I, lines);

%% Arrange Edges to be adjacent in the order of 1234
% line1-line2-line3-line4 are adjacent edges
parallel_threshold = 10;
line1 = lines(1);
line2 = [];
line3 = [];
line4 = [];
if abs(abs(lines(2).theta) - abs(line1.theta)) < parallel_threshold
    line2 = lines(3);
    line3 = lines(2);
    line4 = lines(4);
elseif abs(abs(lines(3).theta) - abs(line1.theta)) < parallel_threshold
    line2 = lines(2);
    line3 = lines(3);
    line4 = lines(4);
else
    line2 = lines(2);
    line3 = lines(4);
    line4 = lines(3);
end
lines = [line1; line2; line3; line4];
% plotEdges(I, lines);

%% Get line equation -> find intersection -> find corners
lineEquations = [];
for j = 1:length(lines)
    line = lines(j);
    x1 = line.point1(1);
    y1 = line.point1(2);
    x2 = line.point2(1);
    y2 = line.point2(2);
    [k, b] = getLineEquation(x1, x2, y1, y2);
    lineEquations = [lineEquations, struct('k', k, 'b', b, 'x', x1, 'y', y1)];
end

corners = [];
corners(1,:) = getIntersection(lineEquations(1), lineEquations(2));
corners(2,:) = getIntersection(lineEquations(2), lineEquations(3));
corners(3,:) = getIntersection(lineEquations(3), lineEquations(4));
corners(4,:) = getIntersection(lineEquations(1), lineEquations(4));

% plot corners
% figure, imshow(BW), title('UnArranged Corners'), hold on
% [a,b] = size(corners);
% for zz = 1:a
%     plot(corners(zz,1), corners(zz,2), 'x','LineWidth',2,'Color','yellow')
%     text(corners(zz,1), corners(zz,2), num2str(zz), 'FontSize', 20, 'Color', 'red');
% end

%% Arrange Corners
% % find centroid
% edgeImage = edge(BW, 'prewitt');
% [rows,cols] = size(edgeImage);
% x = ones(rows,1)*[1:cols];
% y = [1:rows]'*ones(1,cols);
% area = sum(sum(edgeImage));
% meanx = sum(sum(edgeImage.*x))/area;
% meany = sum(sum(edgeImage.*y))/area;
% meanPoint = [meanx, meany];
% % figure, imshow(edgeImage), hold on;
% % plot(meanx,meany,'r+'); %十字标出重心位置
%
% distanceToCentroid(1) = norm(corners(1, :) - meanPoint);
% distanceToCentroid(2) = norm(corners(2, :) - meanPoint);
% distanceToCentroid(3) = norm(corners(3, :) - meanPoint);
% distanceToCentroid(4) = norm(corners(4, :) - meanPoint);
% while min(distanceToCentroid) ~= distanceToCentroid(1)
%     corners = circshift(corners, 1);
%     distanceToCentroid(1) = norm(corners(1, :) - meanPoint);
%     distanceToCentroid(2) = norm(corners(2, :) - meanPoint);
%     distanceToCentroid(3) = norm(corners(3, :) - meanPoint);
%     distanceToCentroid(4) = norm(corners(4, :) - meanPoint);
% end
%
% %plot corners
% figure, imshow(BW), title('Arranged Corners(the shortest edge to be 1-2(counterclock) or 1-4(clockwise))'), hold on
% [a,b] = size(corners);
% for zz = 1:a
%     plot(corners(zz,1), corners(zz,2), 'x','LineWidth',2,'Color','yellow')
%     text(corners(zz,1), corners(zz,2), num2str(zz), 'FontSize', 20, 'Color', 'red');
% end

%% Projective Transform
% make sure line 1(1-4) is the smallest edge
while norm(corners(1, :)-corners(4, :)) > norm(corners(1, :)-corners(2, :))
    % if the corners are rotated counterclockwise
    newCorners = circshift(corners, 1);
    corners = newCorners;
end
% garanteen corners 1-4 rotate clockwise
if norm(corners(1, :)-corners(4, :)) < norm(corners(1, :)-corners(2, :))
    % if the corners are rotated counterclockwise
    newCorners = [corners(1, :); corners(4, :); corners(3, :); corners(2, :)];
    corners = newCorners;
end
card = UprightCard(I, corners);
figure, imshow(card)

%% Crop to Corners(4 corners)
% Crop the image to 2 corners
croppedCorners = CropToCorner(card);
imshow(croppedCorners);
% Detect edges
cornerEdges = edge(croppedCorners, 'prewitt');
% Count edge pixels.
sum1 = sum(cornerEdges(:));

% Arrange corners in couter direction of the above crop image to 2 corners

%% Detection
% Detect edges in the real upper left corner.
% 注意：这里是新的cornerEdges
cornerEdges = edge(croppedCorners, 'canny');

% FFT with the right size.
cornerDft = fft2(cornerEdges, size(rankTemplates{1,1},1), size(rankTemplates{1,1},2));

% Find rankTemplatesConj & suitTemplatesConj
dilateSize = 2;
cornerSize = [105, 67];
for i = 1:length(rankTemplates)
    % Dilate the template
    rankTemplates{i,1} = imdilate(rankTemplates{i,1}, ...
        strel('disk', dilateSize));
    
    % Calculate the FFT. Phase correlation requires removing the
    % mean from the template.
    meanTemplate = mean2(rankTemplates{i,1});
    rankTemplates{i,1} = fft2(rankTemplates{i,1} - ...
        meanTemplate, cornerSize(1), cornerSize(2));
    
    % Conjugate of the FFT.
    currentRankTemplatesConj{i,1} = conj(rankTemplates{i,1});
end
for i = 1:length(suitTemplates)
    % Dilate the template
    suitTemplates{i,1} = imdilate(suitTemplates{i,1}, ...
        strel('disk', dilateSize));
    
    % Calculate the FFT, remove the mean, get the complex conjugate.
    meanTemplate = mean2(suitTemplates{i,1});
    suitTemplates{i,1} = fft2(suitTemplates{i,1} - ...
        meanTemplate, cornerSize(1), cornerSize(2));
    currentSuitTemplatesConj{i,1} = conj(suitTemplates{i,1});
end
rankTemplatesConj = currentRankTemplatesConj;
suitTemplatesConj = currentSuitTemplatesConj;

% if the corner DFT is 0, then our phase correlator will
% have problems (divide by 0). Just set 0s to very small values.
cornerDft(cornerDft == 0) = 0.00001;
bestMatchRank = RunTemplateMatching(cornerDft, rankTemplates(:,1), rankTemplatesConj(:,1));
bestMatchSuit = RunTemplateMatching(cornerDft, suitTemplates(:,1), suitTemplatesConj(:,1));

fprintf('Identified rank: %d suit: %d\n', rankTemplates{bestMatchRank, 2}, suitTemplates{bestMatchSuit, 2});

switch suitTemplates{bestMatchSuit, 2}
    case 1
        display('club');
    case 2
        display('heart');
    case 3
        display('spade');
    case 4
        display('diamond');
    otherwise
        display('error')
end
% If we detect the same card twice...
% Well. Something is wrong. Let's throw the result out.
newRow = [rankTemplates{bestMatchRank, 2} ...
    suitTemplates{bestMatchSuit, 2}];
detectedCards = [detectedCards; newRow];
numValidRegions = 1;

%% Functions
function [k, b] = getLineEquation(x1, x2, y1, y2)
if abs(x1 - x2) < 1e-6
    k = -1;
    b = -1;
else
    kb = [x1 1; x2 1]\[y1;y2];
    k = kb(1);
    b = kb(2);
end
end
%%-------------------
function z = plotEdges(Image, lines)
figure, imshow(Image), hold on
max_len = 0;
for k = 1:length(lines)
    xy = [lines(k).point1; lines(k).point2];
    plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
    % Plot beginnings and ends of lines
    text(xy(1,1),xy(1,2),num2str(k), 'FontSize', 20, 'Color', 'red');
    plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
    plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
    % Determine the endpoints of the longest line segment
    len = norm(lines(k).point1 - lines(k).point2);
    if ( len > max_len)
        max_len = len;
        xy_long = xy;
    end
end
end
%----------------
function [corner] = getIntersection(line1, line2)
x = 1;
y = 1;
% error
if line1.b == -1 && line2.b == -1
    x = 1
end
% line1 is vertical
if line1.b == -1
    x = line1.x;
    y = line2.k * x + line2.b;
    % line2 is vertical
elseif line2.b == -1
    x = line2.x;
    y = line1.k * x + line1.b;
    % Regular case
else
    x = (line2.b - line1.b)/(line1.k-line2.k);
    y = line1.k * x +line1.b;
end
corner = [x y];
end
%----------------
function [isolated_region, filled_region] = getRegionImages(region, im_threshed);
%create an image with just the filled in card.
ox = round(region.BoundingBox(1));
oy = round(region.BoundingBox(2));
ex = round(region.BoundingBox(3) + ox);
ey = round(region.BoundingBox(4) + oy);
filled_region = logical(zeros(size(im_threshed)));
filled_region(oy:ey-1, ox:ex-1) = region.FilledImage;
isolated_region = filled_region & logical(im_threshed);
end
%--------
function [ upright_card ] = UprightCard( original, orig_points )
% UprightCard takes an original image and the four corners of the original
% image and calculates the projective transform that maps each point to an
% upright card. Then this transform is applied, and the result is returned

% We need four corners to do a projective transform
if size( orig_points ) ~= [4 2]
    disp 'Not enough corners!'
end

%Desired aspect ratios
card_height = 700;
card_width = 500;

desired_points = [ 0 0; card_width 0; card_width card_height; 0 card_height];
T = maketform('projective', orig_points, desired_points );

upright_card = imtransform(original, T, 'XData', [1 card_width],'YData', [1 card_height]);

end
%---------------
function [ corner ] = CropToCorner( card )
% CropToCorner takes a 500x700 card and returns the top-left corner
% side-by-side with a 180 degree rotated bottom-left corner.
%

corner_height = 180;
corner_width = 70;

%Get the top-left corner
corner1 = imcrop(card, [0 0 70 180 ]);
%Get the rotated bottom-right corner
corner2 = imcrop(card, [size(card,2)-corner_width+1  size(card,1)-corner_height+1 size(card,2) size(card,1) ]);
corner2 = imrotate(corner2,180);
%Place them next to each other
corner = [corner1 corner2];

end
%----------------------
function [ bestMatch, maxCorr, confidence, indexX, indexY ] = RunTemplateMatching( imageDft, templatesDft, templatesDftConj )
% RunTemplateMatching attempts to match templates against an image. It
% accepts a vector of templates, and returns the best match index, along
% with some ancillary information.
%
%   imageDft should be the DFT of the test image.
%   templatesDft should be a vector of DFTs of templates
%   templatesDftConj should be the conjugates of the templatesDft.
%   Both templatesDft and templatesDftConj should be the same size as the
%   imageDft.
%   The returned maxCorr value is the maximum correlation, indexX and
%   indexY are the location of the maximum correlation.
%   confidence is the maximum correlation detected divided by the second
%   best correlation detected.

% initialize values to 0.
maxCorr = 0;
runnerUpCorr = 0;
indexX = 0;
indexY = 0;
bestMatch = 0;

% Check the correlation for each template.
for k = 1:length(templatesDft)
    [corr deltaX deltaY ] = PhaseCorrelation(imageDft, templatesDft{k}, templatesDftConj{k} );
    
    % New maximum value?
    if corr > maxCorr
        runnerUpCorr = maxCorr;
        maxCorr = corr;
        indexX = deltaX;
        indexY = deltaY;
        bestMatch = k;
    elseif corr > runnerUpCorr
        runnerUpCorr = corr;
    end
end

confidence = maxCorr / runnerUpCorr;

end
%----------------------
function [ correlation, delta_x, delta_y ] = PhaseCorrelation( im1Dft, im2Dft, im2DftConj )
% PhaseCorrelation Calculates the phase correlation between an image and
% a template, returning the maximum correlation and the offset x and y
% locations.
%
%   Images are passed in as DFTs, first the DFT of the image, then the DFT
%   of the template and the complex conjugate of the template DFT. All DFTs
%   should be the same size.

% Make sure we're using two dimensions and that the sizes are the same.
if ndims(im1Dft) ~= 2
    disp('im1 has wrong dimenstions!');
    return;
end
if sum(size(im1Dft) == size(im2DftConj)) ~= 2
    disp('im1 and im2 are not the same size!');
    return;
end

% Calculate the denominator.
H = abs(im1Dft) .* abs(im2Dft);

% Calculate the overall result.
result = ifft2(im1Dft .* im2DftConj ./ H);

% Find the max and its location.
[correlation, maxIdx] = max(result(:));
[delta_y, delta_x] = ind2sub(size(result), maxIdx);

end
%----------------------

%----------------------
