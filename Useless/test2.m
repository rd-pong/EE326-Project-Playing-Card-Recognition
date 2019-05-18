clear;
clc;
%% Some constants used in our algorithms
inputImagePath = 'training/9clubs.jpg';
% test.jpg/8hearts.jpg(3edges)/10hearts.jpg(con-clockwise)/Ahearts.jpg(clockwise)
outputImagePath = 'testOut/';

CONTRAST_GAIN = 2;
BLACK_REGION_SIZE = 1000000;
WHITE_REGION_SIZE = 10000;
EDGE_PROXIMITY = .1;
MAX_DETECTED_CARDS = 5;
ANGLE_THRESHOLD = .3;
INITIAL_EDGE_ALG = 'prewitt';
SECOND_EDGE_ALG = 'canny';

%%
% Read the image. Increase constrast.
inputImage = im2double(rgb2gray(imread(inputImagePath)));
% Threshold using Otsu.
% inputBw = im2bw(inputImage, graythresh(inputImage));
inputBw = imbinarize(inputImage, 0.66);
% Remove small black and white regions.
% inputBw = ~bwareaopen(~inputBw, BLACK_REGION_SIZE);
% inputBw = bwareaopen(inputBw, WHITE_REGION_SIZE);

%%
% Remove regions close to the edges
imageLabeled = bwlabel(inputBw);
[~, badRegions] = FindNonEdgeRegions(imageLabeled, EDGE_PROXIMITY);
for j = 1:length(badRegions)
    inputBw(imageLabeled == badRegions(j)) = 0;
end

% Re-label, find the centroids
[imageLabeled numRegions] = bwlabel(inputBw);
props = regionprops(imageLabeled, 'Centroid', 'BoundingBox');
numValidRegions = 0;

detectedCards = [];
% Iterate over each currently reasonable region
for nRegion = 1:numRegions
    % If we've found the maximum number of cards, quit.
    if numValidRegions == MAX_DETECTED_CARDS
        fprintf('Found %d cards!\n', MAX_DETECTED_CARDS);
        break;
    end
    
    % Create a mask of solely this region
    mask = zeros(size(inputBw));
    mask(imageLabeled == nRegion) = 1;
    
    % Crop the image down so we're working with a smaller
    % image
    mask = imcrop(mask, props(nRegion).BoundingBox);
    imageCropped = imcrop(inputImage, props(nRegion).BoundingBox );
    
    % Shift the calculated centroids, since we cropped the image
    props(nRegion).Centroid(1) =  props(nRegion).Centroid(1) - props(nRegion).BoundingBox(1);
    props(nRegion).Centroid(2) =  props(nRegion).Centroid(2) - props(nRegion).BoundingBox(2);
    
    % Find the corners of the region.
    % Corners are returned in increasing order of angle between the
    % centroid and the corner--so they will be in clockwise order.
    % Returns at most 4 corners.
    [angles, ~, corners] = FindCorners(mask, props(nRegion).Centroid);
    
    % Need at least 4 corners to transform it
    if size(corners,1) < 4
        fprintf('Fewer than 4 corners\n');
        continue;
    end
    
    % angles returns sort increasing order
    % so diff should -pi'ish.
    if (abs(angles(1) - angles(3) + pi()) > ANGLE_THRESHOLD) || ...
            (abs(angles(2) - angles(4) + pi()) > ANGLE_THRESHOLD)
        fprintf('Shape is not card-like\n');
        continue;
    end
    
    % Count the edges in the "upper left" corner. Compare to the "upper
    % right" corner. Call the corner with more edges they actual "upper
    % left" corner.
    
    % Perform the projective transform.
    card = UprightCard(imageCropped, corners);
    % Crop the image so we only have the corner
    corner = CropToCorner(card);
    % Detect edges
    cornerEdges = edge(rgb2gray(corner), INITIAL_EDGE_ALG);
    % Count edge pixels.
    sum1 = sum(cornerEdges(:));
    
    
    % Now, do the same as above for the upper right corner.
    % Rotate the corners by 1 so we get a different corner.
    cornersRotated = [corners(2:end,:); corners(1,:)];
    cardRotated = UprightCard(imageCropped, cornersRotated);
    cornerUpperRight = CropToCorner(cardRotated);
    cornerEdgesUpperRight = edge(rgb2gray(cornerUpperRight), INITIAL_EDGE_ALG);
    sum2 = sum(cornerEdgesUpperRight(:));
    
    % Whichever has more edge pixels is probably the upper left corner.
    if sum2 > sum1
        fprintf('Using upper right corner\n');
        corner = cornerUpperRight;
    end
    
    % Detect edges in the real upper left corner.
    cornerEdges = edge(rgb2gray(corner), SECOND_EDGE_ALG);
    
    % FFT with the right size.
    cornerDft = fft2(cornerEdges, size(rankTemplates{1,1},1), size(rankTemplates{1,1},2));
    % if the corner DFT is 0, then our phase correlator will
    % have problems (divide by 0). Just set 0s to very small values.
    cornerDft(cornerDft == 0) = 0.00001;
    bestMatchRank = RunTemplateMatching(cornerDft, rankTemplates(:,1), rankTemplatesConj(:,1));
    bestMatchSuit = RunTemplateMatching(cornerDft, suitTemplates(:,1), suitTemplatesConj(:,1));
    
    fprintf('Identified %d of %d\n', rankTemplates{bestMatchRank, 2}, suitTemplates{bestMatchSuit, 2});
    
    % If we detect the same card twice...
    % Well. Something is wrong. Let's throw the result out.
    newRow = [rankTemplates{bestMatchRank, 2} ...
        suitTemplates{bestMatchSuit, 2}];
    
    if sum(ismember(detectedCards, newRow, 'rows')) > 0
        fprintf('Found duplicate\n');
    else
        detectedCards = [detectedCards; newRow];
        % Increment number of cards detected
        numValidRegions = numValidRegions + 1;
    end
end

% Write the output file
% Write the detected cards
fid = fopen(outputImagePath, 'w');
for j = 1:numValidRegions
    fprintf(fid, '%d %d|', detectedCards(j, 1), detectedCards(j, 2));
end
if numValidRegions == 0
    disp('None detected!');
    fprintf(fid, 'None detected!');
end
fprintf(fid, '\n');

% If we have enough regions, calculate statistics
if (numValidRegions > 2)
    cards = zeros(1, numValidRegions);
    % Convert the cards into the representation required by
    % GetMostLikelyHands.
    % 0-12 = A, 2...K of clubs
    % 13-25 = A, 2 ... K of diamonds
    % and so on.
    for j = 1:numValidRegions
        cards(j) = (detectedCards(j,1) - 1) + ((detectedCards(j,2) - 1) * 13);
    end
    
    % Get the odds and print them.
    [hands, odds] = GetMostLikelyHands(cards);
    for j = 1:length(odds);
        fprintf(fid, '%s: %g%%\n', handNames{hands(j) + 1}, odds(j)*100);
    end;
end

fclose(fid);
%% Fuctions
%-----------------------------------------------------------------------------
function [ goodRegions, badRegions ] = FindNonEdgeRegions( labeledImage, percentFromEdge )
% FindNonEdgeRegions Takes a labelled image and identifies which regions
% are located near the edge, and which aren't.
%
%   percentFromEdge delimits where "close to the edge" is. If the centroid
%   of the region is located within percentFromEdge of the edge, it is a
%   "bad region", otherwise it is good.

goodRegions = [];
badRegions = [];

% Determine good/bad boundaries.
xMin = percentFromEdge * size(labeledImage, 2);
xMax = (1 - percentFromEdge) * size(labeledImage, 2);

yMin = percentFromEdge * size(labeledImage, 1);
yMax = (1 - percentFromEdge) * size(labeledImage, 1);

% Iterate over each region's centroid, exercising great justice.
props = regionprops(labeledImage, 'Centroid');
for nRegion = 1:length(props)
    cent = props(nRegion).Centroid;
    if (cent(1) > xMin) && (cent(1) < xMax) && ...
            (cent(2) > yMin) && (cent(2) < yMax)
        goodRegions = [goodRegions nRegion];
    else
        badRegions = [badRegions nRegion];
    end
end

end
%-----------------------------------------------------------------------------
function [ angles distances corners ] = FindCorners( cardMask, centroid )
% FindCorners Takes a mask of a card's location, along with the centroid of
% the card and outputs the angles to, distances to, and corner locations.
%   
%   computeCards takes the card mask and dilates it very slightly. This is
%   then XOR'd with the orignal mask to get the outline of the mask. From
%   the outline, we can find the angles to outline pixels from the
%   centroid, as well as the distance. By sorting the distance values using
%   the angle, we can look for local maxima. These will be corners. Limit
%   ourselves to four corners (cards only have four), and return the
%   corners in increasing angle from the centroid.
     
angles = [];
distances = [];
corners = [];

LOCAL_MAXIMA_RANGE = 40;

% find a very narrow edge to the card
maskEroded = imdilate(cardMask,strel('disk',1));
edges = (cardMask ~= maskEroded);
        
% find the x and y coordinates of the edges
[edgesY, edgesX] = ind2sub(size(cardMask), find(edges == 1));

% calculate the distance from the centroid for the edges
% and the angle from the centroid. The angle is useful
% for sorting the distances as one rotates around the edge
dist = zeros(1,length(edgesX));
ang = zeros(1,length(edgesX));
for q = 1:length(edgesX)
    dist(q) = sum(([edgesX(q) edgesY(q)] - centroid).^2);
    ang(q) = atan2(edgesY(q) - centroid(2),edgesX(q) - centroid(1));
end   
        
% combine the angle, distances, and edge points.
% sort by angle.
combined = sortrows([ang; dist; edgesX'; edgesY']', 1);

% This can occur when the whole image is passed to FindCorners
if size(combined,1) < 1
    return;
end;

% find the local maxima
peaks = FindPeakIndices(combined(:,2), LOCAL_MAXIMA_RANGE);        

% If there are more than four, remove the smallest ones.
combined = combined(peaks, :);
combined = sortrows(combined, 2);
if size(combined, 1) > 4
    combined = combined(size(combined,1) - 3:end, :);
end
combined = sortrows(combined, 1);

% Split the values out.
angles = combined(:, 1);
distances = combined(:, 2);
corners = combined(:, 3:4);

end
%-----------------------------------------------------------------------------
function [ indices ] = FindPeakIndices( inVec, delta )
% FindPeakIndices Identifies local maxima within a sorted vector.
%   
%   inVec contains the data from which to find the maxima. inVec needs to
%   be sorted in the appropriate (to the specific data) way, prior to
%   calling FindPeakIndices. delta specifies the number of samples to each
%   side that a pixel must be larger than.

% Create a new vector containing the original vector, plus extra 2*delta
% from the beginning of inVec concatenated to the end.
newVec = [inVec; inVec(1:2*delta,:)];

indices = [];
% For each pixel location (in the new vector) containing delta pixels to
% each side, check if this pixel is the maximum amonst the surrounding
% pixels.
for i = delta+1:length(inVec)+delta
    
    % Get a smaller vector of the data
    subVec = newVec(i-delta:i+delta);
    
    % If the maximum pixel in the sub vector is this one, it is a local
    % maximum.
    if max(subVec) == subVec(delta+1)
        val = i;
        % If the index is past the end of the original vector, then it was
        % a pixel in the first delta pixels.
        if i > length(inVec)
            val = i - length(inVec);
        end
        indices = [indices val];
    end
end

end
%-----------------------------------------------------------------------------
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

%-----------------------------------------------------------------------------

%-----------------------------------------------------------------------------