% crop template
% imwrite(imcrop(cornerEdges,[0,0,67,105]), 'templates/2.jpg')
% imwrite(imcrop(cornerEdges,[0,105,67,73]), 'templates/club.jpg')

%% File read in, set parameters
debug = 0;
NUMBER_OF_CARDS_IN_THE_FIGURE = 1;
I = im2double(rgb2gray(imread('training/2clubs.jpg')));

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
%% Preliminaries
BW = imbinarize(I, 0.66);
% Remove small components of the binary image
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
    BW2=ismember(imLabeled,index(1:NUMBER_OF_CARDS_IN_THE_FIGURE));
end
imLabeled2 = bwlabel(BW2);
stats2 = regionprops(imLabeled2, 'Centroid', 'Area', 'Image', 'FilledImage', 'BoundingBox');
% Reverse complement
% BW2 = 1-BW2;
% figure;
% imshowpair(BW, BW2, 'montage');
% title('Binary Image & Edited Binary Image');

%% Find rankTemplatesConj & suitTemplatesConj
templateCornerSize = [105, 67];
dilateSize = 2;
for i = 1:length(rankTemplates)
    % Dilate the template
    rankTemplates{i,1} = imdilate(rankTemplates{i,1}, strel('disk', dilateSize));
    
    % Calculate the FFT. Phase correlation requires removing the
    % mean from the template.
    meanTemplate = mean2(rankTemplates{i,1});
    rankTemplates{i,1} = fft2(rankTemplates{i,1} - meanTemplate, templateCornerSize(1), templateCornerSize(2));
    
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
        meanTemplate, templateCornerSize(1), templateCornerSize(2));
    currentSuitTemplatesConj{i,1} = conj(suitTemplates{i,1});
end
rankTemplatesConj = currentRankTemplatesConj;
suitTemplatesConj = currentSuitTemplatesConj;
%% Iteration
for nRegion = 1:NUMBER_OF_CARDS_IN_THE_FIGURE
    imcropppedToCardOri = imcrop(imbinarize(I), stats2(nRegion).BoundingBox);
    % ignore small region
    imLabeledCropped = bwlabel(imcropppedToCardOri);
    statsCropped = regionprops(imLabeledCropped,'Centroid', 'Area', 'Image', 'FilledImage', 'BoundingBox');
    [b,index]=sort([statsCropped.Area],'descend');
    if length(statsCropped)<1
        imcropppedToCard=imLabeledCropped;
    else
        imcropppedToCard=ismember(imLabeledCropped,index(1:1));
    end
    corners = FindCorners(I, imcropppedToCard, debug);
    
    %% Projective Transform
    %----------garanteen corners 1-4 rotate clockwise-------------------
    % Calculate cartisian to polar(theta)
    a = corners-statsCropped(1).Centroid; % coordinates of x to the right, y to the down
    a = [a(:, 1), -a(:, 2)]; % coordinates of x to the right, y to the up
    theta = cart2pol(a(:, 1), a(:, 2));
    sortBuffer = sortrows([theta, corners(:, 1), corners(:, 2)], 'descend');
    corners = [sortBuffer(:, 2), sortBuffer(:, 3)];
    
    %     % plot corners arraged to be clockwise
    %     figure, imshow(imcropppedToCard), title('Arranged Corners'), hold on
    %     [a,b] = size(corners);
    %     for zz = 1:a
    %         plot(corners(zz,1), corners(zz,2), 'x','LineWidth',2,'Color','yellow')
    %         text(corners(zz,1), corners(zz,2), num2str(zz), 'FontSize', 20, 'Color', 'red');
    %     end
    %-----------------------------------
    
    %----------make sure line 4(1-2) is the smallest edge
    while norm(corners(1, :)-corners(4, :)) < norm(corners(1, :)-corners(2, :))
        % if the corners are rotated counterclockwise
        newCorners = circshift(corners, 1);
        corners = newCorners;
    end
    
    %     % plot corners 2 line 1(1-4) is the smallest edge
    %     figure, imshow(imcropppedToCard), title('Arranged Corners'), hold on
    %     [a,b] = size(corners);
    %     for zz = 1:a
    %         plot(corners(zz,1), corners(zz,2), 'x','LineWidth',2,'Color','yellow')
    %         text(corners(zz,1), corners(zz,2), num2str(zz), 'FontSize', 20, 'Color', 'red');
    %     end
    %----------
    
    % Perform projection
    card = UprightCard(imcropppedToCardOri, corners);
    
    %figure, imshow(card);
    
    %% Crop to Corners(4 corners)
    % Crop the image to 2 corners
    croppedCorners = CropToCorner(card);
    
    figure, imshow(croppedCorners);
    
    % Detect edges
    cornerEdges = edge(croppedCorners, 'prewitt');
    
    %% Rcognition
    % 注意：这里是新的cornerEdges
    cornerEdges = edge(croppedCorners, 'canny');
    
    % FFT with the right size.
    cornerDft = fft2(cornerEdges, size(rankTemplates{1,1},1), size(rankTemplates{1,1},2));
    % if the corner DFT is 0, then our phase correlator will have problems (divide by 0). Just set 0s to very small values.
    cornerDft(cornerDft == 0) = 0.00001;
    bestMatchRank = RunTemplateMatching(cornerDft, rankTemplates(:,1), rankTemplatesConj(:,1));
    cornerDft = fft2(cornerEdges, size(rankTemplates{1,1},1), size(rankTemplates{1,1},2));
    bestMatchSuit = RunTemplateMatching(cornerDft, suitTemplates(:,1), suitTemplatesConj(:,1));
    
    fprintf('Identified rank: %d suit: %d, ', rankTemplates{bestMatchRank, 2}, suitTemplates{bestMatchSuit, 2});
    
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
    newRow = [rankTemplates{bestMatchRank, 2} suitTemplates{bestMatchSuit, 2}];
    detectedCards = [detectedCards; newRow];
end
