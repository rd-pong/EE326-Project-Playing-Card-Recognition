%% Preliminaries
clear;
clc;
NUMBER_OF_CARDS_IN_THE_FIGURE = 8;
I = im2double(rgb2gray(imread('training/8cards_top.jpg')));
% test.jpg/8hearts.jpg(3edges)/10hearts.jpg(con-clockwise)/Ahearts.jpg(clockwise)
detectedCards = [];
rankTemplates = {...
    imread('templates/template_11_rank_edge.jpg') 1; ...
    imread('templates/template_47_rank_edge.jpg') 2; ...
    imread('templates/template_17_rank_edge.jpg') 3; ...
    imread('templates/template_7_rank_edge.jpg') 4; ...
    imread('templates/template_14_rank_edge.jpg') 5; ...
    imread('templates/template_6_rank_edge.jpg') 6; ...
    imread('templates/template_15_rank_edge.jpg') 7; ...
    imread('templates/template_9_rank_edge.jpg') 8; ...
    imread('templates/template_4_rank_edge.jpg') 9; ...
    imread('templates/template_10_rank_edge.jpg') 10; ...
    imread('templates/template_1_rank_edge.jpg') 11; ...
    imread('templates/template_12_rank_edge.jpg') 12; ...
    imread('templates/template_39_rank_edge.jpg') 13; ...
    };

suitTemplates = { ...
    imread('templates/template_4_suit_edge.jpg') 1; ...
    imread('templates/template_6_suit_edge.jpg') 3; ...
    imread('templates/template_8_suit_edge.jpg') 4; ...
    imread('templates/template_13_suit_edge.jpg') 2; ...
    };

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

%% Iteration
for nRegion = 1:1
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
    corners = FindCorners(I, imcropppedToCard, 0);
    
    %% Projective Transform
    %----------garanteen corners 1-4 rotate clockwise-------------------
    % Calculate cartisian to polar(theta)
    a = corners-statsCropped(1).Centroid; % coordinates of x to the right, y to the down
    a = [a(:, 1), -a(:, 2)]; % coordinates of x to the right, y to the up
    theta = cart2pol(a(:, 1), a(:, 2));
    sortBuffer = sortrows([theta, corners(:, 1), corners(:, 2)], 'descend');
    corners = [sortBuffer(:, 2), sortBuffer(:, 3)];
    
    % plot corners arraged to be clockwise
    figure, imshow(imcropppedToCard), title('Arranged Corners'), hold on
    [a,b] = size(corners);
    for zz = 1:a
        plot(corners(zz,1), corners(zz,2), 'x','LineWidth',2,'Color','yellow')
        text(corners(zz,1), corners(zz,2), num2str(zz), 'FontSize', 20, 'Color', 'red');
    end
    %-----------------------------------
    
    %----------make sure line 4(1-2) is the smallest edge
    while norm(corners(1, :)-corners(4, :)) < norm(corners(1, :)-corners(2, :))
        % if the corners are rotated counterclockwise
        newCorners = circshift(corners, 1);
        corners = newCorners;
    end
    
    % plot corners 2 line 1(1-4) is the smallest edge
    figure, imshow(imcropppedToCard), title('Arranged Corners'), hold on
    [a,b] = size(corners);
    for zz = 1:a
        plot(corners(zz,1), corners(zz,2), 'x','LineWidth',2,'Color','yellow')
        text(corners(zz,1), corners(zz,2), num2str(zz), 'FontSize', 20, 'Color', 'red');
    end
    %----------
    
    
    
%     % garanteen corners 1-4 rotate clockwise
%     if norm(corners(1, :)-corners(4, :)) < norm(corners(1, :)-corners(2, :))
%         % if the corners are rotated counterclockwise
%         newCorners = [corners(1, :); corners(4, :); corners(3, :); corners(2, :)];
%         corners = newCorners;
%     end
    % Perform projection
    card = UprightCard(imcropppedToCardOri, corners);
    figure, imshow(card);
    
end
