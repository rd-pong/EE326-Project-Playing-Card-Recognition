function [ corner ] = CropToCorner( card )
% Input: 500x700 card
% Output: top-left corner & 180 degree rotated bottom-left corner
% Target: crop the uprighted card to 2 corners

corner_height = 180;
corner_width = 70;

%Get the top-left corner
% original: corner1 = imcrop(card, [0 0 70 180 ]);
% Move the crop mask to left down side
a = [5 10 70 180];
corner1 = imcrop(card, a);
%Get the rotated bottom-right corner
corner2 = imcrop(card, [size(card,2)-corner_width-a(1)+1  size(card,1)-a(2)-corner_height+1 70 180]);
corner2 = imrotate(corner2,180);
%Place them next to each other
corner = [corner1 corner2];

end

