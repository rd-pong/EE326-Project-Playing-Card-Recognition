function [upright_card] = UprightCard(original, orig_points)
% Input: original image, coordinators of the corners

% Output: 700 * 500 image

% Target: maps each point to an 700* 500 image

% We need four corners to do a projective transform
if size( orig_points ) ~= [4 2]
    disp 'Not enough corners!'
end

%Desired aspect ratios
card_height = 700;
card_width = 500;

desired_points = [ 0 0; card_width 0; card_width card_height; 0 card_height];
%T = maketform('projective', orig_points, desired_points );
T = fitgeotrans(orig_points, desired_points,'projective');

% upright_card = imtransform(original, T, 'XData', [1 card_width],'YData', [1 card_height]);
RA = imref2d([card_height, card_width], [1 card_width], [1 card_height]);
upright_card = imwarp(original, T, 'OutputView', RA);

end

