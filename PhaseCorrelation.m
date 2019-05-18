function [ correlation, delta_x, delta_y ] = PhaseCorrelation(im1Dft, im2Dft, im2DftConj)
% PhaseCorrelation Calculates the phase correlation between an image and
% a template, returning the maximum correlation and the offset x and y
% locations.
%
%   Images are passed in as DFTs, first the DFT of the image, then the DFT
%   of the template and the complex conjugate of the template DFT. All DFTs
%   should be the same size.

% Make sure we're using two dimensions and that the sizes are the same.
if ndims(im1Dft) ~= 2
    disp('Error@PhaseCorrelation: im1 has wrong dimenstions!');
    return;
end
if sum(size(im1Dft) == size(im2DftConj)) ~= 2
    disp('Error@PhaseCorrelation: im1 and im2 are not the same size!');
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

