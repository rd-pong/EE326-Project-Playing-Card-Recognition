function [correlation, delta_x, delta_y] = PhaseCorrelation(sourceDft, templateDft, templateDftConj)
% Input: the discrete fourier transform of source imamge, template image
% and conjugate of template image

% Output: peak correlation, corrdinator of peak correlation

% Target: calculates the phase correlation between an image and a template,
% returning the maximum correlation and the offset x and y locations.

% Judging the size of 2 input
if ndims(sourceDft) ~= 2
    disp('Error@PhaseCorrelation: im1 has wrong dimenstions!');
    return;
end

if sum(size(sourceDft) == size(templateDftConj)) ~= 2
    disp('Error@PhaseCorrelation: im1 and im2 are not the same size!');
    return;
end

% Calculate the denominator.
H = abs(sourceDft) .* abs(templateDft);

% compute overall coorelation.
result = ifft2(sourceDft .* templateDftConj ./ H);

% Find peak value and its location.
[correlation, maxIdx] = max(result(:));
[delta_y, delta_x] = ind2sub(size(result), maxIdx);

end

