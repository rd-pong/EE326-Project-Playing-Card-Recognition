function [ bestMatch, maxCorr, confidence, indexX, indexY ] = RunTemplateMatching( imageDft, templatesDft, templatesDftConj )
% RunTemplateMatching attempts to match templates against an image. It
% accepts a vector of templates, and returns the best match index, along
% with some ancillary information.
%
%   IMAGEDFT should be the DFT of the test image. TEMPLATEDFT should be a
%   vector of DFTs of templates TEMPLATEDFTCONJ should be the conjugates
%   of the templatesDft. Both templatesDft and templatesDftConj should be
%   the same size as the imageDft. The returned MAXCORR value is the
%   maximum correlation, INDEXX and INDEXY are the location of the maximum
%   correlation. 
%   CONFIDENCE is the maximum correlation detected divided by
%   the second best correlation detected.

% initialize values to 0.
maxCorr = 0;
runnerUpCorr = 0;
indexX = 0;
indexY = 0;
bestMatch = 0;  

% Check the correlation for each template.
for k = 1:length(templatesDft)
    [corr deltaX deltaY ] = PhaseCorrelation( imageDft, templatesDft{k}, templatesDftConj{k} );
    
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

