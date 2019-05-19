function [k, b] = getLineEquation(x1, x2, y1, y2)
% Input: 2 coordinators

% Output: the slope and intersection of y axis

% Target: get 2 parameters for a line equation given 2 points

if abs(x1 - x2) < 1e-6
    k = -1;
    b = -1;
else
    kb = [x1 1; x2 1]\[y1;y2];
    k = kb(1);
    b = kb(2);
end

end