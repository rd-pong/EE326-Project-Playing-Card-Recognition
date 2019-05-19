function [corner] = getIntersection(line1, line2)
% Input: 2 line(2 points on the line)

% Output: an array containing corner coordinators

% Target: find intersection mathematically using 2 line equation(k, b)

x = 1;
y = 1;

% error case
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