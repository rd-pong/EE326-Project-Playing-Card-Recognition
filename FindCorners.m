function corners = FindCorners(I, BW2, debug)
% Input: original image to be shown when debug, bianry image 
% Output: corner coodinators in an double array
% Target: Find margins & their equations -> find intersection -> find corners
if nargin <= 2
    debug = 0;
end
%% Find lines
findingLine = edge(BW2, 'canny');
[H,T,R] = hough(findingLine);
P  = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));
% change parameter if different data set
lines = houghlines(findingLine,T,R,P,'FillGap',3000,'MinLength',200);
% lines = houghlines(findingLine,T,R,P,'FillGap',3000);

if debug
    plotEdges(BW2, lines);
end

%% Arrange margins to be adjacent in the order of 1234
% line1-line2-line3-line4 are adjacent margins
parallel_threshold = 10;
line1 = lines(1);
line2 = [];
line3 = [];
line4 = [];
if abs(lines(2).theta - line1.theta) < parallel_threshold
    line2 = lines(3);
    line3 = lines(2);
    line4 = lines(4);
elseif abs(lines(3).theta - line1.theta) < parallel_threshold
    line2 = lines(2);
    line3 = lines(3);
    line4 = lines(4);
else
    line2 = lines(2);
    line3 = lines(4);
    line4 = lines(3);
end
lines = [line1; line2; line3; line4];
if debug
    plotEdges(BW2, lines);
end
%% Get line equation -> find intersection -> find corners
lineEquations = [];
for j = 1:length(lines)
    line = lines(j);
    x1 = line.point1(1);
    y1 = line.point1(2);
    x2 = line.point2(1);
    y2 = line.point2(2);
    [k, b] = getLineEquation(x1, x2, y1, y2);
    lineEquations = [lineEquations, struct('k', k, 'b', b, 'x', x1, 'y', y1)];
end

corners = [];
corners(1,:) = getIntersection(lineEquations(1), lineEquations(2));
corners(2,:) = getIntersection(lineEquations(2), lineEquations(3));
corners(3,:) = getIntersection(lineEquations(3), lineEquations(4));
corners(4,:) = getIntersection(lineEquations(1), lineEquations(4));

if debug
    % plot corners
    figure, imshow(BW2), title('UnArranged Corners'), hold on
    [a,b] = size(corners);
    for zz = 1:a
        plot(corners(zz,1), corners(zz,2), 'x','LineWidth',2,'Color','yellow')
        text(corners(zz,1), corners(zz,2), num2str(zz), 'FontSize', 20, 'Color', 'red');
    end
end

end
%----------------
function [k, b] = getLineEquation(x1, x2, y1, y2)
if abs(x1 - x2) < 1e-6
    k = -1;
    b = -1;
else
    kb = [x1 1; x2 1]\[y1;y2];
    k = kb(1);
    b = kb(2);
end
end
%%-------------------
function z = plotEdges(Image, lines)
figure, imshow(Image), hold on
max_len = 0;
for k = 1:length(lines)
    xy = [lines(k).point1; lines(k).point2];
    plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
    % Plot beginnings and ends of lines
    text(xy(1,1),xy(1,2),num2str(k), 'FontSize', 20, 'Color', 'red');
    plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
    plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
    % Determine the endpoints of the longest line segment
    len = norm(lines(k).point1 - lines(k).point2);
    if ( len > max_len)
        max_len = len;
        xy_long = xy;
    end
end
end
%----------------
function [corner] = getIntersection(line1, line2)
x = 1;
y = 1;
% error
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
%----------------