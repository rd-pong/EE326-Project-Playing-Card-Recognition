function z = plotEdges(Image, lines)
% Input: image for margin finding, lines array containing start and end
% poing of found lines

% Output: a figure in gui

% Target: plot the margins and its order

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