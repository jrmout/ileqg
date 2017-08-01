function [x,y] = computeEllipse(xc, yc, e1, e2, N)
% ELLIPSE - draw and ellipse on specified or current axis
%
% H = ELLIPSE(...) returns a handle to the plotted ellipse / arc.

%Number of points, arc start*end angles (in radians)
th = [0 2*pi];  %default: full ellipse

% x and y are the center of the eclipse
% e1 and e2 are the eigenvectors

% Compute the rotation t of the ellipse
if norm(e1) > norm(e2)
    t = atan2(e1(2),e1(1));
    rx = norm(e1);
    ry = norm(e2);
else 
    t = atan2(e2(2),e2(1));
    rx = norm(e2);
    ry = norm(e1);
end



%distribute N points between arc start & end angles
th = linspace(th(1),th(2),N);

%calculate x and y points
x = xc + rx*cos(th)*cos(t) - ry*sin(th)*sin(t);
y = yc + rx*cos(th)*sin(t) + ry*sin(th)*cos(t);