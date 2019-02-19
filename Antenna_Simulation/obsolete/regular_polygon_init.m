function [vert, poly ] = regular_polygon_init(N,edge_length)
%initialize poloygon
% 1st polygon

preci = 5; % Decimal digits. Only for the coordinates of the vertices.

poly.cen =[0,0,0]; % define Polygon Center 
poly.N = N;  % define the shape - number of vertices (6-gon, 8-gon, N-gon)
poly.ang_centr = 2*pi / poly.N; % Central Angle of the polygon (radians)
poly.ang_inter = pi-(2*pi/poly.N); % Interior Angle of the polygon
poly.side_length = edge_length;  % define the side length in mm (actually in the chosen unit - see "unit" variable)
poly.r = poly.side_length/(2*sin(pi/poly.N)); % radius of the polygon (see inscribed polygons)
poly.rotation = 0; %rad
cent_from_origin = pi/4;
theta = 0;
poly.cen =cent_from_origin * [ cos(theta),  sin(theta), 0];

for i=1:poly.N
    temp = poly.rotation + (i-1) * poly.ang_centr;
    poly.vert(:,i) = [poly.cen(1), poly.cen(2)] + poly.r * [cos(temp), sin(temp)];
end

poly.vert = round(poly.vert,preci);  % NECESSARY!!! 
% A very small difference between the coordinates of two vertices (due to
% aproximation), will cause an extream discretization of the grid.
% If more precision is required, the second argument may be increased.
vert = poly.vert;
end