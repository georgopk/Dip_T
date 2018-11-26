%%initialize poloygon

poly.cen =[0,0,0]; % define Polygon Center 
poly.N = 16;  % define the shape - number of vertices (6-gon, 8-gon, N-gon)
poly.ang_centr = 2*pi / poly.N; % Central Angle of the polygon (radians)
poly.ang_inter = pi-(2*pi/poly.N); % Interior Angle of the polygon
poly.side_length = 30;  % define the side length in mm (actually in the chosen unit - see "unit" variable)
poly.r = poly.side_length/(2*sin(pi/poly.N)); % radius of the polygon (see inscribed polygons)
poly.rotation = 0.7; % (radians) if this parameter is 0, the first Vertex of the polygon will be on the Y'Y axis
for i=1:poly.N
    temp = poly.rotation + (i-1) * poly.ang_centr;
    poly.vert(:,i) = [poly.cen(1), poly.cen(2)] + poly.r * [cos(temp), sin(temp)];
end

