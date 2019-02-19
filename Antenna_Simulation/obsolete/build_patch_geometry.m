%% initialize complex Geometry
% create a polygon and two transmission lines


preci = 5; % Decimal digits. Only for the coordinates of the vertices.

% initialize Line 1 
% width -> 5 mm,  length -> 15 mm
line(1).length = 17.5; %SPECIAL LENGTH
line(1).width = 5;

% initialize Line 2
line(2).length = 45;
line(2).width = 10;

%initialize 1st poloygon
poly.N = 16;  % define the shape - number of vertices (6-gon, 8-gon, N-gon)
poly.ang_centr = 2*pi / poly.N; % Central Angle of the polygon (radians)
poly.ang_inter = pi-(2*pi/poly.N); % Interior Angle of the polygon
poly.side_length = 30;  % define the side length in mm (actually in the chosen unit - see "unit" variable)
poly.r = poly.side_length/(2*sin(pi/poly.N)); % radius of the polygon (see inscribed polygons)
poly.rotation = 0; %rad
cent_from_origin = line(1).length + line(2).length + poly.r;
theta = pi/2 - 1.5 * poly.ang_centr;
poly.cen =cent_from_origin * [ cos(theta),  sin(theta), 0];

for i=1:poly.N
    temp = poly.rotation + (i-1) * poly.ang_centr;
    poly.vert(:,i) = [poly.cen(1), poly.cen(2)] + poly.r * [cos(temp), sin(temp)];
end

poly.vert = round(poly.vert,preci);  % NECESSARY!!! 
% A very small difference between the coordinates of two vertices (due to
% aproximation), will cause an extream discretization of the grid.
% If more precision is required, the second argument may be increased.


% line 1 placement
% !SPECIAL!
rot = pi/2 - theta;
% start = [0,0]; % special
stop = line(1).length * [cos(pi/2 - rot),sin(pi/2 - rot)];

line(1).vert(:,1) = line(1).width * sqrt(2)/2 * [cos(-pi/4 -rot), sin(-pi/4 - rot)];  % (special)
line(1).vert(:,2) = stop + line(1).width/2 * [cos(-rot), sin(-rot)];
line(1).vert(:,3) = stop + line(1).width/2 * [cos(pi-rot), sin(pi-rot)];
line(1).vert(:,4) = line(1).width * sqrt(2)/2 * [cos(-3*pi/4 -rot), sin(-3*pi/4 - rot) ]  ;% (special)


% line 2 placement
start = stop;  % attached to the first line
stop = start + line(2).length * [cos(pi/2 -rot), sin(pi/2 - rot)]+2;

line(2).vert(:,1) = start + line(2).width/2 * [cos(-rot), sin(-rot)];
line(2).vert(:,2) = stop + line(2).width/2 * [cos(-rot), sin(-rot)];
line(2).vert(:,3) = stop + line(2).width/2 * [cos(pi-rot), sin(pi-rot)];
line(2).vert(:,4) = start + line(2).width/2 * [cos(pi-rot), sin(pi-rot)];

