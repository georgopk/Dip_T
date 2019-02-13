%% Antenna coupling (based on...) % Simple Patch Antenna Tutorial
% Simulates Two identical Patch Antennas (back to back).







%% Initialization

close all
clear
clc

%% Setup the Simulation (Initialize Geometry)
physical_constants;
unit = 1e-3; % all length in mm. (unit only for the geometry)

%distance between Antennas (in units)
backDist = 100; % here in mm

%rotation (rotation of the geometry on XY plane. Useful to "fit" a geometry
%on the grid lines)
% rot = 1.5*2*pi/16;
rot = 0;

%piorities
patchPri = 50;
substratePri = 10;
groundPri = 40;
feedPri = 30;

%ground setup
grnd_pos = -12.3; % ground distance from substrate
grnd_points =  [-185.4,-185.4; -185.4,185.4; 185.4,185.4; 185.4,-185.4]' ;

%substrate setup
sub_freq = 8.5e8; % Frequency to calculate the substrate conductivity for
substrate.epsR   = 3.38;
substrate.kappa  = 1e-3 * 2*pi*sub_freq * EPS0*substrate.epsR; %conductivity 
substrate.thickness = 1.524;
substrate.cells = 4;
substrate.points = [-175,-145;175,-145;175,90;-175,90]';

%setup feeding
feed.pos = [18.85, -115.46]; %feeding position in x,y directions
feed.R = 50;     %feed resistance
feed2.pos = feed.pos;
feed2.R = feed.R;

% size of the simulation and dump box 
SimBox = [600 600 150];
dumpWidth = 400;
dumpLength = 400;

%% Setup FDTD Parameter & Excitation Function
f0 = 8.5e8; % center frequency (Hz)
fc = 4e8; % 20 dB corner frequency (Hz) -----> it determines the bandwidth (keep it less than f0)
FDTD = InitFDTD( 'NrTs', 30000, 'EndCriteria', 1e-5);
FDTD = SetGaussExcite( FDTD, f0, fc );
BC = {'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8'}; % boundary conditions
% FDTD = SetBoundaryCond( FDTD, BC,'PML_Grading','-log(1e-6)*log(2.5)/(2*dl*pow(2.5,W/dl)-1) * pow(2.5, D/dl) / Z' );
FDTD = SetBoundaryCond( FDTD, BC );


%% Setup CSXCAD Geometry & Mesh
CSX = InitCSX();

%initialize the mesh with the "air-box" dimensions
mesh.x = [-SimBox(1)/2 SimBox(1)/2];
mesh.y = [-SimBox(2)/2 SimBox(2)/2];
mesh.z = [-SimBox(3)/2 SimBox(3)/2];

% Create Patch
load('points.mat');                     % load a matrix named "points" which contains the coordinates of the geometry
points = points(:,1:length(points)-1);  % remove the redundant last coordinate (same with the first). this line is not necessary. 
points2 = points;
points = rotate_points(points,rot,1);   % rotate the geometry. The third argument enables plot.
hold on;    % Usefull to plot all the geometries in one diagram. (for debugging)
points = sp_round(points,0.2,[feed.pos, 0]); % "round" the coordinates of some vertices, to reduce the mesh
allpoints = points;
pointInd(1) = size(points,2);
CSX = AddMetal( CSX, 'patch' ); % create a perfect electric conductor (PEC)
CSX = AddPolygon(CSX,'patch',patchPri,2,backDist/2+substrate.thickness,points(1:2,:));
% patch 2
points = points2;
points(1,:) = - points2(1,:);          % mirror the patch
points = rotate_points(points,rot,1);   % rotate the geometry. The third argument enables plot.
points = sp_round(points,0.2,[feed.pos, 0]); % "round" the coordinates of some vertices, to reduce the mesh
allpoints = [allpoints, points];
pointInd(2) = pointInd(1)+size(points,2);
CSX = AddMetal( CSX, 'patch2' ); % create a perfect electric conductor (PEC)
CSX = AddPolygon(CSX,'patch2',patchPri+1,2,-(backDist/2+substrate.thickness),points(1:2,:));


% Create Substrate
CSX = AddMaterial( CSX, 'substrate' );
CSX = SetMaterialProperty( CSX, 'substrate', 'Epsilon', substrate.epsR, 'Kappa', substrate.kappa );
points = substrate.points;
points = rotate_points(points,rot,1);
allpoints = [allpoints, points];
pointInd(3) = pointInd(2)+size(points,2);
CSX = AddLinPoly(CSX, 'substrate', substratePri, 2, backDist/2+ 0,points, substrate.thickness);
% add extra cells to discretize the substrate thickness
mesh.z = [linspace(backDist/2+0,backDist/2 + substrate.thickness,substrate.cells+1) mesh.z];
%substrate 2
points = substrate.points;
points(1,:) = - substrate.points(1,:);          % mirror
points = rotate_points(points,rot,1);
allpoints = [allpoints, points];
pointInd(4) = pointInd(3)+size(points,2);
CSX = AddMaterial( CSX, 'substrate2' );
CSX = SetMaterialProperty( CSX, 'substrate2', 'Epsilon', substrate.epsR, 'Kappa', substrate.kappa );
CSX = AddLinPoly(CSX, 'substrate2', substratePri+1, 2, -(backDist/2+ 0),points, -substrate.thickness);
% add extra cells to discretize the substrate thickness
mesh.z = [linspace(-(backDist/2+0),-(backDist/2 + substrate.thickness),substrate.cells+1) mesh.z];

% Create Ground
CSX = AddMetal( CSX, 'gnd' ); % create a perfect electric conductor (PEC)
points = grnd_points;
points = rotate_points(points,rot,1);
allpoints = [allpoints, points];
pointInd(5) = pointInd(4)+size(points,2);
CSX = AddPolygon(CSX, 'gnd', groundPri, 2, backDist/2 + grnd_pos,points);
%ground 2
points = grnd_points;
points(1,:) = -grnd_points(1,:);
points = rotate_points(points,rot,1);
allpoints = [allpoints, points];
pointInd(6) = pointInd(5)+size(points,2);
CSX = AddMetal( CSX, 'gnd2' ); % create a perfect electric conductor (PEC)
CSX = AddPolygon(CSX, 'gnd2', groundPri+1, 2, -(backDist/2 + grnd_pos),points);

% Apply the Excitation & Resist as a Current Source
feed.pos = (rotate_points(feed.pos',rot,1))' ; % adjust the feeding position to the rotated structure
start = [feed.pos, backDist/2 + grnd_pos];
stop = [feed.pos, backDist/2 + substrate.thickness];
[CSX, port{1}] = AddLumpedPort(CSX, feedPri ,1 ,feed.R, start, stop, [0 0 1], true);
% port2
feed2.pos(1) = -feed2.pos(1);
feed2.pos = (rotate_points(feed2.pos',rot,1))' ; % adjust the feeding position to the rotated structure
start = [feed2.pos, -(backDist/2 + grnd_pos)];
stop = [feed2.pos, -(backDist/2 + substrate.thickness)];
[CSX, port{2}] = AddLumpedPort(CSX, feedPri+1 ,2 ,feed2.R, start, stop, [0 0 1]);



% Finalize the Mesh
% -----------------
% detect all edges except of the patch
mesh = DetectEdges(CSX, mesh,'ExcludeProperty','patch');
% detect and set a special 2D metal edge mesh for the patch
mesh = DetectEdges(CSX, mesh,'SetProperty','patch','2D_Metal_Edge_Res', c0/(f0+fc)/unit/50);
% generate a smooth mesh with max. cell size: (e.g) lambda_min / 20 
% mesh = SmoothMesh(mesh, c0/(f0+fc)/unit/20); 
mesh = SmoothMesh(mesh, c0/(f0+fc)/unit/25,'algorithm',[1 3]); % alternative. Useful when SmoothMesh gets stuck in an infinite loop inside "SmoothMeshLines2" function

CSX = DefineRectGrid(CSX, unit, mesh);
CSX = AddDump(CSX,'Hf', 'DumpType', 11, 'Frequency',[9e8]);
CSX = AddBox(CSX,'Hf',10,[-(dumpWidth/2+10) -(dumpLength/2+10) -10*substrate.thickness],[dumpWidth/2+10,  dumpLength/2+10 10*substrate.thickness]); %assign box

% add a nf2ff calc box; size is 3 cells away from boundary condition
start = [mesh.x(12)     mesh.y(12)     mesh.z(12)];
stop  = [mesh.x(end-11) mesh.y(end-11) mesh.z(end-11)];
[CSX, nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', start, stop);

% just to visualize
CSX = AddDump(CSX,'my_nf2ff', 'DumpType', 11, 'Frequency',[9e8]);
CSX = AddBox(CSX,'my_nf2ff',2,start,stop); %assign box

%% Prepare and Run Simulation
Sim_Path = 'tmp_Patch_Ant_simulation';
Sim_CSX = 'patch_ant_simulation.xml';

% create an empty working directory
[status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
[status, message, messageid] = mkdir( Sim_Path ); % create empty simulation folder

% write openEMS compatible xml-file
WriteOpenEMS( [Sim_Path '/' Sim_CSX], FDTD, CSX );

% show the structure
CSXGeomPlot( [Sim_Path '/' Sim_CSX] );

% run openEMS
RunOpenEMS( Sim_Path, Sim_CSX);

%% Postprocessing & Plots
freq = linspace( max([200e6,f0-fc]), f0+fc, 501 );
port = calcPort(port, Sim_Path, freq);

% %% Smith chart port reflection
% plotRefl_new(port, 'threshold', -10)
% title( 'reflection coefficient' );

% plot feed point impedance
Zin = port{1}.uf.tot ./ port{1}.if.tot;
figure
plot( freq/1e6, real(Zin), 'k-', 'Linewidth', 2 );
hold on
grid on
plot( freq/1e6, imag(Zin), 'r--', 'Linewidth', 2 );
title( 'feed point impedance' );
xlabel( 'frequency f / MHz' );
ylabel( 'impedance Z_{in} / Ohm' );
legend( 'real', 'imag' );

% plot reflection coefficient S11
s11 = port{1}.uf.ref ./ port{1}.uf.inc;
figure
plot( freq/1e6, 20*log10(abs(s11)), 'k-', 'Linewidth', 2 );
grid on
title( 'reflection coefficient S_{11}' );
xlabel( 'frequency f / MHz' );
ylabel( 'reflection coefficient |S_{11}|' );

drawnow

%% NFFF Plots
%find resonance frequncy from s11
f_res_ind = find(s11==min(s11));
f_res = freq(f_res_ind);



% calculate the far field at phi=0 degrees and at phi=90 degrees
disp( 'calculating far field at phi=[0 90] deg...' );

nf2ff = CalcNF2FF(nf2ff, Sim_Path, f_res, [-180:2:180]*pi/180, [0 90]*pi/180);

% display power and directivity
disp( ['radiated power: Prad = ' num2str(nf2ff.Prad) ' Watt']);
disp( ['directivity: Dmax = ' num2str(nf2ff.Dmax) ' (' num2str(10*log10(nf2ff.Dmax)) ' dBi)'] );
disp( ['efficiency: nu_rad = ' num2str(100*nf2ff.Prad./port{1}.P_inc(f_res_ind)) ' %']);

% normalized directivity as polar plot
figure
polarFF(nf2ff,'xaxis','theta','param',[1 2],'normalize',1)

% log-scale directivity plot
figure
plotFFdB(nf2ff,'xaxis','theta','param',[1 2])
% conventional plot approach
% plot( nf2ff.theta*180/pi, 20*log10(nf2ff.E_norm{1}/max(nf2ff.E_norm{1}(:)))+10*log10(nf2ff.Dmax));

drawnow

% Show 3D pattern
disp( 'calculating 3D far field pattern and dumping to vtk (use Paraview to visualize)...' );
thetaRange = (0:2:180);
phiRange = (0:2:360) - 180;
nf2ff = CalcNF2FF(nf2ff, Sim_Path, f_res, thetaRange*pi/180, phiRange*pi/180,'Verbose',1,'Outfile','3D_Pattern.h5');

figure
plotFF3D(nf2ff,'logscale',-20);


E_far_normalized = nf2ff.E_norm{1} / max(nf2ff.E_norm{1}(:)) * nf2ff.Dmax;
DumpFF2VTK([Sim_Path '/3D_Pattern.vtk'],E_far_normalized,thetaRange,phiRange,'scale',1e-3);


%% postproc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for my plots

s11 = port{1}.uf.ref./ port{1}.uf.inc;
s21 = port{2}.uf.ref./ port{1}.uf.inc;
%% plot s-parameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% my plots
figure
plot(freq*1e-6,20*log10(abs(s11)),'k-','Linewidth',2);
xlim([freq(1) freq(end)]*1e-6);
grid on;
hold on;
plot(freq*1e-6,20*log10(abs(s21)),'r--','Linewidth',2);
l = legend('S_{11}','S_{21}','Location','Best');
set(l,'FontSize',12);
ylabel('S-Parameter (dB)','FontSize',12);
xlabel('frequency (MHz) \rightarrow','FontSize',12);