% Simulation of a Slab with an array of SRRs
% 
% -Create a slab made of ground and substrate.
%   The position of the slab can be arbitrary.
% -Calculate the number of SRRs needed to fill the slab (on x axis)
% -Fill the slab with SRRs. (As a patch on the substrate of the slab)
%   (This task uses the function srr_points.m)
% -Create a dipol (lamda/2)
%   The position of the dipol can be arbitrary.
% -Create an auto-adjustable Simulation Box depending on the dimentions of
%   the structure and the dipol, accounting for "extra lamda/4 free space".
% -Store the values of the Electric field (frequency domain) both in .vtr
%   and in .h5 files.
% -Plot the magnitude of the lying on the edge of the slab electric field.
% 
%  ----------------- NOTE: -------------------
%  This Simulation needs the file srr_points.m
%  -------------------------------------------
% 
% SRR parameters explanation:
%        ________________________                 ________________________                 
%   ^   |                        |               |                        |                
%   !   |    ________________    |               |    ________________    |                
%   !   |   |                |   |               |   |                |   |                
%   !   |   | g  ________    |   |               |   | g  ________    |   |                
%   !   |   |<->|  ____  |   |___|    _          |   |<->|  ____  |   |___|     _          
%   !   |   |   |_|    | |             \         |   |   |_|    | |              \         
%   L   |   | S{ _     | |              }-> S    |   | S{ _     | |               }-> S    
%   !   |   |   | |____| |  ->___<-w  _/         |   |   | |____| |  ->___<-w   _/         
%   !   |   |   |________|   |   |               |   |   |________|   |   |                
%   !   |   |                |   |               |   |                |   |                
%   !   |   |________________|   | <==== d ====> |   |________________|   |                
%   !   |                        |               |                        |                
%   -   |________________________|               |________________________|                
%        <--------- L ---------->                 <--------- L ---------->                 



% close all
clear
% clc

%% Setup FDTD Parameter & Excitation Function
f0 = 8.5e8; % center frequency (Hz)
fc = 0.5e8; % 20 dB corner frequency (Hz) -----> it determines the bandwidth (keep it less than f0)
FDTD = InitFDTD( 'NrTs', 300000, 'EndCriteria', 1e-5);
FDTD = SetGaussExcite( FDTD, f0, fc );
BC = {'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8'}; % boundary conditions
% FDTD = SetBoundaryCond( FDTD, BC,'PML_Grading','-log(1e-6)*log(2.5)/(2*dl*pow(2.5,W/dl)-1) * pow(2.5, D/dl) / Z' );
FDTD = SetBoundaryCond( FDTD, BC );

%% Setup the Simulation (Initialize Geometry)
physical_constants;
unit = 1e-3; % all length in mm. (unit only for the geometry)
lamda_meter = c0/f0;
lamda = lamda_meter/unit;
% metamaterial_on = 1;

plot_on = 1;    % useful for debugging

%piorities
patchPri = 50;
substratePri = 10;
groundPri = 40;
feedPri = 30;

%ground setup
grnd.pos = [0;0;0]; % reference point for the ground (x,y)
grnd.thickness = 0;
grnd.x_dim = 370.8;
% grnd.x_dim = 50;
grnd.y_dim = 50;
grnd.points = [-grnd.x_dim/2,-grnd.y_dim/2; 
    -grnd.x_dim/2,grnd.y_dim/2; 
    grnd.x_dim/2,grnd.y_dim/2; 
    grnd.x_dim/2,-grnd.y_dim/2]' ;
grnd.points = grnd.points + grnd.pos(1:2);

%substrate setup
sub_freq = 8.5e8; % Frequency to calculate the substrate conductivity for
substrate.epsR = 4.2;
substrate.tan_delta = 0.025;
substrate.kappa  = substrate.tan_delta * 2*pi*sub_freq * EPS0*substrate.epsR; %conductivity 
substrate.thickness = 1.524;
substrate.cells = 4;
substrate.pos = grnd.pos;
substrate.x_dim = grnd.x_dim;
substrate.y_dim = grnd.y_dim;
substrate.points = [-substrate.x_dim/2,-substrate.y_dim/2; 
    -substrate.x_dim/2,substrate.y_dim/2; 
    substrate.x_dim/2,substrate.y_dim/2; 
    substrate.x_dim/2,-substrate.y_dim/2]' ;
substrate.points = substrate.points + substrate.pos(1:2);

% SRR/metamaterial setup
srr.L = 20;     % length of the outer edge
srr.w = 1;      % width of the strip
srr.g = 2;      % gap between strips
srr.s = 1;      % slot dimention
srr.d = 2;      % distance between two SRRs
srr.d_c = srr.L + srr.d;  % distance betwenn the centers of two (neighbouring) SRRs
srr.nCells = floor (substrate.x_dim/srr.d_c); % Number of SRR cells
srr.first_center = -(srr.nCells-1)*srr.d_c/2;
srr.centers = srr.first_center + (0:(srr.nCells-1)) * srr.d_c;

%setup feeding
feed.pos = [0;150;20];      % reference point (center) for the feeding (x,y,z)
feed.length = lamda/2;      % length of the dipol
feed.dir = [1,0,0];         % direction of antenna. 100->x, 010->y, 001->z
feed.R = 50;                % feed resistance
    % calcualte the vertices of the dipol
feed.start = feed.pos-feed.length/2*feed.dir';
feed.stop = feed.pos+feed.length/2*feed.dir';

% size of the simulation box 
% aggregate all coordinates
all_xy = [grnd.points , substrate.points, feed.start(1:2), feed.stop(1:2)]';
all_z = [grnd.pos(3), substrate.pos(3), feed.start(3), feed.start(3)]';
% find "outlier" coordinates
StructDimen(1,:) = [min(all_xy(:,1)), max(all_xy(:,1))];
StructDimen(2,:) = [min(all_xy(:,2)), max(all_xy(:,2))];
StructDimen(3,:) = [min(all_z), max(all_z)];


%% Setup CSXCAD Geometry & Mesh
CSX = InitCSX();

%initialize the mesh with the "air-box" dimensions
mesh.x = [StructDimen(1,1) - 2/5*lamda, StructDimen(1,2) + 2/5*lamda];
mesh.y = [StructDimen(2,1) - 2/5*lamda, StructDimen(2,2) + 2/5*lamda];
mesh.z = [StructDimen(3,1) - 2/5*lamda, StructDimen(3,2) + 2/5*lamda];


[points1,points2] = srr_points(srr.L,srr.w,srr.g,srr.s,plot_on);               % calculate points for an SRR
 for i=1:srr.nCells
    % Create outer strip
    CSX = AddMetal( CSX, 'patch' );                     % create a perfect electric conductor (PEC) named "patch"
    CSX = AddPolygon(CSX,'patch',patchPri,2,substrate.pos(3)+substrate.thickness,points1' + [srr.centers(i);0]);  % create a polygon of the material "patch"
    %  Create inner strip
    CSX = AddPolygon(CSX,'patch',patchPri+1,2,substrate.pos(3)+substrate.thickness, points2' + [srr.centers(i);0]);    % create a polygon of the material "patch"
 end

% Create Substrate
CSX = AddMaterial( CSX, 'substrate' );                                                              % create a new material named "substrate"
CSX = SetMaterialProperty( CSX, 'substrate', 'Epsilon', substrate.epsR, 'Kappa', substrate.kappa ); % define the properties of the material "substrate"
CSX = AddLinPoly(CSX, 'substrate', substratePri, 2, substrate.pos(3), substrate.points, substrate.thickness);     % create a polygon of the material "substrate" with thickness
% add extra cells to discretize the substrate thickness
mesh.z = [linspace(substrate.pos(3),substrate.pos(3)+ substrate.thickness,substrate.cells+1) mesh.z];

% Create Ground
CSX = AddMetal( CSX, 'gnd' );                           % create a perfect electric conductor (PEC) named "gnd"
CSX = AddLinPoly(CSX, 'gnd', groundPri, 2, grnd.pos(3),grnd.points, grnd.thickness);   % create a polygon of the material "gnd"

% --- Apply the Excitation ---
[CSX, port{1}] = AddLumpedPort(CSX, feedPri ,1 ,feed.R, feed.start, feed.stop, feed.dir, true);

% % detect all edges except of the patch
mesh = DetectEdges(CSX, mesh,'ExcludeProperty','patch');
% mesh = DetectEdges(CSX, mesh);
% detect and set a special 2D metal edge mesh for the patch
mesh = DetectEdges(CSX, mesh,'SetProperty','patch','2D_Metal_Edge_Res', c0/(f0+fc)/unit/50);
% generate a smooth mesh with max. cell size: (e.g) lambda_min / 20 
% mesh = SmoothMesh(mesh, c0/(f0+fc)/unit/20); 
mesh = SmoothMesh(mesh, c0/(f0+fc)/unit/65,'algorithm',[1 3]); % alternative. Useful when SmoothMesh gets stuck in an infinite loop inside "SmoothMeshLines2" function

CSX = DefineRectGrid(CSX, unit, mesh);


% add a nf2ff calc box; size is 3 cells away from boundary condition
start = [mesh.x(12)     mesh.y(12)     mesh.z(12)];
stop  = [mesh.x(end-11) mesh.y(end-11) mesh.z(end-11)];
[CSX, nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', start, stop);

%record E-Field
CSX = AddDump(CSX,'Ef', 'DumpType', 10, 'Frequency',(f0),'FileType',1);
CSX = AddBox(CSX,'Ef',2,start, stop); %assign box
CSX = AddDump(CSX,'Ef_vtr', 'DumpType', 10, 'Frequency',(f0));
CSX = AddBox(CSX,'Ef_vtr',2,start, stop); %assign box

%% Prepare and Run Simulation
Sim_Path = ['tmp_Slab_sim' '_g' num2str(srr.g) '_L' num2str(srr.L) '_s' num2str(srr.s) '_w' num2str(srr.w)];
Sim_CSX = 'Slab_simulation.xml';

% create an empty working directory
[status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
[status, message, messageid] = mkdir( Sim_Path ); % create empty simulation folder

% write openEMS compatible xml-file
WriteOpenEMS( [Sim_Path '/' Sim_CSX], FDTD, CSX );

% show the structure
CSXGeomPlot( [Sim_Path '/' Sim_CSX] );

% run openEMS
RunOpenEMS( Sim_Path, Sim_CSX);


%% Post-processing

[h5field, h5mesh] = ReadHDF5Dump([Sim_Path '/Ef.h5']);
h5ef = squeeze(h5field.FD.values{1}(:,17,18,:));
h5ef_y = abs(h5ef(:,2));
plot(h5mesh.lines{1}/unit,h5ef_y);
