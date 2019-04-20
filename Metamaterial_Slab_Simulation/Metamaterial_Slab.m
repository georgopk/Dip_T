





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

%piorities
substratePri = 10;
groundPri = 40;
feedPri = 30;

%ground setup
grnd.pos = [0;0;0]; % reference point for the ground (x,y)
grnd.thickness = 0;
grnd.x_dim = 370.8;
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
substrate.x_dim = 370.8;
substrate.y_dim = 50;
substrate.points = [-substrate.x_dim/2,-substrate.y_dim/2; 
    -substrate.x_dim/2,substrate.y_dim/2; 
    substrate.x_dim/2,substrate.y_dim/2; 
    substrate.x_dim/2,-substrate.y_dim/2]' ;
substrate.points = substrate.points + substrate.pos(1:2);

%setup feeding
feed.pos = [0;0;20];         % reference point (center) for the feeding (x,y,z)
feed.length = lamda/2; % length of the dipol
feed.dir = [1,0,0];               % direction of antenna. 0->x, 1->y, 2->z
feed.R = 50;                % feed resistance
    % calcualte the vertices of the dipol
feed.start = feed.pos-feed.length/2*feed.dir';
feed.stop = feed.pos+feed.length/2*feed.dir';

% % calculeate the center of the geometry
% center = ;

% size of the simulation and dump box 
SimBox = [2*abs(substrate.pos(1)) + substrate.x_dim + 2*abs(feed.pos(1)) + 4*lamda/5,
    2*abs(substrate.pos(2)) + substrate.y_dim + 2*abs(feed.pos(2)) + 4*lamda/5,
    2*abs(feed.pos(3)) + 2*abs(feed.pos(3)) + 4*lamda/5 ];

%% Setup CSXCAD Geometry & Mesh
CSX = InitCSX();

%initialize the mesh with the "air-box" dimensions
mesh.x = [-SimBox(1)/2 SimBox(1)/2];
mesh.y = [-SimBox(2)/2 SimBox(2)/2];
mesh.z = [-SimBox(3)/2 SimBox(3)/2];


% Create Substrate
CSX = AddMaterial( CSX, 'substrate' );                                                              % create a new material named "substrate"
CSX = SetMaterialProperty( CSX, 'substrate', 'Epsilon', substrate.epsR, 'Kappa', substrate.kappa ); % define the properties of the material "substrate"
CSX = AddLinPoly(CSX, 'substrate', substratePri, 2, substrate.pos(3), substrate.points, substrate.thickness);     % create a polygon of the material "substrate" with thickness
% add extra cells to discretize the substrate thickness
mesh.z = [linspace(substrate.pos(1),substrate.pos(1)+ substrate.thickness,substrate.cells+1) mesh.z];

% Create Ground
CSX = AddMetal( CSX, 'gnd' );                           % create a perfect electric conductor (PEC) named "gnd"
CSX = AddLinPoly(CSX, 'gnd', groundPri, 2, grnd.pos(3),grnd.points, grnd.thickness);   % create a polygon of the material "gnd"

% --- Apply the Excitation ---
[CSX, port{1}] = AddLumpedPort(CSX, feedPri ,1 ,feed.R, feed.start, feed.stop, feed.dir, true);

% % detect all edges except of the patch
% mesh = DetectEdges(CSX, mesh,'ExcludeProperty','patch');
mesh = DetectEdges(CSX, mesh);
% % detect and set a special 2D metal edge mesh for the patch
% mesh = DetectEdges(CSX, mesh,'SetProperty','patch','2D_Metal_Edge_Res', c0/(f0+fc)/unit/50);
% generate a smooth mesh with max. cell size: (e.g) lambda_min / 20 
% mesh = SmoothMesh(mesh, c0/(f0+fc)/unit/20); 
mesh = SmoothMesh(mesh, c0/(f0+fc)/unit/65,'algorithm',[1 3]); % alternative. Useful when SmoothMesh gets stuck in an infinite loop inside "SmoothMeshLines2" function

CSX = DefineRectGrid(CSX, unit, mesh);


% add a nf2ff calc box; size is 3 cells away from boundary condition
start = [mesh.x(12)     mesh.y(12)     mesh.z(12)];
stop  = [mesh.x(end-11) mesh.y(end-11) mesh.z(end-11)];
[CSX, nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', start, stop);

%record E-Field
CSX = AddDump(CSX,'Ef', 'DumpType', 10, 'Frequency',(f0));
CSX = AddBox(CSX,'Ef',2,start, stop); %assign box

%% Prepare and Run Simulation
Sim_Path = 'tmp_Slab_simulation';
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


