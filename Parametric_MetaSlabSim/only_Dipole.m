physical_constants;
unit = 1e-3; % all length in mm. (unit only for the geometry)

%% Setup FDTD Parameter & Excitation Function
if ~exist('input_CSX','var') 
    f0 = 8.5e8; % center frequency (Hz)
    fc = 1e8; % 20 dB corner frequency (Hz) -----> it determines the bandwidth (keep it less than f0)
    FDTD = InitFDTD( 'NrTs', 300000, 'EndCriteria', 1e-5);
    FDTD = SetGaussExcite( FDTD, f0, fc );
    BC = {'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8'}; % boundary conditions
    % FDTD = SetBoundaryCond( FDTD, BC,'PML_Grading','-log(1e-6)*log(2.5)/(2*dl*pow(2.5,W/dl)-1) * pow(2.5, D/dl) / Z' );
    FDTD = SetBoundaryCond( FDTD, BC );
end

%%
SimBox = [300, 300, 300];
dumpWidth = 400;
dumpLength = 400;


    %setup feeding dipole 1
    feed.pos = [0;0;30];      % reference point (center) for the feeding (x,y,z)
    % feed.length = lamda/2;      % length of the dipole
    feed.length = 176.5;      % length of the dipole
    feed.dir = [1,0,0];         % direction of antenna. 100->x, 010->y, 001->z
    feed.R = 50;                % feed resistance
        % calcualte the vertices of the dipole
    feed.start = feed.pos-feed.length/2*feed.dir';
    feed.stop = feed.pos+feed.length/2*feed.dir';

    
%%
CSX = InitCSX();
start = [feed.start];
stop = [feed.stop];
[CSX, port{1}] = AddLumpedPort(CSX, 10 ,1 ,feed.R, start, stop, feed.dir, true);




% Finalize the Mesh
% -----------------
% detect all edges except of the patch
% mesh = DetectEdges(CSX, mesh,'ExcludeProperty',['SRRpatch1', 'SRRpatch2','SRRpatch3','SRRpatch4','patch']);
% detect and set a special 2D metal edge mesh for the patch
% mesh = DetectEdges(CSX, mesh,'SetProperty',['SRRpatch1', 'SRRpatch2','SRRpatch3','SRRpatch4','patch'],'2D_Metal_Edge_Res', c0/(f0+fc)/unit/50);
% generate a smooth mesh with max. cell size: (e.g) lambda_min / 20 
% mesh = SmoothMesh(mesh, c0/(f0+fc)/unit/20); 
mesh = SmoothMesh(mesh, c0/(f0+fc)/unit/30,'algorithm',[1 3]); % alternative. Useful when SmoothMesh gets stuck in an infinite loop inside "SmoothMeshLines2" function
% mesh.x= unique(round(mesh.x,0));
% mesh.y= unique(round(mesh.y,0));

CSX = DefineRectGrid(CSX, unit, mesh);

% add a nf2ff calc box; size is 3 cells away from boundary condition
start = [mesh.x(12)     mesh.y(12)     mesh.z(12)];
stop  = [mesh.x(end-11) mesh.y(end-11) mesh.z(end-11)];
[CSX, nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', start, stop);

%record E-Field
start = [mesh.x(2)     mesh.y(2)     mesh.z(2)];
stop  = [mesh.x(end-1) mesh.y(end-1) mesh.z(end-1)];
CSX = AddDump(CSX,'Ef', 'DumpType', 10, 'Frequency',(f0),'FileType',1);
CSX = AddBox(CSX,'Ef',2,start, stop); %assign box
CSX = AddDump(CSX,'Ef_vtr', 'DumpType', 10, 'Frequency',(f0));
CSX = AddBox(CSX,'Ef_vtr',2,start, stop); %assign box





%% Prepare and Run Simulation
Sim_Path = 'tmp_';
Sim_Path = [Sim_Path, 'Dipole'];
Sim_CSX = [Sim_Path(5:end), '.xml'];

% create an empty working directory
[status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
[status, message, messageid] = mkdir( Sim_Path ); % create empty simulation folder

% write openEMS compatible xml-file
WriteOpenEMS( [Sim_Path '/' Sim_CSX], FDTD, CSX );

% show the structure
CSXGeomPlot( [Sim_Path '/' Sim_CSX] );

% run openEMS
RunOpenEMS( Sim_Path, Sim_CSX);

    
    
    
    