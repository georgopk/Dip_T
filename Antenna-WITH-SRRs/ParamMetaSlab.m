function [freq, port, s11, s21, xy_points, CSX, mesh] = ParamMetaSlab(varargin)
% Simulation of a Slab with two arrays of SRRs
% 
% -Create a slab made of ground and substrate.
%  The ground dimensions and position can be changed parametrically.
% -Calculate the number of SRRs needed to fill the perimeter of the slab.
% -Fill the four edges of the slab with SRRs. (As a patch on the substrate 
%  of the slab). The SRR parameters can be changed parametrically.
%  *** This task uses the function srr_points.m ***
% 
% IF the parameter 'CSX' is not set:
% -Create a dipole (lamda/2)
%   The length of the dipole can be canged parametrically.
% -Create an auto-adjustable Simulation Box depending on the dimensions of
%   the structure and the dipole, accounting for "extra free space".
% -Store the values of the Electric field (frequency domain) both in .vtr
%   and in .h5 files.
% 
% SECOND WAY TO USE THIS FUNCTION:
% --Use the argument 'CSX'--
% The function will return a CSX object containing the slab with SRRs. 
% It also returns a 'mesh' matrix  and a 'xy_points' matrix containing the 
% points of the sructure.
% 
% --In this case the output arguments freq, port, s11, s21 are useless.--
% 
% 
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
% 
% 
% (optional) arguments:
%
%    'L', 'w', 'g', 's', 'd':      SRR parameters.
%    'dipLength':                  Length of the feeding dipole.
%    'grndxDim', 'grndyDim':       Dimensions of the ground.
%    'grndelev':                   Ground elevation. (z position)
%    'CSX':                        A CSX object. If this argument is being 
%                                  used the function creates a ground, 
%                                  fills the edges with SRRs and substrate
%                                  and returns the CSX objecet.
%    'mesh':                       A mesh matrix. Usefull when 'CSX'
%                                  argument is being used to add lines in
%                                  your mesh for each point of the Slab
%                                  structure.
%    'substrate':                  A substrate structure containing epsR,
%                                  tan_delta, kappa, thickness, cells.
%
%   --there are default values for any of the arguments above--    
% 
% 
% 
% 
%  examples:
% -----------
% Ex.1
% Run a simulation with L=40, s=3, dipLength = 150, g = 3, w = 5
% >> ParamMetaSlab('L',40,'s',3,'dipLength',150,'g',3,'w',5)
% 
% Ex.2
% Run a simulation with the default values
% >> ParamMetaSlab()
% 
% Ex.3
% Run a simulation with w=5, L=50 (using the defaults for the other
% variables)
% >> ParamMetaSlab('w',5,'L',50)
% 
% Ex.4
% Call the function using the 'CSX' and 'mesh' arguments, defining L=14
% and grndelev=20;
% >>[~,~,~,~, in_SRR_points, CSX, mesh] = ParamMetaSlab('L',14,'CSX',CSX,'grndelev',20,'mesh',mesh);
% 


%% Store input parameters
for n=1:numel(varargin)/2
    if ~ischar(varargin{2*n-1})
        error([' NOT a variable name: ' varargin{2*n-1}]);
    end
    if strcmp(varargin{2*n-1},'L')
        input_L = varargin{2*n};
    elseif strcmp(varargin{2*n-1},'w')
        input_w = varargin{2*n};
    elseif strcmp(varargin{2*n-1},'g')
        input_g = varargin{2*n};
    elseif strcmp(varargin{2*n-1},'s')
        input_s = varargin{2*n};
    elseif strcmp(varargin{2*n-1},'d')
        input_d = varargin{2*n};
    elseif strcmp(varargin{2*n-1},'dipLength')
        input_dipLength = varargin{2*n};
    elseif strcmp(varargin{2*n-1},'grndxDim')
        input_grndxDim = varargin{2*n};
    elseif strcmp(varargin{2*n-1},'grndyDim')
        input_grndyDim = varargin{2*n};
    elseif strcmp(varargin{2*n-1},'grndelev')
        input_grndelev = varargin{2*n};
    elseif strcmp(varargin{2*n-1},'CSX')
        input_CSX = varargin{2*n};
    elseif strcmp(varargin{2*n-1},'mesh')
        input_mesh = varargin{2*n};
    elseif strcmp(varargin{2*n-1},'substrate')
        input_substrate = varargin{2*n};
    elseif strcmp(varargin{2*n-1},'inv')
        input_inv = varargin{2*n};
    else
        error(['Is this a variable??? ', varargin{2*n-1}]);
    end
end
% 
% % close all
% clear %%%% Important !!!!!!!!!!
% % clc

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

%% Setup the Simulation (Initialize Geometry)


    
physical_constants;
if ~exist('input_CSX','var')                % do not declare these parameters if 'CSX' is being used
    unit = 1e-3; % all length in mm. (unit only for the geometry)
    lamda_meter = c0/f0;
    lamda = lamda_meter/unit;
end

plot_on = 0;    % useful for debuggingd

%piorities
patchPri = 50+5;
substratePri = 10+5;
groundPri = 40+5;
feedPri = 30+5;

%inverse SRRs
inv = 0; 
if exist('input_inv','var')       % if exists, overwrite the default value
    inv = input_inv;
end

%ground setup
grnd.pos = [0;0;0]; % reference point for the ground (x,y)
grnd.thickness = 0;
grnd.x_dim = 370.8;
grnd.y_dim = 370.8;
if exist('input_grndxDim','var')       % if exists, overwrite the default value
    grnd.x_dim = input_grndxDim;
end
if exist('input_grndyDim','var')       % if exists, overwrite the default value
    grnd.y_dim = input_grndyDim;
end
if exist('input_grndelev','var')       % if exists, overwrite the default value
    grnd.pos(3) = input_grndelev;
end
grnd.points = [-grnd.x_dim/2,-grnd.y_dim/2; 
    -grnd.x_dim/2,grnd.y_dim/2; 
    grnd.x_dim/2,grnd.y_dim/2; 
    grnd.x_dim/2,-grnd.y_dim/2]' ;
grnd.points = grnd.points + grnd.pos(1:2);


if ~exist('input_CSX','var')  % do not declare these parameters if 'CSX' is being used
    %setup feeding dipole 1
    feed.pos = [0;0;30];      % reference point (center) for the feeding (x,y,z)
    % feed.length = lamda/2;      % length of the dipole
    feed.length = 176.5;      % length of the dipole
    if exist('input_dipLength','var')       % if exists, overwrite the default value
        feed.length = input_dipLength;
    end
    feed.dir = [1,0,0];         % direction of antenna. 100->x, 010->y, 001->z
    feed.R = 50;                % feed resistance
        % calcualte the vertices of the dipole
    feed.start = feed.pos-feed.length/2*feed.dir';
    feed.stop = feed.pos+feed.length/2*feed.dir';

    %setup "hearing" dipole 2
    feed2.pos = [0;0;-feed.pos(3)];      % reference point (center) for the feeding (x,y,z)
    % feed.length = lamda/2;      % length of the dipole
    feed2.length = feed.length;      % length of the dipole
    feed2.dir = feed.dir;         % direction of antenna. 100->x, 010->y, 001->z
    feed2.R = feed.R;                % feed resistance
        % calcualte the vertices of the dipole
    feed2.start = feed2.pos-feed2.length/2*feed2.dir';
    feed2.stop = feed2.pos+feed2.length/2*feed2.dir';
end



% SRR/metamaterial setup
srr.L = 50;     % length of the outer edge
srr.w = 2;      % width of the strip
srr.g = 2;      % gap between strips
srr.s = 2;      % slot dimension
srr.d = 3;      % distance between two SRRs
if exist('input_L','var')       % if exists, overwrite the default value
    srr.L = input_L;
end
if exist('input_w','var')       % if exists, overwrite the default value
    srr.w = input_w;
end
if exist('input_g','var')       % if exists, overwrite the default value
    srr.g = input_g;
end
if exist('input_s','var')       % if exists, overwrite the default value
    srr.s = input_s;
end
if exist('input_d','var')       % if exists, overwrite the default value
    srr.d = input_d;
end

% calculate parameters for SRRs (on X-axis)
srr.d_c = srr.L + srr.d;  % distance between the centers of two (neighbouring) SRRs
srr.nCellsX = floor (grnd.x_dim/srr.d_c); % Number of SRR cells
srr.first_centerX = -(srr.nCellsX-1)*srr.d_c/2;
srr.centersX = srr.first_centerX + (0:(srr.nCellsX-1)) * srr.d_c;
% calculate parameters for SRRs (on Y-axis)
srr.nCellsY = floor (grnd.y_dim/srr.d_c)-2; % Number of SRR cells
srr.first_centerY = -(srr.nCellsY-1)*srr.d_c/2;
srr.centersY = srr.first_centerY + (0:(srr.nCellsY-1)) * srr.d_c;


%substrate constants
    sub_freq = 8.5e8; % Frequency to calculate the substrate conductivity for
    substrate.epsR = 4.2;
    substrate.tan_delta = 0.025;
    substrate.kappa  = substrate.tan_delta * 2*pi*sub_freq * EPS0*substrate.epsR; %conductivity 
    substrate.thickness = 1.524;
    if(inv==1)
    substrate.thickness=-substrate.thickness;
    end
    substrate.cells = 4;
if exist('input_substrate','var') % if exists, overwrite the default value
    substrate=input_substrate;
end

%substrate 1 setup
substrate1.pos = [grnd.pos(1),grnd.pos(2)-grnd.y_dim/2+srr.L/2+1, grnd.pos(3)];
substrate1.x_dim = grnd.x_dim;
substrate1.y_dim = srr.L+2;
substrate1.points = [-substrate1.x_dim/2,-substrate1.y_dim/2; 
    -substrate1.x_dim/2,substrate1.y_dim/2; 
    substrate1.x_dim/2,substrate1.y_dim/2; 
    substrate1.x_dim/2,-substrate1.y_dim/2]' ;
substrate1.points = substrate1.points + substrate1.pos(1:2)';
%substrate 2 setup
substrate2.pos = [grnd.pos(1),grnd.pos(2)+grnd.y_dim/2-srr.L/2-1, grnd.pos(3)];
substrate2.x_dim = grnd.x_dim;
substrate2.y_dim = srr.L+2;
substrate2.points = [-substrate2.x_dim/2,-substrate2.y_dim/2; 
    -substrate2.x_dim/2,substrate2.y_dim/2; 
    substrate2.x_dim/2,substrate2.y_dim/2; 
    substrate2.x_dim/2,-substrate2.y_dim/2]' ;
substrate2.points = substrate2.points + substrate2.pos(1:2)';

%substrate 3 setup
substrate3.pos = [grnd.pos(1)+grnd.x_dim/2-srr.L/2-1,grnd.pos(2), grnd.pos(3)];
substrate3.x_dim = srr.L+2;
substrate3.y_dim = grnd.y_dim;
substrate3.points = [-substrate3.x_dim/2,-substrate3.y_dim/2; 
    -substrate3.x_dim/2,substrate3.y_dim/2; 
    substrate3.x_dim/2,substrate3.y_dim/2; 
    substrate3.x_dim/2,-substrate3.y_dim/2]' ;
substrate3.points = substrate3.points + substrate3.pos(1:2)';
%substrate 4 setup
substrate4.pos = [grnd.pos(1)-grnd.x_dim/2+srr.L/2+1,grnd.pos(2), grnd.pos(3)];
substrate4.x_dim = srr.L+2;
substrate4.y_dim = grnd.y_dim;
substrate4.points = [-substrate4.x_dim/2,-substrate4.y_dim/2; 
    -substrate4.x_dim/2,substrate4.y_dim/2; 
    substrate4.x_dim/2,substrate4.y_dim/2; 
    substrate4.x_dim/2,-substrate4.y_dim/2]' ;
substrate4.points = substrate4.points + substrate4.pos(1:2)';

% size of the simulation box 
% aggregate all coordinates
all_xy = [grnd.points , substrate1.points,substrate2.points]' ;
if ~exist('input_CSX','var'); all_xy =[all_xy', [feed.start(1:2), feed2.stop(1:2)]]'; end
all_z = [grnd.pos(3), substrate1.pos(3), substrate2.pos(3)]';
if ~exist('input_CSX','var'); all_z = [all_z', feed.start(3), feed2.start(3)]'; end
% find "outlier" coordinates
StructDimen(1,:) = [min(all_xy(:,1)), max(all_xy(:,1))];
StructDimen(2,:) = [min(all_xy(:,2)), max(all_xy(:,2))];
StructDimen(3,:) = [min(all_z), max(all_z)];



%% Setup CSXCAD Geometry & Mesh
CSX = InitCSX();
if exist('input_CSX','var')       % if exists, overwrite the default value
   CSX = input_CSX;
   freq = 0;        %just to assign a value
   port = 0;        %just to assign a value
   s11 = 0;         %just to assign a value
   s21= 0;          %just to assign a value
end

%initialize the mesh with the "air-box" dimensions
if ~exist('input_CSX','var') 
    mesh.x = [StructDimen(1,1) - lamda/6, StructDimen(1,2) + lamda/6];
    mesh.y = [StructDimen(2,1) - lamda/6, StructDimen(2,2) + lamda/6];
    mesh.z = [StructDimen(3,1) - lamda/6, StructDimen(3,2) + lamda/6];
end
if exist('input_mesh','var')       % if exists, overwrite the default value
    mesh = input_mesh;
end


[points1,points2] = srr_points(srr.L,srr.w,srr.g,srr.s,plot_on);               % calculate points for an SRR
all_SRR_points =[];
 
% SRRs on X-axis
CSX = AddMetal( CSX, ['SRRpatch1',num2str(inv)] );                     % create a perfect electric conductor (PEC) named "patch"
 for i=1:srr.nCellsX
    % Create outer strip    
    tmp_SRR_points = [points1' + [srr.centersX(i);substrate1.pos(2)] ];
    CSX = AddPolygon(CSX,['SRRpatch1',num2str(inv)],patchPri,2,substrate1.pos(3) + substrate.thickness,tmp_SRR_points);  % create a polygon of the material "patch"
    all_SRR_points = [all_SRR_points, tmp_SRR_points]; % store SRR points
    %  Create inner strip
    tmp_SRR_points = [ points2' + [srr.centersX(i);substrate1.pos(2)]];
    CSX = AddPolygon(CSX,['SRRpatch1',num2str(inv)],patchPri+1,2,substrate1.pos(3)+substrate.thickness,tmp_SRR_points);    % create a polygon of the material "patch"
    all_SRR_points = [all_SRR_points, tmp_SRR_points]; % store SRR points
 end
CSX = AddMetal( CSX, ['SRRpatch2',num2str(inv)] );                     % create a perfect electric conductor (PEC) named "patch"
for i=1:srr.nCellsX
    % Create outer strip 2
    tmp_SRR_points = points1' + [srr.centersX(i);substrate2.pos(2)];
    CSX = AddPolygon(CSX,['SRRpatch2',num2str(inv)],patchPri,2,substrate2.pos(3)+substrate.thickness,tmp_SRR_points);  % create a polygon of the material "patch"
    all_SRR_points = [all_SRR_points, tmp_SRR_points]; % store SRR points
    %  Create inner strip 2
    tmp_SRR_points = points2' + [srr.centersX(i);substrate2.pos(2)];
    CSX = AddPolygon(CSX,['SRRpatch2',num2str(inv)],patchPri+1,2,substrate2.pos(3)+substrate.thickness,tmp_SRR_points );    % create a polygon of the material "patch"
    all_SRR_points = [all_SRR_points, tmp_SRR_points]; % store SRR points
end

% SRRs on Y-axis 
points1 = sp_round(rotate_points(points1',pi/2),0.2)';
points2 = sp_round(rotate_points(points2',pi/2),0.2)';
points1=round(points1,1);
points2=round(points2,1);
CSX = AddMetal( CSX, ['SRRpatch3',num2str(inv)] );                     % create a perfect electric conductor (PEC) named "patch"
for i=1:srr.nCellsY
    % Create outer strip    
    tmp_SRR_points = [points1' + [substrate3.pos(1);srr.centersY(i)] ];
    CSX = AddPolygon(CSX,['SRRpatch3',num2str(inv)],patchPri,2,substrate3.pos(3) + substrate.thickness,tmp_SRR_points);  % create a polygon of the material "patch"
    all_SRR_points = [all_SRR_points, tmp_SRR_points]; % store SRR points
    %  Create inner strip
    tmp_SRR_points = [ points2' + [substrate3.pos(1);srr.centersY(i)]];
    CSX = AddPolygon(CSX,['SRRpatch3',num2str(inv)],patchPri+1,2,substrate3.pos(3)+substrate.thickness,tmp_SRR_points);    % create a polygon of the material "patch"
    all_SRR_points = [all_SRR_points, tmp_SRR_points]; % store SRR points
 end
CSX = AddMetal( CSX, ['SRRpatch4',num2str(inv)] );                     % create a perfect electric conductor (PEC) named "patch"
for i=1:srr.nCellsY
    % Create outer strip    
    tmp_SRR_points = [points1' + [substrate4.pos(1);srr.centersY(i)] ];
    CSX = AddPolygon(CSX,['SRRpatch4',num2str(inv)],patchPri,2,substrate4.pos(3) + substrate.thickness,tmp_SRR_points);  % create a polygon of the material "patch"
    all_SRR_points = [all_SRR_points, tmp_SRR_points]; % store SRR points
    %  Create inner strip
    tmp_SRR_points = [ points2' + [substrate4.pos(1);srr.centersY(i)]];
    CSX = AddPolygon(CSX,['SRRpatch4',num2str(inv)],patchPri+1,2,substrate4.pos(3)+substrate.thickness,tmp_SRR_points);    % create a polygon of the material "patch"
    all_SRR_points = [all_SRR_points, tmp_SRR_points]; % store SRR points
end
 


xy_points = [all_xy', all_SRR_points];
  
% Create Substrate 1
CSX = AddMaterial( CSX, ['SRRsub',num2str(inv)] );                                                              % create a new material named "substrate"
CSX = SetMaterialProperty( CSX, ['SRRsub',num2str(inv)], 'Epsilon', substrate.epsR, 'Kappa', substrate.kappa ); % define the properties of the material "substrate"
CSX = AddLinPoly(CSX, ['SRRsub',num2str(inv)], substratePri, 2, substrate1.pos(3), substrate1.points, substrate.thickness);     % create a polygon of the material "substrate" with thickness
% add extra cells to discretize the substrate thickness
mesh.z = [linspace(substrate1.pos(3),substrate1.pos(3)+ substrate.thickness,substrate.cells+1) mesh.z];
% Create Substrate 2
% CSX = AddMaterial( CSX, 'SRRsub' );                                                              % create a new material named "substrate"
CSX = SetMaterialProperty( CSX, ['SRRsub',num2str(inv)], 'Epsilon', substrate.epsR, 'Kappa', substrate.kappa ); % define the properties of the material "substrate"
CSX = AddLinPoly(CSX, ['SRRsub',num2str(inv)], substratePri, 2, substrate2.pos(3), substrate2.points, substrate.thickness);     % create a polygon of the material "substrate" with thickness
% add extra cells to discretize the substrate thickness
mesh.z = [linspace(substrate2.pos(3),substrate2.pos(3)+ substrate.thickness,substrate.cells+1) mesh.z];

% Create Substrate 3
% CSX = AddMaterial( CSX, 'SRRsub' );                                                              % create a new material named "substrate"
CSX = SetMaterialProperty( CSX, ['SRRsub',num2str(inv)], 'Epsilon', substrate.epsR, 'Kappa', substrate.kappa ); % define the properties of the material "substrate"
CSX = AddLinPoly(CSX, ['SRRsub',num2str(inv)], substratePri, 2, substrate3.pos(3), substrate3.points, substrate.thickness);     % create a polygon of the material "substrate" with thickness
% add extra cells to discretize the substrate thickness
mesh.z = [linspace(substrate3.pos(3),substrate3.pos(3)+ substrate.thickness,substrate.cells+1) mesh.z];
% Create Substrate 4
% CSX = AddMaterial( CSX, 'SRRsub' );                                                              % create a new material named "substrate"
CSX = SetMaterialProperty( CSX, ['SRRsub',num2str(inv)], 'Epsilon', substrate.epsR, 'Kappa', substrate.kappa ); % define the properties of the material "substrate"
CSX = AddLinPoly(CSX, ['SRRsub',num2str(inv)], substratePri, 2, substrate4.pos(3), substrate4.points, substrate.thickness);     % create a polygon of the material "substrate" with thickness
% add extra cells to discretize the substrate thickness
mesh.z = [linspace(substrate4.pos(3),substrate4.pos(3)+ substrate.thickness,substrate.cells+1) mesh.z];

% Create Ground
CSX = AddMetal( CSX, ['SRRgnd',num2str(inv)] );                           % create a perfect electric conductor (PEC) named "SRRgnd"
CSX = AddLinPoly(CSX, ['SRRgnd',num2str(inv)], groundPri, 2, grnd.pos(3),grnd.points, grnd.thickness);   % create a polygon of the material "SRRgnd"

if ~exist('input_CSX','var') 
% --- Add dipole 1 ---
[CSX, port{1}] = AddLumpedPort(CSX, feedPri ,1 ,feed.R, feed.start, feed.stop, feed.dir, true);

% % --- Add dipole 2 ---
[CSX, port{2}] = AddLumpedPort(CSX, feedPri ,2 ,feed2.R, feed2.start, feed2.stop, feed2.dir);
end

if ~exist('input_CSX','var')       % if the function is used to return the CSX to an other function, DO NOT execute the following lines


        % % detect all edges except of the patch
        mesh = DetectEdges(CSX, mesh,'ExcludeProperty','patch');
        % mesh = DetectEdges(CSX, mesh);
        % detect and set a special 2D metal edge mesh for the patch
        mesh = DetectEdges(CSX, mesh,'SetProperty','patch','2D_Metal_Edge_Res', c0/(f0+fc)/unit/50);
        % generate a smooth mesh with max. cell size: (e.g) lambda_min / 20 
        % mesh = SmoothMesh(mesh, c0/(f0+fc)/unit/20); 
        mesh = SmoothMesh(mesh, c0/(f0+fc)/unit/65,'algorithm',[1 3]); % alternative. Useful when SmoothMesh gets stuck in an infinite loop inside "SmoothMeshLines2" function

        CSX = DefineRectGrid(CSX, unit, mesh);

        % % --- Add port ---
        % start = [-substrate.points(1,1), substrate.points(2,1), substrate.pos(3) - substrate.thickness/4 ];
        % stop = [substrate.points(1,1)+0.5, substrate.points(2,1)+0.5, substrate.thickness];
        % [CSX, port{2}] = AddRectWaveGuidePort(CSX,0,2,start,stop, 'y',(stop(2)-start(2)),(stop(3)-start(3)),'TE11', 0);


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
        Sim_Path = ['tmp_Slab_sim' '_s' num2str(srr.s) '_g' num2str(srr.g) '_L' num2str(srr.L) '_w' num2str(srr.w)];
        Sim_CSX = 'Slab_simulation.xml';

        % create an empty working directory
        [status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
        [status, message, messageid] = mkdir( Sim_Path ); % create empty simulation folder

        % write openEMS compatible xml-file
        WriteOpenEMS( [Sim_Path '/' Sim_CSX], FDTD, CSX );

        % % show the structure
        CSXGeomPlot( [Sim_Path '/' Sim_CSX] );

        % run openEMS
        RunOpenEMS( Sim_Path, Sim_CSX);


        %% Post-processing

        freq = linspace( f0-fc, f0+fc, 501 );
        port = calcPort(port, Sim_Path, freq);

        Zin = port{1}.uf.tot ./ port{1}.if.tot;
        s11 = port{1}.uf.ref ./ port{1}.uf.inc;
        s21 = port{2}.uf.ref ./ port{1}.uf.inc;
        P_in = 0.5 * port{1}.uf.inc .* conj( port{1}.if.inc ); % antenna feed power
        % % % 
        % plot feed point impedance
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
        figure
        plot( freq/1e6, 20*log10(abs(s11)), 'k-', 'Linewidth', 2 );
        grid on
        title( 'reflection coefficient S_{11}' );
        xlabel( 'frequency f / MHz' );
        ylabel( 'reflection coefficient |S_{11}|' );
        
        % plot coefficient S21
        figure
        plot( freq/1e6, 20*log10(abs(s21)), 'r', 'Linewidth', 2 );
        grid on
        title( 'S_{21}' );
        xlabel( 'frequency f / MHz' );
        ylabel( '|S_{21}|' );



end
end