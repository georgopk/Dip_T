%% Antenna coupling (based on...) % Simple Patch Antenna Tutorial
% Simulates one patch Antenna or Two Patch Antennas (back to back).
%-------------------------
%
%  ### SINGLE ANTENNA MODE ### 
% If sec_antenna != 1 (e.g. 0 )
%
% Creates:
%  -a patch antenna
%  -a Rectangle substrate
%  -a Rectangle ground plane
%
% Calculates:
% - the feed point impedance
% - the frequency of the minimum impedance
% (for the frequency range [f0-fc , f0+fc] ) 
% (a lower frequency boundary is set)
%     
% Calculates:
% - the radiation pattern at phi = 0 and phi = 90
% - the 3d representation of the radiation pattern
% (for the frequency of minimum impedance )
% 
%   NOTE:
% In "single antenna mode" it calculates the geometry of the second antenna
% but it does not create this geometry.
% Do NOT remove the lines for the second antenna!!! It causes erros!
% 
% 
%
% ### DOUBLE ANTENNA MODE (BACK TO BACK) ###
% If sec_antenna == 1
% 
% -Creates two identical patch antennas (back to back). 
%   The distance between the antennas is "backDist" (in mm)
% -Creates a simulation Box based on "backDist"
% -rotates the geometry.
%   The rotation angle is "rot" (in rad)
% -"rounds" the coordinates of the points to reduce the mesh
%   keep the coordinates of the ports unchanged and change 
%   other coordinates to "fit" on ports (if necessary)
%
% 
% ### SRRs MODE ###
% If add_srrs == 1
% 
% -Replaces the Ground of the antenna with a new one containing substrate
% with SRRs.
% *** This task uses ParamMetaSlab.m and srr_points.m ***
% To run the simulation without SRRS set 'add_srrs' to 0 
% 
% 
% ---- results ----
% plot S11, S21...
%
% Uses sp_round.m, points.mat, rotate_points.m, ParamMetaSlab.m,
% srr_points.m
%
% 
%  ------ Previous Version ------
%    Complete_Patch_antenna.m
% 
% == based on Simple Patch Antenna Tutorial ==
% %  -------Tutorial info----------
% % Describtion at:
% % <http://openems.de/index.php/Tutorial:_Simple_Patch_Antenna>
% %
% % Tested with
% %  - Matlab 2013a / Octave 4.0
% %  - openEMS v0.0.35
% %
% % (C) 2010-2017 Thorsten Liebig <thorsten.liebig@uni-due.de>
% 



%% Initialization

% close all
clear
% clc

%% Setup the Simulation (Initialize Geometry)
physical_constants;
unit = 1e-3; % all length in mm. (unit only for the geometry)

sec_antenna = 1;
add_srrs = 1;
same_grnd_size = 1;
on_server = 0;

%rotation (rotation of the geometry on XY plane. Useful to "fit" a geometry
%on the grid lines)
% rot = 1.5*2*pi/16;        % 1.5 * (16gon_central_angle)
rot = 0;

%piorities
patchPri = 50;
substratePri = 10;
groundPri = 40;
feedPri = 30;

% srr setup
srr.L = 50;

%ground setup
grnd_pos = -12.3; % ground distance from substrate
% grnd.points =  [-185.4,-185.4; -185.4,185.4; 185.4,185.4; 185.4,-185.4]' ;
grnd.xdim = 370.8 + 4*srr.L;
grnd.ydim = 370.8 + 2*srr.L;
% grnd.xdim = 370.8;
% grnd.ydim = 370.8;
grnd.points = [-grnd.xdim/2,-grnd.ydim/2;
     -grnd.xdim/2,grnd.ydim/2;
     grnd.xdim/2,grnd.ydim/2;
     grnd.xdim/2,-grnd.ydim/2]';
if (same_grnd_size == 1)
    grnd2 = grnd; 
elseif(same_grnd_size == 0)
    grnd2.xdim = 370.8;
    grnd2.ydim = 370.8;
    grnd2.points = [-grnd2.xdim/2,-grnd2.ydim/2;
         -grnd2.xdim/2,grnd2.ydim/2;
         grnd2.xdim/2,grnd2.ydim/2;
         grnd2.xdim/2,-grnd2.ydim/2]';
else
    error('Check the "same_grnd_size" value!!!');
end

%distance between Antennas (in units). (actually, distance between substrates)
backDist = 20 + 2 * abs(grnd_pos); %24.6 here in mm

%substrate setup
sub_freq = 8.5e8; % Frequency to calculate the substrate conductivity for
% substrate.epsR   = 3.38;
substrate.epsR = 4.2;
substrate.tan_delta = 0.025;
substrate.kappa  = substrate.tan_delta * 2*pi*sub_freq * EPS0*substrate.epsR; %conductivity 
substrate.thickness = 1.524;
substrate.cells = 4;
substrate.points = [-175,-145;175,-145;175,90;-175,90]';

%setup feeding
feed.pos = [18.85, -115.46]; %feeding position in x,y directions
feed.R = 50;     %feed resistance
feed2.pos = feed.pos;
feed2.R = feed.R;

% size of the simulation and dump box 
SimBox = [300 + grnd.xdim, 300 + grnd.ydim, 300 + backDist ];
dumpWidth = 400;
dumpLength = 400;

%% Setup FDTD Parameter & Excitation Function
f0 = 8.5e8; % center frequency (Hz)
fc = 1e8; % 20 dB corner frequency (Hz) -----> it determines the bandwidth (keep it less than f0)
FDTD = InitFDTD( 'NrTs', 1000000, 'EndCriteria', 1e-5);
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

% --- Calculate Points ---
% Calculate patch
load('points.mat');                     % load a matrix named "points" which contains the coordinates of the geometry
points = points(:,1:length(points)-1);  % remove the redundant last coordinate (same with the first). this line is not necessary. 
points2 = points;                       % store the points for patch2
points = rotate_points(points,rot,1);   % rotate the geometry. The third argument enables plot (for debugging).
hold on;                                % Usefull to plot all the geometries in one diagram. (for debugging)
allpoints = points;                     % store the points to "round" the coordinates later
pointInd(1) = size(points,2);           % necessary to recall the points after "sp_round"
% Calculate patch2
points = points2;
points(1,:) = - points2(1,:);               % mirror the patch
points = rotate_points(points,rot,1);       % rotate the geometry. The third argument enables plot.
allpoints = [allpoints, points];            % store the points to "round" the coordinates later
pointInd(2) = pointInd(1)+size(points,2);   % necessary to recall the points after "sp_round"

% Calculate substrate
points = substrate.points;                  
points = rotate_points(points,rot,1);       % rotate the geometry. The third argument enables plot.
allpoints = [allpoints, points];            % store the points to "round" the coordinates later
pointInd(3) = pointInd(2)+size(points,2);   % necessary to recall the points after "sp_round"
% Calculate substrate2
points = substrate.points;
points(1,:) = - substrate.points(1,:);      % mirror
points = rotate_points(points,rot,1);       % rotate the geometry. The third argument enables plot.
allpoints = [allpoints, points];            % store the points to "round" the coordinates later
pointInd(4) = pointInd(3)+size(points,2);   % necessary to recall the points after "sp_round"

% Calculate ground
points = grnd.points;                       
points = rotate_points(points,rot,1);       % rotate the geometry. The third argument enables plot.
allpoints = [allpoints, points];            % store the points to "round" the coordinates later
pointInd(5) = pointInd(4)+size(points,2);   % necessary to recall the points after "sp_round"
% Calculate ground2
points = grnd2.points;                       
points(1,:) = -grnd2.points(1,:);            % mirror
points = rotate_points(points,rot,1);       % rotate the geometry. The third argument enables plot.
allpoints = [allpoints, points];            % store the points to "round" the coordinates later
pointInd(6) = pointInd(5)+size(points,2);   % necessary to recall the points after "sp_round"

% Calculate feed
feed.pos = (rotate_points(feed.pos',rot,1))' ;      % adjust the feeding position to the rotated structure
% Calculate feed2
feed2.pos(1) = -feed2.pos(1);                       % mirror
feed2.pos = (rotate_points(feed2.pos',rot,1))' ;    % adjust the feeding position to the rotated structure


% import ground with SRRs calling ParamMetaSlab
in_SRR_points = [];
in_SRR_points2 = [];
if ((add_srrs == 1) || (add_srrs == 2))
    [~,~,~,~, in_SRR_points, CSX,mesh] = ParamMetaSlab('L',srr.L,'CSX',CSX,'grndelev',grnd_pos+backDist/2,'mesh',mesh,'grndxDim',grnd.xdim,'grndyDim',grnd.ydim,'inv',0);    % import SRRs with ground
    if (add_srrs == 2 && sec_antenna == 1)
        [~,~,~,~, in_SRR_points2, CSX,mesh] = ParamMetaSlab('L',srr.L,'CSX',CSX,'grndelev',-grnd_pos-backDist/2,'mesh',mesh,'grndxDim',grnd2.xdim,'grndyDim',grnd2.ydim,'inv',1);    % import SRRs with ground
    end
elseif (add_srrs == 0)
% --- Nothing special
else
    error('Check the "add_srrs" value!!!');
end

%ROUND ALL POINTS
% allpoints = round(allpoints,2); %bad practice
allpoints = sp_round(allpoints,0.5,[in_SRR_points,in_SRR_points2,[feed2.pos]' , [feed.pos]']); % "round" the coordinates of some vertices, to reduce the mesh. Don't change the coordinates of SRRs or feeds



% --- CREATE ---
% Create Patch
points = allpoints(:,1:pointInd(1));                % recall the (rounded) points
CSX = AddMetal( CSX, 'patch' );                     % create a perfect electric conductor (PEC) named "patch"
CSX = AddPolygon(CSX,'patch',patchPri,2,backDist/2+substrate.thickness,points(1:2,:));  % create a polygon of the material "patch"
% % Create patch 2
if (sec_antenna == 1)
points = allpoints(:,pointInd(1)+1:pointInd(2));    % recall the (rounded) points
CSX = AddMetal( CSX, 'patch2' );                    % create a perfect electric conductor (PEC) named "patch2"
CSX = AddPolygon(CSX,'patch2',patchPri+1,2,-(backDist/2+substrate.thickness),points(1:2,:));    % create a polygon of the material "patch2"
end

% Create Substrate
points = allpoints(:,pointInd(2)+1:pointInd(3));                                                    % recall the (rounded) points
CSX = AddMaterial( CSX, 'substrate' );                                                              % create a new material named "substrate"
CSX = SetMaterialProperty( CSX, 'substrate', 'Epsilon', substrate.epsR, 'Kappa', substrate.kappa ); % define the properties of the material "substrate"
CSX = AddLinPoly(CSX, 'substrate', substratePri, 2, backDist/2+ 0,points, substrate.thickness);     % create a polygon of the material "substrate" with thickness
% add extra cells to discretize the substrate thickness
mesh.z = [linspace(backDist/2+0,backDist/2 + substrate.thickness,substrate.cells+1) mesh.z];
% Create substrate 2
if (sec_antenna == 1)
points = allpoints(:,pointInd(3)+1:pointInd(4));                                                    % recall the (rounded) points
CSX = AddMaterial( CSX, 'substrate2' );                                                             % create a new material named "substrate2"
CSX = SetMaterialProperty( CSX, 'substrate2', 'Epsilon', substrate.epsR, 'Kappa', substrate.kappa );% define the properties of the material "substrate2"
CSX = AddLinPoly(CSX, 'substrate2', substratePri+1, 2, -(backDist/2+ 0),points, -substrate.thickness);% create a polygon of the material "substrate2" with thickness
% add extra cells to discretize the substrate thickness
mesh.z = [linspace(-(backDist/2+0),-(backDist/2 + substrate.thickness),substrate.cells+1) mesh.z];
end

if (add_srrs == 0)
% Create Ground
points = allpoints(:,pointInd(4)+1:pointInd(5));        % recall the (rounded) points
CSX = AddMetal( CSX, 'gnd' );                           % create a perfect electric conductor (PEC) named "gnd"
CSX = AddLinPoly(CSX, 'gnd', groundPri, 2, backDist/2 + grnd_pos,points, -1);   % create a polygon of the material "gnd"
end
% Create ground 2
if (sec_antenna == 1 && (add_srrs == 1 || add_srrs == 0))
points = allpoints(:,pointInd(5)+1:pointInd(6));        % recall the (rounded) points
CSX = AddMetal( CSX, 'gnd2' );                          % create a perfect electric conductor (PEC)
CSX = AddLinPoly(CSX, 'gnd2', groundPri+1, 2, -(backDist/2 + grnd_pos),points, 1);% create a polygon of the material "gnd2"
end


% --- Apply the Excitation & Resist as a Current Source ---
start = [feed.pos, backDist/2 + grnd_pos];
stop = [feed.pos, backDist/2 + substrate.thickness];
[CSX, port{1}] = AddLumpedPort(CSX, feedPri ,1 ,feed.R, start, stop, [0 0 1], true);
% port2
start = [feed2.pos, -(backDist/2 + grnd_pos)];
stop = [feed2.pos, -(backDist/2 + substrate.thickness)];
[CSX, port{2}] = AddLumpedPort(CSX, feedPri+1 ,2 ,feed2.R, start, stop, [0 0 1]);



% Finalize the Mesh
% -----------------
% detect all edges except of the patch
mesh = DetectEdges(CSX, mesh,'ExcludeProperty',['SRRpatch1', 'SRRpatch2','SRRpatch3','SRRpatch4','patch']);
% detect and set a special 2D metal edge mesh for the patch
mesh = DetectEdges(CSX, mesh,'SetProperty',['SRRpatch1', 'SRRpatch2','SRRpatch3','SRRpatch4','patch'],'2D_Metal_Edge_Res', c0/(f0+fc)/unit/50);
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
if(sec_antenna == 1)
    Sim_Path = [Sim_Path, 'Antennas'];
elseif (sec_antenna ==0)
    Sim_Path = [Sim_Path, 'SingleAntenna'];
else
    error('Check the "sec_antenna" value!!!');    
end
if(add_srrs == 2)
    Sim_Path = [Sim_Path, '_with_SRRs', '_L',num2str(srr.L),'_grndX',num2str(grnd.xdim),'_twoSides'];
elseif(add_srrs == 1)
    Sim_Path = [Sim_Path, '_with_SRRs', '_L',num2str(srr.L),'_grndX',num2str(grnd.xdim),'_oneSide'];
elseif(add_srrs == 0)
    Sim_Path = [Sim_Path, '_without_SRRs'];
else
    error('Check the "add_srrs" value!!!');
end

if(sec_antenna == 1)
    if(same_grnd_size == 1)   
        Sim_Path = [Sim_Path, '_sameSize'];    
    end
end

Sim_CSX = [Sim_Path(5:end), '.xml'];
% Sim_CSX = 'Ant_with_SRRs_simulation.xml';

% create an empty working directory
[status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
[status, message, messageid] = mkdir( Sim_Path ); % create empty simulation folder

% write openEMS compatible xml-file
WriteOpenEMS( [Sim_Path '/' Sim_CSX], FDTD, CSX );

% show the structure
CSXGeomPlot( [Sim_Path '/' Sim_CSX] );

% run openEMS
if (on_server == 1)
    Settings.SSH.Putty.Path = 'C:/"Program Files"/PuTTY';
    Settings.SSH.Putty.Key = 'C:/.ssh/openEmsPuttyKey_priv.ppk';
%     Settings.SSH.host = %%%%%% ;
    Settings.SSH.bin='/home/lamda/opt/openEMS/bin/openEMS.sh';
    RunOpenEMS( Sim_Path, Sim_CSX,'',Settings);
else
    RunOpenEMS( Sim_Path, Sim_CSX);
end

%% Postprocessing & Plots
freq = linspace((f0-fc), (f0+fc), 501 );
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
%f_res_ind = find(s11==min(s11));
f_res_ind = 295;
f_res = freq(f_res_ind);



% calculate the far field at phi=0 degrees and at phi=90 degrees
disp( 'calculating far field at phi=[0 90] deg...' );

nf2ff = CalcNF2FF(nf2ff, Sim_Path, f_res, (-180:2:180)*pi/180, [0 90]*pi/180);

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

% % Show 3D pattern
% disp( 'calculating 3D far field pattern and dumping to vtk (use Paraview to visualize)...' );
% thetaRange = (0:2:180);
% phiRange = (0:2:360) - 180;
% nf2ff = CalcNF2FF(nf2ff, Sim_Path, f_res, thetaRange*pi/180, phiRange*pi/180,'Verbose',1,'Outfile','3D_Pattern.h5');
% 
% figure
% plotFF3D(nf2ff,'logscale',-20);
% 
% 
% E_far_normalized = nf2ff.E_norm{1} / max(nf2ff.E_norm{1}(:)) * nf2ff.Dmax;
% DumpFF2VTK([Sim_Path '/3D_Pattern.vtk'],E_far_normalized,thetaRange,phiRange,'scale',1e-3);


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
tit = Sim_Path(5:end);
tit(tit == '_') = ' ';
title (tit);


% store workspace variables
save([Sim_Path, '/', Sim_Path] );