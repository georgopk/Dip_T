%% Simple polygon patch antenna (based on...) % Simple Patch Antenna Tutorial
%   Simulates a patch antenna with the shape of a regular polygon.
% --------------------------------
% Creates:
%  -a patch antenna with the shape of a regular polygon
%  -a Rectangle substrate
%  -a Rectangle ground plane
%
% Calculates:
% - the feed point impedance
% - the frequency of the minimum impedance
% (for the frequency range [f0-fc , f0+fc] ) 
% (a lower frequency boundary is set)
%     
%     
% Calculates:
% - the radiation pattern at phi = 0 and phi = 90
% - the 3d representation of the radiation pattern
% (for the frequency of minimum impedance )
% 
% 
% % -------Tutorial info----------
% % Describtion at:
% % <http://openems.de/index.php/Tutorial:_Simple_Patch_Antenna>
% %
% % Tested with
% %  - Matlab 2013a / Octave 4.0
% %  - openEMS v0.0.35
% %
% % (C) 2010-2017 Thorsten Liebig <thorsten.liebig@uni-due.de>
% 



%%

close all
clear
clc

%% Setup the Simulation
physical_constants;
unit = 1e-3; % all length in mm

% patch width in x-direction
patch.width  = 32; % resonant length
% patch length in y-direction
patch.length = 40;

%substrate setup
substrate.epsR   = 3.38;
substrate.kappa  = 1e-3 * 2*pi*2.45e9 * EPS0*substrate.epsR;
substrate.width  = 340;
substrate.length = 170;
substrate.thickness = 1.524;
substrate.cells = 4;

%setup feeding
feed.pos = [58, 49]; %feeding position in x-direction
feed.R = 50;     %feed resistance

% size of the simulation box
SimBox = [600 300 100];

%% Setup FDTD Parameter & Excitation Function
f0 = 9e8; % center frequency
fc = 3e8; % 20 dB corner frequency -----> it determines the bandwidth (keep it less than f0)
FDTD = InitFDTD( 'NrTs', 30000, 'EndCriteria', 1e-5);
FDTD = SetGaussExcite( FDTD, f0, fc );
BC = {'MUR' 'MUR' 'MUR' 'MUR' 'MUR' 'MUR'}; % boundary conditions
FDTD = SetBoundaryCond( FDTD, BC );

%% initialize complex Geometry
%initialize poloygon
poly.cen =[0,0,0]; % define Polygon Center 
poly.N = 16;  % define the shape - number of vertices (6-gon, 8-gon, N-gon)
poly.ang_centr = 2*pi / poly.N; % Central Angle of the polygon (radians)
poly.ang_inter = pi-(2*pi/poly.N); % Interior Angle of the polygon
poly.side_length = 30;  % define the side length in mm (actually in the chosen unit - see "unit" variable)
poly.r = poly.side_length/(2*sin(pi/poly.N)); % radius of the polygon (see inscribed polygons)
poly.rotation = 0.7; %rad
for i=1:poly.N
    temp = poly.rotation + (i-1) * poly.ang_centr;
    poly.vert(:,i) = [poly.cen(1), poly.cen(2)] + poly.r * [cos(temp), sin(temp)];
end

poly2.vert = [(poly.vert(1,:) +148) ; poly.vert(2,:)];

%% Setup CSXCAD Geometry & Mesh
CSX = InitCSX();

%initialize the mesh with the "air-box" dimensions
mesh.x = [-SimBox(1)/2 SimBox(1)/2];
mesh.y = [-SimBox(2)/2 SimBox(2)/2];
mesh.z = [-SimBox(3)/2 SimBox(3)/2];



% Create Patch
CSX = AddMetal( CSX, 'patch' ); % create a perfect electric conductor (PEC)
CSX = AddPolygon(CSX,'patch',10,2,0,poly.vert); % add a box-primitive to the metal property 'patch'

CSX = AddPolygon(CSX,'patch',11,2,0,poly2.vert); % add a 2nd box-primitive to the metal property 'patch'

% Create Substrate
CSX = AddMaterial( CSX, 'substrate' );
CSX = SetMaterialProperty( CSX, 'substrate', 'Epsilon', substrate.epsR, 'Kappa', substrate.kappa );
start = [-substrate.width/2+74 -substrate.length/2 -substrate.thickness];
stop  = [ substrate.width/2+74  substrate.length/2 0];
CSX = AddBox( CSX, 'substrate', 0, start, stop );

% add extra cells to discretize the substrate thickness
mesh.z = [linspace(0,substrate.thickness,substrate.cells+1) mesh.z];

% Create Ground same size as substrate
CSX = AddMetal( CSX, 'gnd' ); % create a perfect electric conductor (PEC)
start(3)=-substrate.thickness;
stop(3) =-substrate.thickness;
CSX = AddBox(CSX,'gnd',10,start,stop);

% Apply the Excitation & Resist as a Current Source
start = [feed.pos(1) feed.pos(2) 0];
stop  = [feed.pos(1) feed.pos(2) substrate.thickness];
[CSX port{1}] = AddLumpedPort(CSX, 5 ,1 ,feed.R, start, stop, [0 0 1], true);


% Finalize the Mesh
% -----------------
% detect all edges except of the patch
mesh = DetectEdges(CSX, mesh,'ExcludeProperty','patch');
% detect and set a special 2D metal edge mesh for the patch
mesh = DetectEdges(CSX, mesh,'SetProperty','patch','2D_Metal_Edge_Res', c0/(f0+fc)/unit/50);
% generate a smooth mesh with max. cell size: lambda_min / 20 
mesh = SmoothMesh(mesh, c0/(f0+fc)/unit/20);
CSX = DefineRectGrid(CSX, unit, mesh);

CSX = AddDump(CSX,'Hf', 'DumpType', 11, 'Frequency',[9e8]);
CSX = AddBox(CSX,'Hf',10,[-(74+substrate.width/2+10) -(substrate.length/2+10) -10*substrate.thickness],[74+substrate.width/2+10 +substrate.length/2+10 10*substrate.thickness]); %assign box

% add a nf2ff calc box; size is 3 cells away from MUR boundary condition
start = [mesh.x(4)     mesh.y(4)     mesh.z(4)];
stop  = [mesh.x(end-3) mesh.y(end-3) mesh.z(end-3)];
[CSX nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', start, stop);

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
% s21 = port{2}.uf.ref./ port{1}.uf.inc;
% %% plot s-parameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % my plots
% figure
% plot(freq*1e-6,20*log10(abs(s11)),'k-','Linewidth',2);
% xlim([freq(1) freq(end)]*1e-6);
% grid on;
% hold on;
% plot(freq*1e-6,20*log10(abs(s21)),'r--','Linewidth',2);
% l = legend('S_{11}','S_{21}','Location','Best');
% set(l,'FontSize',12);
% ylabel('S-Parameter (dB)','FontSize',12);
% xlabel('frequency (MHz) \rightarrow','FontSize',12);