% close all
clear
% clc

%% setup the simulation
physical_constants;
unit = 1e-3; % all length in mm

%% Materials

% copper
copper.conductivity = 56e6;

% substrate
substrate.epsR   = 3.38;
substrate.kappa  = 1e-3 * 2*pi*2.45e9 * EPS0*substrate.epsR;

%% Dimensions
% val.g = 0.3; % gap width of slow in rings
val.g = 1;
val.h = 0.25; % substrate height
% val.l1 = 2.2; % side length of outer ring
val.l1 = 30;
val.ls = 2.5; % side length of substrate
val.lw = 0.14; % width of wire underneath substrate
% val.s = 0.15; % separation between rings
val.s = 2;
val.t = 0.017; % thickness of rings (above substrate; PEC with 0 ok for now).
% val.w = 0.2; % width of wire underneath substrate
val.w = 1;
val.l2 = val.l1-val.w-val.s; % side length of innter ring


% [points1,points2] = srr_points(srr.L,srr.w,srr.g,srr.s,plot_on);

ring.a.gap = val.g;
ring.a.t = val.t;
ring.a.w = val.w;
ring.a.length = val.l2;

ring.b.gap = val.g;
ring.b.t = val.t;
ring.b.w = val.w;
ring.b.length = ring.a.length + 2*(val.s) + 2*(ring.b.w); %val.l1;

substrate.h = val.h;
substrate.length = 1.2*val.l1;
% substrate.length = val.ls;

wire.w = val.lw;
wire.t = val.t;
wire.length = substrate.length;

%% Structures
% inner ring
% x axis stops
ring.a.x3p = 0.5*ring.a.length;
ring.a.x2p = (ring.a.x3p - ring.a.w);
ring.a.x1p = 0.5*ring.a.gap;
ring.a.x1n = -ring.a.x1p;
ring.a.x2n = -ring.a.x2p;
ring.a.x3n = -ring.a.x3p;
% y axis stops
ring.a.y2p = 0.5*ring.a.length;
ring.a.y1p = (ring.a.y2p - ring.a.w);
ring.a.y1n = -ring.a.y1p;
ring.a.y2n = -ring.a.y2p;

ring.a.p(1,1) = (ring.a.x3p);   ring.a.p(2,1) = (ring.a.y2p);
ring.a.p(1,2) = (ring.a.x3n);   ring.a.p(2,2) = (ring.a.y2p);
ring.a.p(1,3) = (ring.a.x3n);   ring.a.p(2,3) = (ring.a.y2n);
ring.a.p(1,4) = (ring.a.x1n);   ring.a.p(2,4) = (ring.a.y2n);
ring.a.p(1,5) = (ring.a.x1n);   ring.a.p(2,5) = (ring.a.y1n);
ring.a.p(1,6) = (ring.a.x2n);   ring.a.p(2,6) = (ring.a.y1n);
ring.a.p(1,7) = (ring.a.x2n);   ring.a.p(2,7) = (ring.a.y1p);
ring.a.p(1,8) = (ring.a.x2p);   ring.a.p(2,8) = (ring.a.y1p);
ring.a.p(1,9) = (ring.a.x2p);   ring.a.p(2,9) = (ring.a.y1n);
ring.a.p(1,10) = (ring.a.x1p);  ring.a.p(2,10) = (ring.a.y1n);
ring.a.p(1,11) = (ring.a.x1p);  ring.a.p(2,11) = (ring.a.y2n);
ring.a.p(1,12) = (ring.a.x3p);  ring.a.p(2,12) = (ring.a.y2n);
ring.a.p(1,13) = (ring.a.x3p);  ring.a.p(2,13) = (ring.a.y2p);

% outer ring
% x axis stops
ring.b.x3p = 0.5*ring.b.length;
ring.b.x2p = (ring.b.x3p - ring.b.w);
ring.b.x1p = 0.5*ring.b.gap;
ring.b.x1n = -ring.b.x1p;
ring.b.x2n = -ring.b.x2p;
ring.b.x3n = -ring.b.x3p;
% y axis stops
ring.b.y2p = 0.5*ring.b.length;
ring.b.y1p = (ring.b.y2p - ring.b.w);
ring.b.y1n = -ring.b.y1p;
ring.b.y2n = -ring.b.y2p;

ring.b.p(1,1) = (ring.b.x3p);   ring.b.p(2,1) = (ring.b.y2p);
ring.b.p(1,2) = (ring.b.x1p);   ring.b.p(2,2) = (ring.b.y2p);
ring.b.p(1,3) = (ring.b.x1p);   ring.b.p(2,3) = (ring.b.y1p);
ring.b.p(1,4) = (ring.b.x2p);   ring.b.p(2,4) = (ring.b.y1p);
ring.b.p(1,5) = (ring.b.x2p);   ring.b.p(2,5) = (ring.b.y1n);
ring.b.p(1,6) = (ring.b.x2n);   ring.b.p(2,6) = (ring.b.y1n);
ring.b.p(1,7) = (ring.b.x2n);   ring.b.p(2,7) = (ring.b.y1p);
ring.b.p(1,8) = (ring.b.x1n);   ring.b.p(2,8) = (ring.b.y1p);
ring.b.p(1,9) = (ring.b.x1n);   ring.b.p(2,9) = (ring.b.y2p);
ring.b.p(1,10) = (ring.b.x3n);   ring.b.p(2,10) = (ring.b.y2p);
ring.b.p(1,11) = (ring.b.x3n);   ring.b.p(2,11) = (ring.b.y2n);
ring.b.p(1,12) = (ring.b.x3p);   ring.b.p(2,12) = (ring.b.y2n);
ring.b.p(1,13) = (ring.b.x3p);   ring.b.p(2,13) = (ring.b.y2p);

% % (L,w,g,s,plot_on)
% [ring.b.p, ring.a.p] = srr_points(val.l1,val.w,val.g,val.s,1);

ring.thickness = val.t;

% substrate
substrate.width  = substrate.length;
substrate.thickness = substrate.h;
substrate.cells = 4;

% size of the simulation box
SimBox = [1,1,1] + substrate.length + 20;

%% setup FDTD parameter & excitation function
f0 = 8.5e8; % center frequency
fc = 4e8; % 20 dB corner frequency
FDTD = InitFDTD('NrTS', 400000 );
FDTD = SetGaussExcite( FDTD, f0, fc );
%BC = {'MUR' 'MUR' 'MUR' 'MUR' 'MUR' 'MUR'}; % boundary conditions
BC = {'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8'}; % boundary conditions
FDTD = SetBoundaryCond( FDTD, BC );

%% setup CSXCAD geometry & mesh
max_res = c0 / (f0+fc) / unit / 20; % cell size: lambda/20
CSX = InitCSX();

%create fixed lines for the simulation box, substrate and port
mesh.x = [-SimBox(1)/2 SimBox(1)/2 -substrate.width/2 substrate.width/2 ring.b.x3n ring.b.x3p ring.b.x2n ring.b.x2p ring.a.x3n ring.a.x3p];
mesh.x = SmoothMeshLines( mesh.x, max_res, 1.4); % create a smooth mesh between specified fixed mesh lines
mesh.y = [-SimBox(2)/2 SimBox(2)/2 -substrate.length/2 substrate.length/2 ring.b.y2n ring.b.y2p ring.b.y1n ring.b.y1p];
mesh.y = SmoothMeshLines( mesh.y, max_res, 1.4 );
%create fixed lines for the simulation box and given number of lines inside the substrate
mesh.z = [-SimBox(3)/2 linspace(0,substrate.thickness,substrate.cells) SimBox(3)/2 ];
mesh.z = SmoothMeshLines( mesh.z, max_res, 1.4 );
CSX = DefineRectGrid( CSX, unit, mesh );

%% create patch
%CSX = AddConductingSheet( CSX, 'copper', copper.conductivity, 0.0004);
CSX = AddMetal( CSX, 'copper');
CSX = AddLinPoly( CSX, 'copper', 0, 2, substrate.thickness, ring.a.p, 0, 'CoordSystem',0);
CSX = AddLinPoly( CSX, 'copper', 0, 2, substrate.thickness, ring.b.p, 0, 'CoordSystem',0);

%% create substrate
CSX = AddMaterial( CSX, 'substrate' );
CSX = SetMaterialProperty( CSX, 'substrate', 'Epsilon', substrate.epsR, 'Kappa', substrate.kappa );
start = [-substrate.width/2 -substrate.length/2 0];
stop  = [ substrate.width/2  substrate.length/2 substrate.thickness];
CSX = AddBox( CSX, 'substrate', 0, start, stop );

%% create ground (same size as substrate)
%CSX = AddMetal( CSX, 'gnd' ); % create a perfect electric conductor (PEC)
%start = [-substrate.width/2 -substrate.length/2 0];
%stop  = [ substrate.width/2  substrate.length/2 0];
%CSX = AddBox(CSX,'gnd',10,start,stop);

%% create wire
start = [-wire.w/2 -substrate.length/2 -wire.t];
stop  = [ wire.w/2  substrate.length/2 0];
CSX = AddBox(CSX,'copper',10,start,stop);

%% apply the excitation
%CSX = AddExcitation(CSX,'excitation',0,[0 1 0 ]);
%start = [-substrate.length/2 -substrate.length/2 -wire.t];
%stop = [-substrate.length/2 substrate.length/2 2*substrate.thickness + wire.t];
%CSX = AddBox(CSX,'excitation',0,start,stop);

%CSX = AddExcitation(CSX,'excitation',0,[1 0 0 ]);
%start = [substrate.length/2 -substrate.length/2 -wire.t];
%stop = [substrate.length/2 substrate.length/2 2*substrate.thickness + wire.t];
%CSX = AddBox(CSX,'excitation',0,start,stop);

%CSX = AddDump(CSX,'Et','DumpMode',0);
%start = [-substrate.length/2 0 -wire.t];
%stop = [substrate.length/2 0 2*substrate.thickness + wire.t];
%CSX = AddBox(CSX,'Et',0,start,stop);

%% Add ports
%waveguide TE-mode definition
TE_mode = 'TE10';

%start = [(-substrate.length/2)+0.0 -substrate.length/2 -wire.t];
%stop = [(-substrate.length/2)+0.1 substrate.length/2 2*substrate.thickness + wire.t];
%[CSX, port{1}] = AddRectWaveGuidePort(CSX,0,1,start,stop, 'x',(stop(2)-start(2))*unit,(stop(3)-start(3))*unit,TE_mode, 1);
start = [(-substrate.width/2)+0.0 -substrate.length/2 -wire.t];
stop = [(-substrate.width/2)+0.1 substrate.length/2 2*substrate.thickness + wire.t];
[CSX, port{1}] = AddRectWaveGuidePort(CSX,0,1,start,stop, 'x',(stop(2)-start(2))*unit,(stop(3)-start(3))*unit,TE_mode, 1);
start = [-start(1), start(2), start(3)];
stop = [-stop(1), stop(2), stop(3)];
[CSX, port{2}] = AddRectWaveGuidePort(CSX,0,2,start,stop, 'x',(stop(2)-start(2))*unit,(stop(3)-start(3))*unit,TE_mode, 0);

%% add a nf2ff calc box
%SimBox = SimBox - max_res * 4; %reduced SimBox size for nf2ff box
%[CSX nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', -SimBox/2, SimBox/2);

%% prepare simulation folder
Sim_Path = 'tmp_SRR2';
Sim_CSX = 'patch_ant.xml';

[status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
[status, message, messageid] = mkdir( Sim_Path ); % create empty simulation folder

%% write openEMS compatible xml-file
WriteOpenEMS( [Sim_Path '/' Sim_CSX], FDTD, CSX );

%% show the structure
CSXGeomPlot( [Sim_Path '/' Sim_CSX] );

%% run openEMS
%RunOpenEMS( Sim_Path, Sim_CSX);
RunOpenEMS( Sim_Path, Sim_CSX ,'-v');

%% postprocessing & do the plots
freq = linspace( max([1e9,f0-fc]), f0+fc, 501 );
port = calcPort(port, Sim_Path, freq);

Zin = port{1}.uf.tot ./ port{1}.if.tot;
s11 = port{1}.uf.ref ./ port{1}.uf.inc;
s21 = port{2}.uf.ref ./ port{1}.uf.inc;
P_in = 0.5 * port{1}.uf.inc .* conj( port{1}.if.inc ); % antenna feed power

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


% plot reflection coefficient S21
figure
plot( freq/1e6, 20*log10(abs(s21)), 'r--', 'Linewidth', 2 );
grid on
title( 'S_{21}' );
xlabel( 'frequency f / MHz' );
ylabel( '|S_{21}|' );

drawnow