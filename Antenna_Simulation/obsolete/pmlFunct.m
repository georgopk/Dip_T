%           D  = depth in the pml in meter
% 			dl = mesh delta inside the pml in meter %---------> skin depth???  tan(delta) ????
% 			W  = width (length) of the pml in meter
% 			N  = number of cells for the pml
% 			Z  = wave impedance at the current depth and position



 Z = 50;



dl = 0.01;
W = 0.01;
N = 8;
D=linspace (0,W,100);

Y = -log(1e-6)*log(2.5)/(2*dl*power(2.5,W/dl)-1) * power(2.5, D/dl) / Z;

plot(D,Y)