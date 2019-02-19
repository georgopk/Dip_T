%% Warning!!! This script starts clearing your workspace
% Evaluates the range for the forward link and the range for the reverse
% link


clear

%general
f=915 * 10^6; % Hz              % frequency
c=299792458; % m/s              % speed of light
wlength=c/f; % m                % wavelength 

%% forward link
%transmitter (Reader)
Ptx_Reader= 1; % Watt                   % transmission power
Ptx_dBm= w2dBm(Ptx_Reader); % dBm       % transmission power (dBm)
GReader_dBi= 6; % dBi                   % Antena Gain (dBi)
GReader=dB2factor(GReader_dBi); %noUnits % Antena Gain (factor)
%----------------------------
%amplifier
amp_dist = 5; % m                          % repeater from reader (m)
GampRx_dBi = 5; % dBil                      % Receiving antena Gain (dB) linear polarization
GampRx= dB2factor(GampRx_dBi); %noUnits     % Receiving antena Gain
GampTx_dBi = 5; % dBil                      % Transmission antena Gain (dB) linear polarization
GampTx= dB2factor(GampTx_dBi); %noUnits     % Transmission antena Gain

Ampl_dB = 30; % dB                          % Amplification (gain) (dB)
Ampl = dB2factor(Ampl_dB); %noUnits         % Amplification (gain) (factor)
Amp_feedback=0; %                           % Amplif output to amplif input

PampRx = Ptx_Reader*GReader*GampRx* (wlength/(4*pi*amp_dist))^2;%watt %Power received from repeater (Watt)
PampRx_dBm = w2dBm(PampRx); % dBm           % Power received from repeater (dBm)

PampTx = PampRx*Ampl; % Watt                % Repeater transmission power (watt)
PampTx_dBm = w2dBm(PampTx); %dBm            % Repeater transmission power (dBm)


%----------------------------
%receiver (Tag)
% Ae_tag= wlength^2 / (4*pi); % m^2           % effective aperture
Prx_min = 10^(-4); % Watt                   % minimum receive power required
Prx_min_dBm= w2dBm(Prx_min); %dBm           % minimum receive power required (dBm)



R_frwrd=amp_dist + sqrt((PampTx*GampTx /Prx_min)* (wlength/(4*pi))^2 );


%% reverse link

%Transmitter (Tag)
Tag_cons = 0;% Watt                                 % Tag consumption
conv_ef = 0.3;% percetage                           % conversion efficiency. From CW to DC (percentage)
Ptx_tag = Prx_min * conv_ef - Tag_cons;% Watt       % Backscatter power


%Receiver (Reader)
Prx_back = Ptx_tag * GReader * (wlength / (4*pi* R_frwrd))^2;       % Power received from reader
Prx_back_dBm = w2dBm(Prx_back);                                     % Power received from reader (dBm)
Prx_Reader_min_dBm = -75; % dBm                                     % Minimum power received from reader (dBm)
Prx_Reader_min= dB2factor(Prx_Reader_min_dBm) /1000; % Watt         % Minimum power received from reader
reverse_range = wlength/(4*pi) * sqrt(Ptx_tag * GReader / Prx_Reader_min);

%% Monitoring

message = ' Transmission Power (reader): %f dBm\n Receiving Power (Tag): %f dBm\n Forward Range: %f m\n Power received back: %f dBm\n Reverse Range: %f m';
sprintf(message,Ptx_dBm,Prx_min_dBm,R_frwrd,Prx_back_dBm,reverse_range)



