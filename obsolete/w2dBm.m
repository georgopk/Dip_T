function [ dBm ] = w2dBm( watt )
% converts Watts to dBm

dBm = 10*log10(watt*1000);


end

