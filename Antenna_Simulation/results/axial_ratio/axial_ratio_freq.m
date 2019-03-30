function axial_ratio_freq(freq, nf2ff)

k=1;

axr_f_mag_phi0 = zeros(1, 501) ;
axr_f_phase_phi0 = zeros(1, 501) ;
axr_f_mag_phi90 = zeros(1, 501) ;
axr_f_phase_phi90 = zeros(1, 501) ;

A = cell2mat(nf2ff.E_theta) ;
B = cell2mat(nf2ff.E_phi) ;

for i=1:2:(2*length(freq))
    axr_f_mag_phi0(1, k) = 20*log10(abs(A(92, i)/B(92, i))) ;
    axr_f_phase_phi0(1, k) = angle(A(92, i)/B(92, i))*180/pi ; 
    axr_f_mag_phi90(1, k) = 20*log10(abs(A(92, 1+i)/B(92, 1+i))) ;
    axr_f_phase_phi90(1, k) = angle(A(92, 1+i)/B(92, 1+i))*180/pi ;
    k = k + 1;
end
 
figure
plot(freq, axr_f_mag_phi0);
grid on
title('axial ratio(freq) magnitude(dB) phi0');

figure
plot(freq, axr_f_phase_phi0);
grid on
title('axial ratio(freq) phase phi0');

figure
plot(freq, axr_f_mag_phi90);
grid on
title('axial ratio(freq) magnitude(dB) phi90');

figure
plot(freq, axr_f_phase_phi90);
grid on
title('axial ratio(freq) phase phi90');
