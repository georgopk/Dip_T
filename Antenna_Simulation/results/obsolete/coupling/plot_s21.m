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