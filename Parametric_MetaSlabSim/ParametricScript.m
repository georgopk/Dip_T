clear
% declare the parameters
dplen = 140:10:180;
disp (dplen);   % for convenience

for i = 1: length(dplen)
    
    % Simulation    
    disp(['Dipole length = ',num2str(dplen(i))]);   % for convenience
    [tfreq, tport, ts11, ts21 ] = ParamMetaSlab('dipLength',dplen(i));  % Run the simulation and store the results
    outVal{i,:} = {tfreq, tport, ts11, ts21};       % aggregate all the results 
    
    
    %%% plot during the process
    % plot reflection coefficient S11
    freq = cell2mat(outVal{i}(1));
    s11 = cell2mat(outVal{i}(3));
    figure(1)
    hold on
    plot( freq/1e6, 20*log10(abs(s11)), 'k-', 'Linewidth', 2 );
    grid on
    title( 'reflection coefficient S_{11}' );
    xlabel( 'frequency f / MHz' );
    ylabel( 'reflection coefficient |S_{11}|' );

    % plot coefficient S21
    freq = cell2mat(outVal{i}(1));
    s21 = cell2mat(outVal{i}(4));
    figure(2)
    hold on
    plot( freq/1e6, 20*log10(abs(s21)), 'r', 'Linewidth', 2 );
    grid on
    title( 'S_{21}' );
    xlabel( 'frequency f / MHz' );
    ylabel( '|S_{21}|' );
end


% % %% OR plot after the process
% % 
% % for i = 1: length(dplen)
% %     % plot reflection coefficient S11
% %     freq = cell2mat(outVal{i}(1));
% %     s11 = cell2mat(outVal{i}(3));
% %     figure(1)
% %     hold on
% %     plot( freq/1e6, 20*log10(abs(s11)), 'k-', 'Linewidth', 2 );
% %     grid on
% %     title( 'reflection coefficient S_{11}' );
% %     xlabel( 'frequency f / MHz' );
% %     ylabel( 'reflection coefficient |S_{11}|' );
% % 
% %     % plot coefficient S21
% %     freq = cell2mat(outVal{i}(1));
% %     s21 = cell2mat(outVal{i}(4));
% %     figure(2)
% %     hold on
% %     plot( freq/1e6, 20*log10(abs(s21)), 'r', 'Linewidth', 2 );
% %     grid on
% %     title( 'S_{21}' );
% %     xlabel( 'frequency f / MHz' );
% %     ylabel( '|S_{21}|' );
% % end
