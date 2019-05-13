clear
% declare the parameters
dplen = 160;
disp (dplen);   % for convenience
L=[50,70];

tic
for j = 1:length(L)
    disp('*')
    disp('*')
    disp('*')
    disp('*')
    disp('*')
    disp('*')
    disp('*')
    disp('*')
    disp('*')
    disp('*')
    disp('*')
    disp('*')
    disp('*')
    disp([ 'Simulating @L=' , num2str(L(j)) ])
    for i = 1: length(dplen)

        % Simulation    
        disp(['Dipole length = ',num2str(dplen(i))]);   % for convenience
        [tfreq, tport, ts11, ts21 ] = ParamMetaSlab('dipLength',dplen(i),'L',L(j));  % Run the simulation and store the results
        outVal{i,:} = {tfreq, tport, ts11, ts21};       % aggregate all the results 

        %plot during the process
        myplot(i,j,outVal,L(j))
    end
end
toc
disp(toc/3600)
% % % OR plot after the process
% % for i = 1: length(dplen)
% %    myplot(i,outVal,L)
% % end


% functions
function myplot(n,m,outVal,L)
    freq = cell2mat(outVal{n}(1));
    s11 = abs(cell2mat(outVal{n}(3)));
    
    tmp1 = find(s11 == min(abs(s11)));
    if length(tmp1)>1
        warning('The resonance value apears in more than one position!!!')
    end
    ind = freq-freq(tmp1) > -350e6 & freq-freq(tmp1) < 350e6;
      
    % plot reflection coefficient S11
    figure(2*m-1)
    hold on
    plot( freq(ind)/1e6, 20*log10(abs(s11(ind))), 'k-', 'Linewidth', 2 );
    grid on
    title(['reflection coefficient S_{11} @ L=', num2str(L)]);
    xlabel( 'frequency f / MHz' );
    ylabel( 'reflection coefficient |S_{11}|' );

    % plot coefficient S21
    freq = cell2mat(outVal{n}(1));
    s21 = cell2mat(outVal{n}(4));
    figure(2*m)
    hold on
    plot( freq(ind)/1e6, 20*log10(abs(s21(ind))), 'r', 'Linewidth', 2 );
    grid on
    title( ['S_{21}    @ L=', num2str(L)] );
    xlabel( 'frequency f / MHz' );
    ylabel( '|S_{21}|' );
end