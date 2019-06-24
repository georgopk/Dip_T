L=[45:53];


% Load saved figures
for i = 1:size(L,2)
    figure(i)=openfig(['L-',num2str(L(i)), '_S11.fig' ]);
    % k=hgload('MySecondFigure.fig');
end

% Prepare subplots
figure
for i = 2:size(L,2)
h(i)=findobj(i,'type','line');
% h(2)=subplot(1,2,2);
end

% Paste figures on the subplots
for i = 2:size(L,2)
copyobj(h(i),findobj(1,'type','axes'));
% copyobj(allchild(get(k,'CurrentAxes')),h(2));
end

% % Add legends
% l(1)=legend(h(1),'LegendForFirstFigure')
% l(2)=legend(h(2),'LegendForSecondFigure')