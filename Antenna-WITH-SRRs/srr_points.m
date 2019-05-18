function [points_out, points_in] = srr_points(L,w,g,s,plot_on)
% Calculate points for an SRR (with 2 strips)
%
%
%    ________________________                         
%   |                        |
%   |    ________________    |                   
%   |   |                |   |
%   |   | g  ________    |   |     
%   |   |<->|  ____  |   |___|      ___     
%   |   |   |_|    | |                 \ 
%   |   | S{ _     | |                  }-> S
%   |   |   | |____| |  ->___<-w    ___/ 
%   |   |   |________|   |   |     
%   |   |                |   |
%   |   |________________|   |                    
%   |                        |
%   |________________________|                         
%    <--------- L ---------->      
% 
% 
% L: length
% w: width
% g: gap
% s: slot length
% plot_on: flag to plot
% 
% 



% % slot direction
% s_dir = 0; % 0->x, 1->y

%out
points(1,:)= [L/2, s/2];
points(2,:)=[L/2, L/2];
points(3,:)=[-L/2, L/2];
points(4,:)=[-L/2, -L/2];
points(5,:)=[L/2, -L/2];
points(6,:)=[L/2, -s/2];

points(7,:)=[L/2-w,-s/2];
points(8,:)=[L/2-w,-(L/2-w)];
points(9,:)=[-(L/2-w),-(L/2-w)];
points(10,:)=[-(L/2-w),L/2-w];
points(11,:)=[L/2-w,L/2-w];
points(12,:)=[L/2-w,s/2];

points_out = points;
if ( plot_on == 1)
    plot(points(:,1),points(:,2))
    xlim([-(L+2*w), L+2*w]/2);
    ylim([-(L+2*w), L+2*w]/2);
    hold on;
end

%In
L = L-2*w-2*g;
points(1,:)= [L/2, s/2];
points(2,:)=[L/2, L/2];
points(3,:)=[-L/2, L/2];
points(4,:)=[-L/2, -L/2];
points(5,:)=[L/2, -L/2];
points(6,:)=[L/2, -s/2];

points(7,:)=[L/2-w,-s/2];
points(8,:)=[L/2-w,-(L/2-w)];
points(9,:)=[-(L/2-w),-(L/2-w)];
points(10,:)=[-(L/2-w),L/2-w];
points(11,:)=[L/2-w,L/2-w];
points(12,:)=[L/2-w,s/2];
ang = pi;
rot_mat = [cos(ang), -sin(ang); sin(ang), cos(ang)];    % evaluate rotation matrix
points = (rot_mat * (points'))';                        % apply the rotation 
points = round(points,2);                               % necessary, because the rotation causes a small calculation error
points_in = points;
if ( plot_on == 1)
    plot(points(:,1),points(:,2))
    hold off;
end


end

