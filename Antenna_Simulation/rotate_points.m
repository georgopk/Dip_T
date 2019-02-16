function pts2 = rotate_points(points, ang, plot_on)
% This function rotates given points on XY plane counterclockwise and plots
% the result.
% 
% -------Arguments-------
% Points:   a matrix containing the coordinates of the points. 
%           The first row is suposed to be X and the second row is suposed
%           to be Y. Any other row will be ignored.
% ang:      angle of rotation in rad
%-------Optional-------
% plot      Using a value in this argument, it plots the result. Useful for
%           debugging.
%
%

rot_mat = [cos(ang), -sin(ang); sin(ang), cos(ang)];    % evaluate rotation matrix
pts2 = rot_mat * (points(1:2,:));       % apply the rotation 

% plot
if exist('plot_on','var')
    plot(pts2(1,:), pts2(2,:))
    xlim([-250, 250])
    ylim([-250, 250])
    grid on;
    grid minor;
end

end