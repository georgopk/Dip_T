function points = sp_round(points,eee,feed)
% This function finds coordinates (in "points" matrix) with difference less 
% than "eee" and replaces them with a unique value in order to be 
% identical.
% 
% 
% points -> A matrix with the coordinates of the points. If the matrix has
%           more than 2 rows, nothing will change in them.
% eee    -> Threshold. If the difference between two values in a row is
%           less than "eee", they will be replaced in order to be
%           identical.
% feed   -> A vector 1xN (or Nx1) containing the coordinates of a point 
%           which has to be unchanged, where N is the number of rows of the 
%           matrix "points" (Useful for feeding points) Other points with 
%           coordinates similar to "feed" will change in order to fit    
%           This argument is optional.
%
% 
% It may not work properly enough for a row in which the values increase 
% in small steps.
% --- Worst case example ---
% >> pts = [4.001,5,5.01,6,4;2,2.9,3.4,3.7,3.9]
% pts =
% 
%     4.0010    5.0000    5.0100    6.0000    4.0000
%     2.0000    2.9000    3.4000    3.7000    3.9000
% >> pts2 = sp_round(pts,1) 
% pts2 =
% 
%     4.5027    4.5027    4.5027    6.0000    4.5027
%     3.1800    3.1800    3.1800    3.1800    3.1800
% --> each element of the second row has been replaced with the same value
%     (3.18). Use smaller "eee" to avoid it.


for i=1:length (points) % check every cell of a row............
%   variable "temp" is used to store the indices of the values with
%   difference less than "eee"
    temp(1:size(points,1),1:size(points,2)) = false; 
    for j=1:length (points) %............with any other cell in the row !!!
        for k = 1 : size (points,1) % for each row in "points"

            if exist('feed','var')  % if 'feed' argument is set. Check if the value is similar to any of the corresponding "feed-coordinates"
                for m=1:size(feed,2)
                    if (abs(points(k,j) - feed(k,m))< eee)
                    % if a value is similar to the corresponding
                    % "feed-coordinate" overwrite it!
                    points(k,j) = feed(k,m);
                    end
                end
            end
            if(abs(points(k,j) - points(k,i))< eee)
                % if the difference is less than "eee" switch the
                % corresponding cell in "temp" to "TRUE"
                temp(k,j) = true; 
            end
            
        end
    end
    for k = 1 : size (points,1) % for each row in "points"
        % replace the similar values with the mean of those values
        points(k,temp(k,:))= mean(points(k,temp(k,:))) ;
       
    end

end

 
 
 