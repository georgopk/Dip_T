function points = sp_round(points,eee,feed)


for i=1:length (points)
    temp = false; 
    for j=1:length (points)
        if (abs(points(1,j) - feed(1))< eee)
            points(1,j) = feed(1);
        elseif(abs(points(1,j) - points(1,i))< eee)
            temp(1,j) = true; 
        end
        if (abs(points(2,j) - feed(2))< eee)
            points(1,j) = feed(1);
        elseif (abs(points(2,j) - points(2,i))< eee)
            temp(2,j) = true;
        end
    end
     
    points(1,temp(1,:))= mean(points(1,temp(1,:))) ;
    points(2,temp(2,:))= mean(points(2,temp(2,:))) ;
end