function points = FindPropulsion(accelometer, HS , TO)
% HS
% TO

% if(HS(1) < TO(1))
    points = ones(1,length(HS));
    for k = 1:length(HS)
        bool = 1;
        num_positive = 0;
        HSn = HS(k);
        if (isnan(round(((HS(k)+TO(k))/2))) == 0)
            for i =  round(((HS(k)+TO(k))/2)):3:TO(k)
                
                if bool
                    first_sample = [i:i+2];
                    second_sample = [i+3:i+5];
                    temp_mean1 = mean(accelometer(first_sample));
                    temp_STD1 = std(accelometer(first_sample));
                    temp_mean2 = mean(accelometer(second_sample));
                    temp_STD2 = std(accelometer(second_sample));
                    %                 result(round(i/3)+1,:) = [i temp_mean2-temp_mean1/s temp_mean2-temp_mean1 temp_STD2-temp_STD1 temp_mean1 temp_mean2 temp_STD1 temp_STD2];
                    temp_mean2-temp_mean1;
                    if temp_mean2-temp_mean1 > .01
                        num_positive = num_positive + 1;
                    else num_positive =0;
                    end
                    if num_positive > 3
                        point_interest = i - 9;
                        points(k) = FindMin(accelometer,point_interest);
                        bool = 0;
                    end
                end
            end
        end
        
    end
end


function m = FindMin(x,poin)
    [value, m] = min(x(poin-5:poin+5));
    m = poin - 5 + m;
end