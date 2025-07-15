function points = FindPropulsionXAccel(DG_X, HS , TO)
% DG_X = DG_x
% HS = [res_St.HsR]
% TO = [res_St.ToR]

% if(HS(1) < TO(1))
    points = ones(1,length(HS));
    for k = 1:length(HS)
        bool = 1;
        num_positive = 0;
        HSn = HS(k);
        if (isnan(round(((HS(k)+TO(k))/2))) == 0)
            points(k) = find(DG_X(HS(k):TO(k)) == max(DG_X(HS(k):TO(k)))) + HS(k);
        end
        
    end
end

k=4
plot(DG_X(HS(k):TO(k)))
plot(DG_X)
hold on 
plot(p,'*r')
hold off

function m = FindMin(x,poin)
    [value, m] = min(x(poin-5:poin+5));
    m = poin - 5 + m;
end