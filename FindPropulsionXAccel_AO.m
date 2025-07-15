function points = FindPropulsionXAccel(DG_X, HS , TO)
% % Modified by Dr. Abderrahman Ouattas on 03/27/2024
% To undestand what specific mehtod was used to idenfity the start of
% propulsion
% Email: Abderrahman.Ouattas@bcm.edu

% % The specific method used is relatively simple and straighforward; 
% They identified the start of propulsion phase the local maximum of shin
% angular velocity between HS and TO. 

% % However That does not correspond to the mid-stance time that was
% referenced in the drafted paper (cited by Winter et al., 1990) 
% Here is how they defined the start of propulsion phase: 
% The propulsion phase was defined as the beginning point of the 
% second rise in the pressure, the time that force switch to the forefoot

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

% k=5
% plot(DG_X(HS(k):TO(k)))
% plot(DG_X)
% hold on 
% plot(p,'*r')
% hold off

function m = FindMin(x,poin)
    [value, m] = min(x(poin-5:poin+5));
    m = poin - 5 + m;
end