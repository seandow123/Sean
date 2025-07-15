function [pro_point L] = FindPropulsionPressure(Pressure, HS , TO)
% HS
% TO
        HS = HS - 30;
        TO = TO + 30;
        % if(HS(1) < TO(1))
        pro_point = ones(1,length(HS));

        for k = 1:length(HS)
            k
            B = Pressure(HS(k):TO(k));
            L = findpeaks(B,'MinPeakDistance',size(B,1)/2);
            S1 = find(B == L(1));
            S2 = find(B == L(2))
            K = min(B(S1:S2));
            Local_min = find(B(S1:S2) == K)
            Bdiff = diff(B);
            Bdiff_max = max(Bdiff(Local_min(1):S2));
            i = Local_min(1);
            while i < S2
                temp = mean(Bdiff(i:i+2,1));
                if(temp > Bdiff_max/4)
                    pro_point(k) = i + HS(k)
                    i = S2;
                end
                i = i + 2;
            end
        end
end


