
function res = GaitAnalyze(ss, h)

% function res = GaitAnalyze(ss, h)

% res = GaitAnalyze(s, h);
% Does a gait analysis of the ss signal. h parameter is the height of the
% subject. The order of the channels is:
% s(:,1) = right thigh
% s(:,2) = right shank
% s(:,3) = left shank
% s(:,4) = left thigh



clear global shanks thighs;
%%%%%%%%%%%%%%%%%%%%%%% drift cancelation %%%%%%%%%%%%%%%%%%%%
global shanks thighs;
shanks = [];
thighs = [];

[m,n] = size(ss)

if n == 2
    shanks(:,1) = RemoveDrift(ss(:,1));
    shanks(:,2) = RemoveDrift(ss(:,2));
    thighs = [];
elseif n == 3
    shanks(:,1) = RemoveDrift(ss(:,2));
    shanks(:,2) = RemoveDrift(ss(:,3));
    thighs(:,1) = RemoveDrift(ss(:,1));
    thighs(:,2) = 0;
elseif n == 4
    shanks(:,1) = RemoveDrift(ss(:,2));
    shanks(:,2) = RemoveDrift(ss(:,3));
    thighs(:,1) = RemoveDrift(ss(:,1));
    thighs(:,2) = RemoveDrift(ss(:,4));
end

C = CalcCycles;

if n == 2
    res = CalcAll2(C);
elseif n == 3
    L1=.245*h;
    L2=.246*h;
    L1 = L1 / 100;
    L2 = L2 / 100;
    res = CalcAll3(C, L1, L2);
elseif n == 4
    L1=.245*h;
    L2=.246*h;
    L1 = L1 / 100;
    L2 = L2 / 100;
    res = CalcAll4(C, L1, L2);
end

function Cycle = CalcAll2(C);
global shanks thighs;
Cycle = C;

for i=1:length(Cycle)
    %%%%%%%%%%%%%%%%%%%%%%% Temporal parameters %%%%%%%%%%%%%%%%%%%%%
    Cycle(i).GctR = Cycle(i).end - Cycle(i).HsR;
    
    % Rajout pour avoir le Gct de la jambe gauche
%     Cycle(i).GctL = Cycle(i).endL - Cycle(i).HsL;
    if i < length(Cycle)
        Cycle(i).GctL = Cycle(i+1).HsL - Cycle(i).HsL;
    else
        Cycle(i).GctL = NaN;
    end
    Cycle(i).GctL = Cycle(i).GctL / 200;
    % %%%%%
    
    Cycle(i).SwingR = 100 * (Cycle(i).end - Cycle(i).ToR) / Cycle(i).GctR;
    Cycle(i).SwingL = 100 * (Cycle(i).HsL - Cycle(i).ToL) / Cycle(i).GctR;
    Cycle(i).StanceR = 100 * (Cycle(i).ToR - Cycle(i).HsR) / Cycle(i).GctR;
    Cycle(i).StanceL = 100 * (Cycle(i).end - Cycle(i).HsL + Cycle(i).ToL - Cycle(i).HsR) / Cycle(i).GctR;
    Cycle(i).IDS = 100 * (Cycle(i).ToL - Cycle(i).HsR) / Cycle(i).GctR;
    Cycle(i).TDS = 100 * (Cycle(i).ToR - Cycle(i).HsL) / Cycle(i).GctR;
    Cycle(i).DS = Cycle(i).IDS + Cycle(i).TDS;
    Cycle(i).Limp = abs(Cycle(i).IDS - Cycle(i).TDS);
    Cycle(i).GctR = Cycle(i).GctR / 200;
    
    %%%%%%%%%%%%%%%%%%%%%%% Range of Angels %%%%%%%%%%%%%%%%%%%%%
    Cycle(i).ThighR = NaN;
    Cycle(i).ThighL = NaN;
    Cycle(i).KneeR = NaN;
    Cycle(i).KneeL = NaN;
    if isnan(Cycle(i).HsR) | isnan(Cycle(i).end)
        Cycle(i).ShankR = NaN;
        Cycle(i).ShankL = NaN;
    else
        Cycle(i).ShankR = range(naninteg(2, Cycle(i).HsR, Cycle(i).end));
        Cycle(i).ShankL = range(naninteg(3, Cycle(i).HsR, Cycle(i).end));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%% Spatial parameters %%%%%%%%%%%%%%%%%%%
    
    Cycle(i).StrideR = .491*sqrt(2*(1-cos(Cycle(i).ShankR*pi/180))); %%% stride length
    Cycle(i).SpeedR = Cycle(i).StrideR / Cycle(i).GctR; %% stride speed
    
    Cycle(i).StrideL = .491*sqrt(2*(1-cos(Cycle(i).ShankL*pi/180)));
    Cycle(i).SpeedL = Cycle(i).StrideL / Cycle(i).GctR;
        
    %%%%%%%%%%%%%%%%%%%%%%% Midswing Speeds %%%%%%%%%%%%%%%%%%%%
    if isnan(Cycle(i).HsR) | isnan(Cycle(i).end)
        Cycle(i).PeakSwingSpeedR = NaN;
        Cycle(i).PeakSwingSpeedL = NaN;
    else
        Cycle(i).PeakSwingSpeedR = max(shanks(Cycle(i).HsR:Cycle(i).end, 1));
        Cycle(i).PeakSwingSpeedL = max(shanks(Cycle(i).HsR:Cycle(i).end, 2));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Cycle = CalcAll3(C, L1, L2);
global shanks thighs;
Cycle = C;

for i=1:length(Cycle)
    %%%%%%%%%%%%%%%%%%%%%%% Temporal parameters %%%%%%%%%%%%%%%%%%%%%
    Cycle(i).GctR = Cycle(i).end - Cycle(i).HsR;
    
    % Rajout pour avoir le Gct de la jambe gauche
    Cycle(i).GctL = Cycle(i).endL - Cycle(i).HsL;
    Cycle(i).GctL = Cycle(i).GctL / 200;
    % %%%%%
    
    Cycle(i).SwingR = 100 * (Cycle(i).end - Cycle(i).ToR) / Cycle(i).GctR;
    Cycle(i).SwingL = 100 * (Cycle(i).HsL - Cycle(i).ToL) / Cycle(i).GctR;
    Cycle(i).StanceR = 100 * (Cycle(i).ToR - Cycle(i).HsR) / Cycle(i).GctR;
    Cycle(i).StanceL = 100 * (Cycle(i).end - Cycle(i).HsL + Cycle(i).ToL - Cycle(i).HsR) / Cycle(i).GctR;
    Cycle(i).IDS = 100 * (Cycle(i).ToL - Cycle(i).HsR) / Cycle(i).GctR;
    Cycle(i).TDS = 100 * (Cycle(i).ToR - Cycle(i).HsL) / Cycle(i).GctR;
    Cycle(i).DS = Cycle(i).IDS + Cycle(i).TDS;
    Cycle(i).Limp = abs(Cycle(i).IDS - Cycle(i).TDS);
    Cycle(i).GctR = Cycle(i).GctR / 200;
    
    %%%%%%%%%%%%%%%%%%%%%%% Spatial parameters %%%%%%%%%%%%%%%%%%%
    
    %%% Right leg %%%
    %%%%%% right swing 
    alpha = naninteg(1, Cycle(i).ToR, Cycle(i).end);
    beta = naninteg(2, Cycle(i).ToR, Cycle(i).end);               
    d1d2 = SolveGaitModel(alpha(end) - alpha(1), beta(end)-beta(1), L1, L2);
    %%%%%% left stance
    if i == 1       
        Cycle(i).StrideR = NaN;
        Cycle(i).SpeedR = NaN;
    else
        %alpha = naninteg(4, Cycle(i-1).HsL, Cycle(i).ToL);
        alpha = naninteg(1, Cycle(i).HsR, Cycle(i).ToR);
        beta = naninteg(3, Cycle(i-1).HsL, Cycle(i).ToL);
        d3 = SolveGaitModel(beta(end)-beta(1), alpha(end) - alpha(1), L2, L1);
        Cycle(i).StrideR = d1d2 + d3;
        Cycle(i).SpeedR = Cycle(i).StrideR / Cycle(i).GctR;
    end
%     %%% Left leg %%%
%     %%%%%% Left swing
%     %alpha = naninteg(4, Cycle(i).ToL, Cycle(i).HsL);
%     alpha = naninteg(1, Cycle(i).ToR, Cycle(i).end);
%     beta = naninteg(3, Cycle(i).ToL, Cycle(i).HsL);
%     d1d2 = SolveGaitModel(alpha(end) - alpha(1), beta(end)-beta(1), L1, L2);
%     %%%%%% Rigth stance
%     alpha = naninteg(1, Cycle(i).HsR, Cycle(i).ToR);
%     beta = naninteg(2, Cycle(i).HsR, Cycle(i).ToR);
%     d3 = SolveGaitModel(beta(end)-beta(1), alpha(end) - alpha(1), L2, L1);
%     Cycle(i).StrideL = d1d2 + d3;
%     Cycle(i).SpeedL = Cycle(i).StrideL / Cycle(i).GctR;
    Cycle(i).StrideL = NaN;
    Cycle(i).SpeedL = NaN;
    
    %%%%%%%%%%%%%%%%%%%%%%% Range of Angels %%%%%%%%%%%%%%%%%%%%%
    if isnan(Cycle(i).HsR) | isnan(Cycle(i).end)
        Cycle(i).ShankR = NaN;
        Cycle(i).ShankL = NaN;
        Cycle(i).ThighR = NaN;
        Cycle(i).ThighL = NaN;
        Cycle(i).KneeR = NaN;
        Cycle(i).KneeL = NaN;
    else
        Cycle(i).ShankR = range(naninteg(2, Cycle(i).HsR, Cycle(i).end));
        Cycle(i).ShankL = range(naninteg(3, Cycle(i).HsR, Cycle(i).end));
        Cycle(i).ThighR  = range(naninteg(1, Cycle(i).HsR, Cycle(i).end));
        Cycle(i).ThighL  = NaN;
        Cycle(i).KneeR = range(naninteg(1, Cycle(i).HsR, Cycle(i).end) - naninteg(2, Cycle(i).HsR, Cycle(i).end));
        Cycle(i).KneeL = NaN;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%% Midswing Speeds %%%%%%%%%%%%%%%%%%%%
    if isnan(Cycle(i).HsR) | isnan(Cycle(i).end)
        Cycle(i).PeakSwingSpeedR = NaN;
        Cycle(i).PeakSwingSpeedL = NaN;
    else
        Cycle(i).PeakSwingSpeedR = max(shanks(Cycle(i).HsR:Cycle(i).end, 1));
        Cycle(i).PeakSwingSpeedL = max(shanks(Cycle(i).HsR:Cycle(i).end, 2));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Cycle = CalcAll4(C, L1, L2);
global shanks thighs;
Cycle = C;

for i=1:length(Cycle)
    %%%%%%%%%%%%%%%%%%%%%%% Temporal parameters %%%%%%%%%%%%%%%%%%%%%
    Cycle(i).GctR = Cycle(i).end - Cycle(i).HsR;
    
    % Rajout pour avoir le Gct de la jambe gauche
    Cycle(i).GctL = Cycle(i).endL - Cycle(i).HsL;
    Cycle(i).GctL = Cycle(i).GctL / 200;
    % %%%%%
    
    Cycle(i).SwingR = 100 * (Cycle(i).end - Cycle(i).ToR) / Cycle(i).GctR;
    Cycle(i).SwingL = 100 * (Cycle(i).HsL - Cycle(i).ToL) / Cycle(i).GctR;
    Cycle(i).StanceR = 100 * (Cycle(i).ToR - Cycle(i).HsR) / Cycle(i).GctR;
    Cycle(i).StanceL = 100 * (Cycle(i).end - Cycle(i).HsL + Cycle(i).ToL - Cycle(i).HsR) / Cycle(i).GctR;
    Cycle(i).IDS = 100 * (Cycle(i).ToL - Cycle(i).HsR) / Cycle(i).GctR;
    Cycle(i).TDS = 100 * (Cycle(i).ToR - Cycle(i).HsL) / Cycle(i).GctR;
    Cycle(i).DS = Cycle(i).IDS + Cycle(i).TDS;
    Cycle(i).Limp = abs(Cycle(i).IDS - Cycle(i).TDS);
    Cycle(i).GctR = Cycle(i).GctR / 200;
    
    %%%%%%%%%%%%%%%%%%%%%%% Spatial parameters %%%%%%%%%%%%%%%%%%%
    
    %%%%%% right swing 
    alpha = naninteg(1, Cycle(i).ToR, Cycle(i).end);
    beta = naninteg(2, Cycle(i).ToR, Cycle(i).end);               
    d1d2 = SolveGaitModel(alpha(end) - alpha(1), beta(end)-beta(1), L1, L2);
    %%%%%% left stance
    if i == 1       
        Cycle(i).StrideR = NaN;
        Cycle(i).SpeedR = NaN;
    else
        alpha = naninteg(4, Cycle(i-1).HsL, Cycle(i).ToL);
        beta = naninteg(3, Cycle(i-1).HsL, Cycle(i).ToL);
        d3 = SolveGaitModel(beta(end)-beta(1), alpha(end) - alpha(1), L2, L1);
        Cycle(i).StrideR = d1d2 + d3;
        Cycle(i).SpeedR = Cycle(i).StrideR / Cycle(i).GctR;
    end
    %%%%%% Left swing
    alpha = naninteg(4, Cycle(i).ToL, Cycle(i).HsL);
    beta = naninteg(3, Cycle(i).ToL, Cycle(i).HsL);
    d1d2 = SolveGaitModel(alpha(end) - alpha(1), beta(end)-beta(1), L1, L2);
    %%%%%% Rigth stance
    alpha = naninteg(1, Cycle(i).HsR, Cycle(i).ToR);
    beta = naninteg(2, Cycle(i).HsR, Cycle(i).ToR);
    d3 = SolveGaitModel(beta(end)-beta(1), alpha(end) - alpha(1), L2, L1);
    Cycle(i).StrideL = d1d2 + d3;
    Cycle(i).SpeedL = Cycle(i).StrideL / Cycle(i).GctL;
    
    %%%%%%%%%%%%%%%%%%%%%%% Range of Angels %%%%%%%%%%%%%%%%%%%%%
    if isnan(Cycle(i).HsR) | isnan(Cycle(i).end)
        Cycle(i).ShankR = NaN;
        Cycle(i).ShankL = NaN;
        Cycle(i).ThighR = NaN;
        Cycle(i).ThighL = NaN;
        Cycle(i).KneeR = NaN;
        Cycle(i).KneeL = NaN;
    else
        Cycle(i).ShankR = range(naninteg(2, Cycle(i).HsR, Cycle(i).end));
        Cycle(i).ShankL = range(naninteg(3, Cycle(i).HsR, Cycle(i).end));
        Cycle(i).ThighR  = range(naninteg(1, Cycle(i).HsR, Cycle(i).end));
        Cycle(i).ThighL  = range(naninteg(4, Cycle(i).HsR, Cycle(i).end));
        Cycle(i).KneeR = range(naninteg(1, Cycle(i).HsR, Cycle(i).end) - naninteg(2, Cycle(i).HsR, Cycle(i).end));
        Cycle(i).KneeL = range(naninteg(4, Cycle(i).HsR, Cycle(i).end) - naninteg(3, Cycle(i).HsR, Cycle(i).end));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%% Midswing Speeds %%%%%%%%%%%%%%%%%%%%
    if isnan(Cycle(i).HsR) | isnan(Cycle(i).end)
        Cycle(i).PeakSwingSpeedR = NaN;
        Cycle(i).PeakSwingSpeedL = NaN;
    else
        Cycle(i).PeakSwingSpeedR = max(shanks(Cycle(i).HsR:Cycle(i).end, 1));
        Cycle(i).PeakSwingSpeedL = max(shanks(Cycle(i).HsR:Cycle(i).end, 2));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = SolveGaitModel(alpha, beta, l1, l2)
if isnan(alpha) | isnan(beta)
    d = NaN;
else 
    alpha = abs(alpha * pi / 180);
    beta = abs(beta * pi / 180);
    
    gamma = (pi - alpha) / 2;
    phi = (pi - 2*beta + alpha) / 2;
    
    d_ = l1 * sqrt(2*(1 - cos(alpha)));
    
    M1 = sin(phi) * d_ / sin(beta);
    M2 = sin(gamma) * d_ / sin(beta);
    
    d = sqrt((l2+M2)^2 + (l2+M1)^2 - 2*(l2+M2)*(l2+M1)*cos(beta));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Cycle = CalcCycles
global shanks;
[HsR,ToR] = DetectHsTo(shanks(:,1));
[HsL,ToL] = DetectHsTo(shanks(:,2));

% --------------- Added by bijan 2010
% if HsL==HsR
%     HsL = ToR - 30;
%     ToL = HsR +30;
% end

% PlotHsTo(shanks, ToR, HsR, ToL, HsL);

%%%
% Rajout pour éviter des NaN qui n'ont pas lieu d'etre
HsR(isnan(HsR))=[];
ToL(isnan(ToL))=[];
HsL(isnan(HsL))=[];
ToR(isnan(ToR))=[];
%%%

for i=1:length(HsR)-1
    Cycle(i).HsR = HsR(i);
    Cycle(i).ToR = FindFirstAfter(Cycle(i).HsR, ToR);
    
    Cycle(i).ToL = FindFirstAfter(Cycle(i).HsR, ToL);
    if Cycle(i).ToL > Cycle(i).ToR
        Cycle(i).ToL = NaN;
    end
    
    Cycle(i).HsL = FindFirstAfter(Cycle(i).ToL, HsL);
    if Cycle(i).HsL > Cycle(i).ToR
        Cycle(i).HsL = NaN;
    end
    
    Cycle(i).end = HsR(i+1);
    if Cycle(i).HsR == Cycle(i).end
        Cycle(i).end = NaN;
    end
    
%     % Ancien Rajout au programme pour avoir le gct pour la jambe gauche
%     Cycle(i).endL2 = FindFirstAfter((Cycle(i).HsL+1), HsL);
%     if Cycle(i).HsL == Cycle(i).endL
%         Cycle(i).endL2 = NaN;
%     end
%     % fin Rajout %%%%%%  
end

%%% Nouveau Rajout au programme pour avoir le gct pour la jambe gauche %%%
endL = [Cycle(2:end).HsL];
%endL = [Cycle(1:end).HsL];
abc=endL(end);
bcd=find(HsL>abc);
if isempty(bcd)
    endL=[endL NaN];
else
    if HsL(bcd(1))>Cycle(end).end
        endL=[endL HsL(bcd(1))];
    else
        endL=[endL NaN];
    end
end
for j=1:length(Cycle)
    Cycle(j).endL=endL(j);
end


%%% patch %%%%
if isnan(Cycle(end).ToR) & isnan(Cycle(end).ToL) & isnan(Cycle(end).HsL)
    Cycle(end) = [];
end

%%% filter %%%
k = 1;
while k < length(Cycle)
    if Cycle(k).ToR == Cycle(k+1).ToR;
        Cycle(k) = [];
    else
        k = k + 1;
    end
end
% pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function i = naninteg(n, a, b)
global shanks thighs;

if isnan(a) | isnan(b)
    i = NaN;
else
    switch n
        case 1 
            i = cumtrapz(thighs(a:b, 1))/200;
        case 2
            i = cumtrapz(shanks(a:b, 1))/200;
        case 3
            i = cumtrapz(shanks(a:b, 2))/200;
        case 4
            i = cumtrapz(thighs(a:b, 2))/200;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function PlotHsTo(shank, ToR, HsR, ToL, HsL)
% figure
% subplot(211); plotpeak(shank(:,1), ToR, HsR);
% title('Shank Angular Rate [DEG/SEC]')
% ylabel('Right')
% 
% subplot(212); plotpeak(shank(:,2), ToL, HsL);
% ylabel('Left');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function plotpeak(x, p, p2, fs)
% 
% n = find(isnan(p));
% p(n) = [];
% 
% if nargin == 2
%     fs = 200;
%     plot((1:length(x))/fs, x, p/fs, x(p), 'o')
%     xlabel('sec');    
% else
%     if nargin == 3
%         fs = 200;
%     end
%     n = find(isnan(p2));
%     p2(n) = [];
%     plot((1:length(x))/fs, x, p/fs, x(p), 'o', p2/fs, x(p2), 'x')
%     xlabel('sec');
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Hs,To]=DetectHsTo(gyr);
% 	[To,Hs]=DetectHsToF(gyr)
% Detects the Hs and To peaks in the shank's signal gyr

[min_peak, max_peak] = FindMinMax(gyr);
%%%%% detect MidSwings %%%%%%%%%%%%
 max_peak = max_peak(gyr(max_peak) > 40);
%max_peak = max_peak(gyr(max_peak) > 90);
n = length(max_peak);
i = 1;
MidSwing = [];
while i<=n 
    k = i;
    while k <= n-1
        if max_peak(k+1) - max_peak(k) < 100
            k = k + 1;
        else
            break;
        end
    end
    MidSwing = [MidSwing; max_peak(fix((i+k)/2))];
    i = k+1;
end

%%%%% detect Heel-strike %%%%%%%%%%%
n = length(MidSwing);
nmp = length(min_peak);
Hs = zeros(length(MidSwing),1);
for i=1:n
    [peak, idx] = FindFirstAfter(MidSwing(i), min_peak);
    Hs(i) = NaN;
    if idx > 0
        while(idx <= nmp)
            if gyr(min_peak(idx)) >= 0 
                idx = idx + 1;
            elseif min_peak(idx) - MidSwing(i) < 200
                Hs(i) = min_peak(idx);
                break;
            end
        end
    end
end

%%%%% detect Toe-Off %%%%%%%%%%%
% load D:\MATLAB6p5\work\AE\lowpass30
Num=[1.3347068e-003	  1.9320946e-002	  1.3286959e-003	 -4.3606951e-003	 -8.1735044e-003	 -4.1616412e-003	  5.0501299e-003	  1.0787990e-002	  6.2937331e-003	 -5.7450590e-003	 -1.4144170e-002	 -9.3237067e-003	  6.3980351e-003	  1.8637193e-002	  1.3717633e-002	 -6.9751517e-003	 -2.5097091e-002	 -2.0654645e-002	  7.4576492e-003	  3.5605168e-002	  3.3188356e-002	 -7.8069825e-003	 -5.7388784e-002	 -6.3566699e-002	  8.0281057e-003	  1.4117799e-001	  2.7127322e-001	  3.2522741e-001	  2.7127322e-001	  1.4117799e-001	  8.0281057e-003	 -6.3566699e-002	 -5.7388784e-002	 -7.8069825e-003	  3.3188356e-002	  3.5605168e-002	  7.4576492e-003	 -2.0654645e-002	 -2.5097091e-002	 -6.9751517e-003	  1.3717633e-002	  1.8637193e-002	  6.3980351e-003	 -9.3237067e-003	 -1.4144170e-002	 -5.7450590e-003	  6.2937331e-003	  1.0787990e-002	  5.0501299e-003	 -4.1616412e-003	 -8.1735044e-003	 -4.3606951e-003	  1.3286959e-003	  1.9320946e-002	  1.3347068e-003];
sf = filtfilt(Num, 1, gyr);
[min_peak, max_peak] = FindMinMax(sf);

n = length(MidSwing);
To = zeros(n,1);
for i=1:n
    [peak, idx] = FindFirstBefore(MidSwing(i), min_peak);
    To(i) = NaN;
    while(idx > 1)
        if sf(min_peak(idx)) >= -20 
            idx = idx - 1;
        else
            if MidSwing(i)-min_peak(idx) < 200
                To(i) = min_peak(idx);
            end
            break;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mi, ma] = FindMinMax(gyr, thr)

d = diff(gyr);
f = find( d(2:end) .* d(1:end-1) <= 0);
f = f + 1;

mi = f(d(f)>=0);
ma = f(d(f)<0);

if nargin == 2
    ma = ma(gyr(ma) > thr);
    mi = mi(gyr(mi) < -thr);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t, k] = FindFirstBefore(x, a)
k = 0;
t = NaN;
if isnan(x) | isempty(a)
    return;
elseif a(1) >= x
    return;
else
    k = binarySearch(x, a);
    if k == 0
        return
    else
        t = a(k);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t, i] = FindFirstAfter(x, a)
%   [t, i] = FindFirstAfter(x, a)
% Finds the first elemet of a, bigger than or equal to x
i = 0;
t = NaN;
if isnan(x) | isempty(a)
    return;
elseif a(end) <= x
    return;
else
    k = binarySearch(x, a);
    if k == 0
        i = k;
        return
    end
    for i=k:length(a)
        if a(i) >= x
            t = a(i);
            return;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = binarySearch(x, a)
%       m = binarySearch(x, a)
% Returns the lower bound of index m in the sorted array a where you can find x
n1 = 1;
n2 = length(a);

if a(1) > x
    m = 0;
end

while(n1 < n2)
    m = fix((n1+n2) / 2);
    if m == n1
        if a(n2) < x
            m = n2;
        end
        return
    elseif a(m) > x
        n2 = m;
    elseif a(m) == x
        return
    else
        n1 = m;
    end
end
m = n1;        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = RemoveDrift(a)
% s = filtfilt([1,-1],[1,-0.995],a);
s=a;
% s = a;