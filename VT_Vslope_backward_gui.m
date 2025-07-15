function [VT1best,VT1m1,VT1m2, VT1b1, VT1b2, kbest,kbestMSE,RSSmin, MSEmin, R2best,MSE] = VT_Vslope_backward_gui(VO2,VCO2,options)

    
    
    % Fit the two compartment model to estimate VT1 in the data set 
    % =======================================================================
    % The real-time version of this program will need to work breath-by-breath,  
    % or, when a mixing box is used, sample-by-sample. Here we begin developing
    % that approach using the V-slope method, in which VT1 is identified as
    % the point where the ratio of VCO2/VO2 suddenly increases. 
    % (cf. Gaskill et al. MSSE 33(11) 1841-1848, 2001). 
    %=========================================================================
    % The V-slope method examines the VCO2 vs VO2 data relationship, which is
    % linear (VCO2=m1*VO2+b1) below VT1 and linear (VCO2=m2*VO2+b2) above VT1,
    % with slopes m2 > m1.
    %--------------------------------------------------------------------------
    % The algorithm used creates a set of potential VT1 points by sequentially
    % (breath-by-breath/sample-by-sample) fitting two straight lines to the data
    % using a least squares techniques (polyfit). The VT1 candidate for each  
    % pair of lines fit is identified as the join point between those two
    % lines. When all the fits have been made, the estimated VT1 is determined 
    % as the one that minimizs the mean square error of the overall fit.
    %--------------------------------------------------------------------------
    % Note: this assumes the full data set is available (post hoc analysis)
    % but that it is truncated below VT2 (or there is no VT2 for this trial)
    % these lines will need to be modified for breath-by-breath analysis
    %--------------------------------------------------------------------------
    
    % first set limits for locating VT1. Leave enough data above and below
    % to successfully fit lines to those regions.
    iVT1_first = 2;             % start after 2 samples of selected arrays
    iVT1_last = length(VO2)-2;  % end 2 samples before end of selected arrays
    
    % now initialize some curve fit parameters
    m1 = zeros(length(VO2), 1); % slopes below VT1
    b1 = zeros(length(VO2), 1); % intercepts below VT1
    m2 = zeros(length(VO2), 1); % slopes above VT1
    b2 = zeros(length(VO2), 1); % intercepts above VT1
    VT1 = zeros(length(VO2), 1); % join points--VO2 coordinate
    VCO2atVT1 = zeros(length(VO2), 1); % join points--VCO2 coordinate
%     xplot = linspace (min(VO2), max(VO2), length(VO2)); 
    VCO2fit = zeros(length(VO2), length(VO2)); 
    yresid = zeros(length(VO2), length(VO2)); 
    RSS = zeros(length(VO2), 1); % sum of squared errors of residuals
    TSS = zeros(length(VO2), 1); % total sum of square
    R2 = zeros(length(VO2), 1); % goodness of fit
    
    MSE = zeros(length(VO2), 1); % mean square errors of the fitted data
  
    switch(options)
        case 1
            % and set up a loop to find VT1 and MSE for each breath/sample in the range 
            VO2  = VO2';
            for k = iVT1_first: iVT1_last
                
                % find the slope and intercept of the least squares best fit lines
                % through the VCO2 vs VO2 data below and above the current VT1
                
                P1 = polyfit(VO2(1:k),VCO2(1:k), 1);
                m1(k) = P1(1);  
                b1(k) = P1(2);

                P2 = polyfit(VO2(k:length(VO2)),VCO2(k:length(VO2)), 1); 
                m2(k) = P2(1);
                b2(k) = P2(2);

                % now find and save the current VT1 (join point between best fit lines)    
                VT1(k) = (b2(k) - b1(k))/(m1(k) - m2(k));
                VCO2atVT1(k) = m1(k)*VT1(k)+b1(k);

                % and concatenate the curves to create a best fit model
                for j = (1:length(VO2));
        %             if polyval(P1, VO2(j)) < VT1(k);
        %                 VCO2fit(k,j) = polyval(P1, VO2(j));
        %             else
        %                 VCO2fit(k,j) = polyval(P2, VO2(j));
        %             end
                    if j < k;
                        VCO2fit(k,j) = polyval(P1, VO2(j));
                    else
                        VCO2fit(k,j) = polyval(P2, VO2(j));
                    end
                end

                % now compute the fit errors (see Jones et al 1984 for details)
                yresid(k,:) = VCO2' - VCO2fit(k,:);
                RSS(k) = sum(yresid(k,:).^2);  % residual sum of squares
                MSE(k)=RSS(k)./length(VCO2);
                %LP
                meanVCO2=mean(VCO2);
                TSS(k) = sum((VCO2'-meanVCO2).^2);  % total sum of squares
                R2(k)=1-RSS(k)./TSS(k);
                n=length(VO2);
                p=1;
                adjR2(k)=1-(1-R2(k))*(n-1)/(n-p-1);

        %         figure(1)
        %         clf;
        %         plot( VCO2fit(k,:),'b.')
        %         hold on
        %         plot( VCO2,'r')
        %         hold on
        %         plot(k,VCO2(k),'k+')
        %         k

            end
                
        case 2
            % and set up a loop to find VT1 and MSE for each breath/sample in the range 
            for k = iVT1_first: iVT1_last

                % Alternative approach: we force the linear fit to go through the
                % candidate VT point
                %% m1
                opt = optimset('LargeScale','off','Display','off');
                x=VO2(1:k);
                y=VCO2(1:k);
                x0 = VO2(k);
                y0 = VCO2(k);
                x = x(:); %reshape the data into a column vector
                y = y(:);
                n = 1; % Degree of polynomial to fit
                V=[];
                V(:,n+1) = ones(length(x),1,class(x));
                for j = n:-1:1
                    V(:,j) = x.*V(:,j+1);
                end
                C = V;  % 'C' is the Vandermonde matrix for 'x'
                d = y;  % 'd' is the vector of target values, 'y'.
                A = []; % There are no inequality constraints in this case
                b = [];
                Aeq = x0.^(n:-1:0); % We use linear equality constraints to force the curve to hit the required point. In this case, 'Aeq' is the Vandermoonde matrix for 'x0' 
                beq = y0; % and 'beq' is the value the curve should take at that point
                P1 = lsqlin( C, d, A, b, Aeq, beq,[],[],[],opt);
                m1(k) = P1(1);  
                b1(k) = P1(2);
                P1(1) = VO2(1);
                P1(2) = VCO2(1);
                

                %% m2
                x=VO2(k:length(VO2));
                y=VCO2(k:length(VO2));
                x0 = VO2(k);
                y0 = VCO2(k);
                x = x(:); %reshape the data into a column vector
                y = y(:);
                n = 1; % Degree of polynomial to fit
                V=[];
                V(:,n+1) = ones(length(x),1,class(x));
                for j = n:-1:1
                    V(:,j) = x.*V(:,j+1);
                end
                C = V;  % 'C' is the Vandermonde matrix for 'x'
                d = y;  % 'd' is the vector of target values, 'y'.
                A = []; % There are no inequality constraints in this case
                b = [];
                Aeq = x0.^(n:-1:0); % We use linear equality constraints to force the curve to hit the required point. In this case, 'Aeq' is the Vandermoonde matrix for 'x0' 
                beq = y0; % and 'beq' is the value the curve should take at that point
                P2 = lsqlin( C, d, A, b, Aeq, beq,[],[],[],opt);
                m2(k) = P2(1);  
                b2(k) = P2(2);
                P2(1) = VO2(end);
                P2(2) = VCO2(end);

                % now find and save the current VT1 (join point between best fit lines)    
                VT1(k) = (b2(k) - b1(k))/(m1(k) - m2(k));
                VCO2atVT1(k) = m1(k)*VT1(k)+b1(k);

                % and concatenate the curves to create a best fit model
                for j = (1:length(VO2));
                    if j < k;
                        VCO2fit(k,j) = polyval(P1, VO2(j));
                    else
                        VCO2fit(k,j) = polyval(P2, VO2(j));
                    end
                end

                % now compute the fit errors (see Jones et al 1984 for details)
                yresid(k,:) = VCO2' - VCO2fit(k,:);
                RSS(k) = sum(yresid(k,:).^2);  % residual sum of squares
                MSE(k)=RSS(k)./length(VCO2);
                meanVCO2=mean(VCO2);
                TSS(k) = sum((VCO2'-meanVCO2).^2);  % total sum of squares
                R2(k)=1-RSS(k)./TSS(k);
                n=length(VO2);
                p=1;
                adjR2(k)=1-(1-R2(k))*(n-1)/(n-p-1);

        %         figure(1)
        %         clf;
        %         plot( VCO2fit(k,:),'b.')
        %         hold on
        %         plot( VCO2,'r')
        %         hold on
        %         plot(k,VCO2(k),'k+')
        %         k

            end
    
    end
    RSS(iVT1_first);
    RSSmin = min(RSS(iVT1_first:iVT1_last));
    kbest = find(RSS == RSSmin,1);
    
    MSE(iVT1_first)
    MSEmin = min(MSE(iVT1_first:iVT1_last));
    portion = MSE(iVT1_first)/MSEmin
%     if portion < 4.0
%         MSEmin = MSE(iVT1_first);
%     end
    kbestMSE = find(MSE == MSEmin,1);

    R2best=R2(kbest);
    find(R2==max(R2),1);
    find(adjR2==max(adjR2),1);
    
    error = (RSS-RSSmin)/MSE;
    
    VT1best = VT1(kbest);
%     VEO2best=VEO2(kbest); %ADDED BY RG ON 8/12/2014
    
%     figure(20); % RSS evolution
%     clf;
%     plot(RSS); % RSS
%     axis([1 length((iVT1_first:iVT1_last)) 0 3]);
%     grid;   
%     hold on
%     plot(kbest, RSSmin, 'r+', 'MarkerSize', 36);
%     
    
    
    
%     figure(20)
%     clf;
%     [AX,H1,H2] = plotyy(1:length(VO2),RSS,1:length(VO2),MSE);
%     set(H1,'Color','r')
%     set(H2,'Color','b')
%     set(AX(1),'XLim',[1 220],'YLim',[0 3])
%     set(AX(2),'XLim',[1 220],'YLim',[0 0.01])
%     grid;   
%     hold on
%     plot(kbest, RSSmin, 'r+', 'MarkerSize', 36);
%     grid;   
%     hold on
%     plot(kbestMSE, MSEmin, 'b+', 'MarkerSize', 36);
    
    
    
    
    yVT1best = VCO2atVT1(kbest);
    VT1m1=m1(kbest);
    VT1m2=m2(kbest);
    VT1b1=b1(kbest);
    VT1b2=b2(kbest);
    yfitbestlow = m1(kbest).*VO2 + b1(kbest);
    yfitbesthigh = m2(kbest).*VO2 + b2(kbest);

  
    
    %axes(handles.axes_vslope); % VCO2 vs VO2 (V-Slopes plot)
%     cla(handles.axes_vslope);
%     plot(handles.axes_vslope,VO2, VCO2, 'b.'); % VCO2 vs VCO2 
%     axis(handles.axes_vslope,[0 VO2max+0.5 0 VO2max+0.5]);
%     xlabel(handles.axes_vslope,'VO2 (l/min)');
%     ylabel(handles.axes_vslope,'VCO2 (l/min)');
%     grid(handles.axes_vslope,'on');   
%     hold(handles.axes_vslope,'on')
%     plot(handles.axes_vslope,VO2, yfitbestlow, 'k', 'Linewidth', 1.0);
%     plot(handles.axes_vslope,VO2, yfitbesthigh, 'k', 'Linewidth', 1.0);
%     plot(handles.axes_vslope,VT1best, yVT1best, 'r+', 'MarkerSize', 36);
%     title(handles.axes_vslope,'VT1 estimation')
 
%     figure(23); % error plot: MSE vs VO2
%     clf;
%     plot(VO2(iVT1_first:iVT1_last),error(iVT1_first:iVT1_last), 'ko')
%     axis([0 max(VO2)  0 100]);
%     title(sprintf('V-Slopes Error Plot: data from Subject File %s', current_file));
%     xlabel('VO2 (L/min)');
%     ylabel('normalized error [(RSS-RSSmin)/MSE)]');
%     grid;   
%     hold on
%     %plot(VT1best, MSEbest, 'r+', 'MarkerSize', 24);
    
    
    

    

end

