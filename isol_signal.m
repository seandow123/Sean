function isol = isol_signal(signal)
% isolate the signal for analysis

len = length(signal);
filteredSignal = medfilt1(abs(signal), 1000);

quietParts = filteredSignal < 1;

startingBlockIndexes = find(diff(quietParts) < 0);
endingBlockIndexes   = find(diff(quietParts) > 0);

if (isempty(startingBlockIndexes))
    startingBlockIndexes=10;
end

if (isempty(endingBlockIndexes))
    endingBlockIndexes =len-10;
end
startingBlockIndexes;
endingBlockIndexes;

isol = signal(startingBlockIndexes(1):endingBlockIndexes(1));
% Test plot
%     close all
%     figure; hold on
%     plot(filteredSignal,'-k');
%     plot(100*quietParts,'-r');
%     plot(isol,'--k');
end