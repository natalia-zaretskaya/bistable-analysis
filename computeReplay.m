function [ replay ] = computeReplay(rivalryFileName, replay, d,  plotFlag)

% computes the vector that indicates the transition stage between two
% images for every video frame depending on the trial type gives either
% smooth or abrupt function
% INPUTS:
% replayFileName - file name to be reaplayed
% replay - replay structure with replay-specific info
% d - desing structure for replay experiment
% ifi - interframe interval of the screen the replay is ran on (can be
%   different from the screen where rivalry was collected)
% plotFlag
%
% 2017.10.09:
% new - option to shuffle the individual percepts around, so that the the
% distribution is the same, but the subject can't remember the order

replay.rivalryFileName = rivalryFileName;

load(deblank(replay.rivalryFileName), 'log', 'design', 'ptb');
ifi = ptb.scrn.ifi;
[~, percepts] = analyzeBistableKeysForReplay(log.exp, log.key, 'plotFlag', plotFlag);

if isfield(d, 'shuffleFlag') 
    if d.shuffleFlag
        fprintf('shuffling percept durations for replay ... \n')
        % do not shuffle the order, just the duration and occasionally swap
        % percept ids with each other
        
        % transform on-off-id format into dur-id format
        percepts2 = [percepts(:,2)-percepts(:,1) percepts(:,2)];
        
        % shuffle durations
        nonMixedIdx = find(percepts(:,1)~=0);
        nonMixedIdxSchuffled = nonMixedIdx(randperm(length(nonMixedIdx)));
        percepts2(nonMixedIdx,1) = percepts2(nonMixedIdxSchuffled,1); % this shuffles the durations;
        
        % shuffle mixed durations
        mixedIdx = find(percepts(:,1)==0);
        mixedIdxSchuffled = mixedIdx(randperm(length(mixedIdx)));
        percepts2(mixedIdx,1) = percepts2(mixedIdxSchuffled,1); % this shuffles the durations;
        
        % transform back into x by 3 on-off-id format
        percepts3 = [ cumsum(percepts2(:,1)) - percepts2(1,1) cumsum(percepts2(:,1)) percepts2(:,2)];
       
        % randomly swap or not-swap the percept identities
        tmp = randperm(2);
        if tmp(1) ==1;
            uniquePercepts = unique(percepts(:,3));
            uniquePercepts(uniquePercepts==0) = [];
            idxPercept1 = percepts3(:,3) ==uniquePercepts(1);
            idxPercept2 = percepts3(:,3) ==uniquePercepts(2);
            
            percepts3(idxPercept1,3) =  uniquePercepts(2);
            percepts3(idxPercept2,3) =  uniquePercepts(1);
        end
        
        percepts = percepts3;
    end
end


replay.transparencyVector = cell(1, size(percepts,2));

if plotFlag; figure; hold on; end

for trl = 1:size(percepts,2)
    
    if design.trialSequence(trl)
        
        % create transparency timeline:
        timeValues = 0:ifi:design.trialDuration(trl);
        perceptionValuesAbrupt = zeros(2, length(timeValues));
        
        for ptype = 1:2
            
            % substract trial start time
            percepts(ptype, trl).onset = percepts(ptype, trl).onset-log.exp.trialStartTime(trl);
            
            for p = 1:length(percepts(ptype,trl).onset)
                startIndex = find(timeValues...
                    >percepts(ptype,trl).onset(p), 1, 'first');
                stopIndex = find(timeValues...
                    > (percepts(ptype,trl).onset(p)+percepts(ptype,trl).duration(p)), 1, 'first');
                
                % for abrupt replay:
                perceptionValuesAbrupt(ptype,startIndex:stopIndex)= ptype;
                
            end % percept
            
        end % percept type
        
        
        % abrupt replay
        perceptionValuesAbrupt = sum(perceptionValuesAbrupt);
        perceptionValuesAbrupt(perceptionValuesAbrupt==3) = 1.5; % overlapping percepts = mixed
        perceptionValuesAbrupt(perceptionValuesAbrupt==0) = 1.5; % real mixed
        perceptionValuesAbrupt = (perceptionValuesAbrupt-min(perceptionValuesAbrupt(:))); % scale to go from 0 to 1
        
        perceptionValuesAbrupt = perceptionValuesAbrupt.*2-1;
        
        % calculate smooth replay from abrupt:
        perceptChangeOnsets = [1 find(diff(perceptionValuesAbrupt)) length(perceptionValuesAbrupt)];
        perceptPeaks = zeros(length(perceptChangeOnsets)-1,2);
        
        for p = 1:length(perceptChangeOnsets)-1
            
            % get the center index of current percept
            centerIdx = round((perceptChangeOnsets(p+1)-perceptChangeOnsets(p))/2)+perceptChangeOnsets(p);
            perceptPeaks(p,:) = [centerIdx perceptionValuesAbrupt(centerIdx)];
            
            
        end
        
        % add first and last values to perceptPeaks for interpolation
        perceptPeaks = [1 0; perceptPeaks; length(timeValues) 0];
        
        
        [knownTimeIndices]= ismember( 1:length(timeValues), perceptPeaks(:,1));
        unknownTimeIndices = ~knownTimeIndices;
        
        perceptionValuesSmooth=zeros(size(timeValues));
        perceptionValuesSmooth(knownTimeIndices) = perceptPeaks(:,2);
        
        [ Vq ] = interp1(find(knownTimeIndices), perceptPeaks(:,2), find(unknownTimeIndices), 'pchip', 'extrap');
        perceptionValuesSmooth(unknownTimeIndices) = Vq;
        
        % --- depending on whether smooth of abrupt replay requried --- %
        
        if d.replayTransitionPattern(trl) == -1 % abrupt
            finalPerceptionValues = perceptionValuesAbrupt;
            
        elseif d.replayTransitionPattern(trl) ==1 % smooth
            
            finalPerceptionValues = perceptionValuesSmooth;
        else
            error('undefined temporal transition type (should be 0=abrupt or 1=smooth)')
        end
        
        
        % --- invert temporal order if replay backwards --- %
        
        if d.replayDirection(trl)==1
            replay.transparencyVector{trl} = finalPerceptionValues;
        elseif d.replayDirection(trl)==-1
            replay.transparencyVector{trl} = flip(finalPerceptionValues);
            fprintf('inverting time course ... \n')
        else
            error('undefined replay direction (should be 1 or -1)')
        end
        
        % save this for analysis:
        replay.perceptChangeOnsets{trl} = [(perceptChangeOnsets(1:end-1)'-1).*ifi perceptionValuesAbrupt(perceptChangeOnsets(1:end-1)+1)'];
        replay.perceptPeaks{trl} = [perceptPeaks(:,1)*ifi perceptPeaks(:,2)];
        
        if plotFlag
            subplot(size(percepts,2), 1, trl)
            plot(perceptionValuesAbrupt, 'b', 'LineWidth', 0.5);
            hold on
            plot(perceptionValuesSmooth, 'b', 'LineWidth', 0.5);
            % also plot the one that is actuallu used in this trial in thick
            plot(replay.transparencyVector{trl}, 'r', 'LineWidth', 2)
            plot(perceptPeaks(:,1), perceptPeaks(:,2), '*')
            plot(perceptChangeOnsets(1:end-1), perceptionValuesAbrupt(perceptChangeOnsets(1:end-1)+1), 'g*')
            title(sprintf('Trial %u', trl));
            xlabel('frame no')
            ylabel('percept')
        end % plotFlag
        
    else
        timeValues = 0:ifi:design.trialDuration(trl);
        replay.transparencyVector{trl} = zeros(size(timeValues));
    end % if non-baseline trial
    
    
    
end % trial
end


