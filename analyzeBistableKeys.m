function [sortedPercepts, allPercepts] = analyzeBistableKeys(exp, key, varargin)
%
% analyze a continuous bistability period
%
% INPUTS:
% ========
%
% EXP - structure with fields .trialStartTime; .trialEndTime of one or multipe periods of
% stimulus on the screen
%
% KEY - structure with fields .idUp .idDown .timeUp .timeDown containing id of each up
% and down key episode
%
% OPTIONAL INPUTS:
% ========
% s - sturcture with fields
%   .responseType - instruction given to the subject on how to respond
%   'press-hold'
%   ------------- time
%    d--u d---u   perpcept
%
%   'press-release'
%   ------------- time
%    du     du    perpcept
%
%   .rejectLastPercept - this deals with the fact that the last percept is interrupted at trial end
%   ------------- time
%      d--u   d-- perpcept
%
%   1: regect
%   ------------- time
%      d--u       perpcept
%
%   0: perceptEnd = trialEndTime
%   ------------- time
%      d--u   d-u perpcept
%
%   .ignoreSameKeyPressedAgain = 0; what to do if same key is pressed reprearedly
% 	.perceptKeys = [37 79; 39 80]; % [nPercepts nKeys] only these keys will be considered to infer percepts
%       -  beware of differences on windows and mac!
%
% plotFlag - 0/1 whether to plot the data
%
%
% OUTPUTS:
% ========
% sortedPercepts - structure of size n_percepts, nTrials with fields
%   .keyId % key id of that sortedPercepts; 0 codes for mixed
%   .onset
%   .duration
%
% allPercepts - matrix with size n percepts x 3
%   colums: 1-percept start time 2-percept end time 3-percept key
%   (format useful for e.g. eyetracker data analysis)
%
% call the function from your scripts e.g. like this:
%   analyzeBistableKeys(log.exp, log.key, ?settings?, s, ?plotFlag?, 1);
%
% or using defaults below:
%   analyzeBistableKeys(log.exp, log.key);
%

% defaults:
s.responseType              = 'press-hold';
s.rejectLastPercept         = 0;
s.ignoreSameKeyPressedAgain = 1; % only applicable to press-release
s.perceptKeys               = [114 80 54 97 37 5; 115 79 55 99 39 6]; % beware of differences on windows and mac!
s.debug                     = 0;
s.excludeOverlaps           = 1;

plotFlag                    = 1;

% check optional arguments
for i = 1:length(varargin)
    
    if strcmp(varargin{i}, 'plotFlag')
        plotFlag = varargin{i+1};
    end
    
    if strcmp(varargin{i}, 'settings')
        s = varargin{i+1};
        if ~isfield(s, 'printStats') % make sure this setting is defined by the user
            error('please specify s.printStats setting (0,1) \n')
        end
    end
    
end


% get rid of keys that are not subject keys (e.g. scanner trigger):

isSubjectKeyDown = ismember(key.idDown, s.perceptKeys(:));
isSubjectKeyUp = ismember(key.idUp, s.perceptKeys(:));

key.idDown = key.idDown(isSubjectKeyDown);
key.timeDown = key.timeDown(isSubjectKeyDown);

key.idUp = key.idUp(isSubjectKeyUp);
key.timeUp = key.timeUp(isSubjectKeyUp);


% fix the bug: two near-identical keydown times with the same id (probably a joystick problem)
tmp = find(abs(diff(key.timeDown))<=0.001);

if key.idDown(tmp)==key.idDown(tmp+1)
    warning('difference between two neighbouring keydowns is smaller than 1ms')
    %             pause
    key.idDown(tmp)=[];
    key.timeDown(tmp)=[];
end


tmp = find(abs(diff(key.timeUp))<=0.001);

if key.idUp(tmp)==key.idUp(tmp+1)
    warning('difference between two neighbouring keyUps is smaller than 1ms')
    %             pause
    key.idUp(tmp)=[];
    key.timeUp(tmp)=[];
end

% determine the number of unique keys (usedKeys)
usedKeys = unique([key.idUp; key.idDown]);

% sort keys chronologically
[sorted, order] = sort(key.timeUp);
key.timeUp = sorted;
key.idUp = key.idUp(order);

[sorted, order] = sort(key.timeDown);
key.timeDown = sorted;
key.idDown = key.idDown(order);

sortedPercepts = struct([]);
allPercepts = cell(1,1);



% for each trial
for trl = 1:length(exp.trialStartTime)
    
    % find key events within this trial
    tmp = intersect( find(key.timeDown>exp.trialStartTime(trl)),...
        find(key.timeDown<exp.trialEndTime(trl)) );
    
    trialKeyTimeDown = key.timeDown(tmp);
    trialKeyIdDown = key.idDown(tmp);
    
    tmp = intersect( find(key.timeUp>exp.trialStartTime(trl)),...
        find(key.timeUp<exp.trialEndTime(trl)) );
    
    trialKeyIdUp = key.idUp(tmp);
    trialKeyTimeUp = key.timeUp(tmp);
    
    
    
    if ~isempty(trialKeyTimeDown) % if there is at least one percept in this trial
                
        if strcmp(s.responseType, 'press-release')
            
            if ~s.rejectLastPercept % add an artificial key press at the end of the trial
                trialKeyTimeDown = [trialKeyTimeDown; exp.trialEndTime(trl)];
                trialKeyIdDown = [trialKeyIdDown; trialKeyIdDown(end)];
                trialKeyTimeUp = [trialKeyTimeUp; exp.trialEndTime(trl)];
                trialKeyIdUp = [trialKeyIdUp; trialKeyIdDown(end)];
            end
            
            if s.ignoreSameKeyPressedAgain
                sameKeyRepeats = find(diff(trialKeyIdDown)==0);
                trialKeyIdDown(sameKeyRepeats+1) = [];
                trialKeyTimeDown(sameKeyRepeats+1) = [];
            end
            
            trialKeyTimeUp = trialKeyTimeDown(2:end);
            trialKeyTimeDown = trialKeyTimeDown(1:end-1);
            trialKeyIdDown = trialKeyIdDown(1:end-1);
            trialKeyIdUp = trialKeyIdDown;
            
        end % if press-release
        
        
        cleanTrialKeyTimeDown = [];
        cleanTrialKeyTimeUp = [];
        cleanTrialKeyIdDown = [];
        cleanTrialKeyIdUp = [];
        
        % for each percept type
        for p = 1:size(s.perceptKeys,1)
            
            theKey = intersect(s.perceptKeys(p,:), usedKeys);
            perceptStartTime{p} = trialKeyTimeDown(trialKeyIdDown==theKey);
            perceptEndTime{p} = trialKeyTimeUp(trialKeyIdUp==theKey);
            
            
            % new: do the key-cleanup for each percept separately
            if ~isempty(perceptStartTime{p}) && ~isempty(perceptEndTime{p})
                if perceptStartTime{p}(1)>perceptEndTime{p}(1)
                    perceptEndTime{p}(1) = [];
                    fprintf('deleted first upKey \n');
                end
                
                if ~isempty(perceptStartTime{p}) && ~isempty(perceptEndTime{p}) % if still not empty
                    if perceptStartTime{p}(end)>perceptEndTime{p}(end)
                        if s.rejectLastPercept
                            perceptStartTime{p}(end) = [];
                            if s.debug
                                fprintf('deleteing the last percept at trial end \n')
                            end
                        else
                            perceptEndTime{p}(end+1) = exp.trialEndTime(trl);
                            if s.debug
                                fprintf('adding the last keuUp at trial end \n')
                            end
                        end
                    end
                else
                    perceptStartTime{p} = [];
                    perceptEndTime{p} = [];
                end % if still not empty
            elseif ~isempty(perceptStartTime{p}) && isempty(perceptEndTime{p}) % if start isnt't empty, but stop is
                
                if ~s.rejectLastPercept % if we do not reject the last percept
                    perceptEndTime{p}(end+1) = exp.trialEndTime(trl);
                else
                    perceptStartTime{p} = [];
                    perceptEndTime{p} = [];
                end
            else
                
                perceptStartTime{p} = [];
                perceptEndTime{p} = [];
                
            end % if there is at least one percept
            
        end % for each sortedPercepts type
        
        
        % new2: additionally exclude overlaping key presses
        % you have to cross-check one percept agains the other
        if s.excludeOverlaps
            % check percept 2
            for i = 1:length(perceptStartTime{1})
                overlapIdx = intersect( find(perceptStartTime{2}>=perceptStartTime{1}(i)), find(perceptEndTime{2}<=perceptEndTime{1}(i))  );
                if ~isempty(overlapIdx)
                    fprintf('=> removing overlapping percept \n')
                    perceptStartTime{2}(overlapIdx) = [];
                    perceptEndTime{2}(overlapIdx) = [];
                end
            end
            
            % check percept 1
            for i = 1:length(perceptStartTime{2})
                overlapIdx = intersect( find(perceptStartTime{1}>=perceptStartTime{2}(i)), find(perceptEndTime{1}<=perceptEndTime{2}(i))  );
                if ~isempty(overlapIdx)
                    fprintf('=> removing overlapping percept \n')
                    perceptStartTime{1}(overlapIdx) = [];
                    perceptEndTime{1}(overlapIdx) = [];
                end
            end
        end
        
        
        
        
        for p = 1:length(perceptStartTime)
            theKey = usedKeys(p);
            sortedPercepts(p,trl).keyId = theKey;
            sortedPercepts(p,trl).onset = perceptStartTime{p} - exp.trialStartTime(trl);
            % note for replay keep the original timing:
            %             sortedPercepts(p,trl).onset = perceptStartTime;
            sortedPercepts(p,trl).duration = perceptEndTime{p}(:)-perceptStartTime{p}(:);
            
            cleanTrialKeyTimeDown = [cleanTrialKeyTimeDown; perceptStartTime{p}(:)];
            cleanTrialKeyTimeUp = [cleanTrialKeyTimeUp; perceptEndTime{p}(:)];
            cleanTrialKeyIdDown = [cleanTrialKeyIdDown; perceptStartTime{p}(:).*0+theKey];
            cleanTrialKeyIdUp = [cleanTrialKeyIdUp; perceptEndTime{p}(:).*0+theKey];
        end
        
        % reassamble trialKey Matrices after cleaning
        trialKeyTimeDown = cleanTrialKeyTimeDown;
        trialKeyTimeUp  = cleanTrialKeyTimeUp;
        trialKeyIdDown = cleanTrialKeyIdDown;
        trialKeyIdUp = cleanTrialKeyIdUp;
        
        % re-sort
        [~, order] = sort(trialKeyTimeDown);
        trialKeyTimeDown = trialKeyTimeDown(order);
        trialKeyIdDown = trialKeyIdDown(order);
        %         [~, order] = sort(trialKeyTimeUp);
        trialKeyTimeUp = trialKeyTimeUp(order);
        %         trialKeyIdUp = trialKeyIdUp(order); % should be the same as trialKeyIdUp
        
        
        % infer mixed perception duration from times when nothing was pressed (valid only for press-hold)
        if strcmp(s.responseType, 'press-hold')
            
            
            if ~isempty(trialKeyTimeUp(1:end-1)) && ~isempty(trialKeyTimeDown(2:end))
                
                sortedPercepts(p+1,trl).keyId = 0;
                % note, for replay use a different timeline:
                % sortedPercepts(p+1,trl).onset = trialKeyTimeUp(1:end-1);% - exp.trialStartTime(trl);
                sortedPercepts(p+1,trl).onset = trialKeyTimeUp(1:end-1) - exp.trialStartTime(trl);
                sortedPercepts(p+1,trl).duration = trialKeyTimeDown(2:end) - trialKeyTimeUp(1:end-1);
                
                % add additional onsets and offsets of the mixed percepts
                trialKeyIdDown2 = [trialKeyIdDown; zeros(size(trialKeyTimeUp(1:end-1)))];
                trialKeyTimeDown2 = [trialKeyTimeDown; trialKeyTimeUp(1:end-1)];
                trialKeyTimeUp2 = [trialKeyTimeUp(:); trialKeyTimeDown(2:end)];
            else
                trialKeyIdDown2 = trialKeyIdDown;
                trialKeyTimeDown2 = trialKeyTimeDown;
                trialKeyTimeUp2 = trialKeyTimeUp;
            end
            
        elseif strcmp(s.responseType, 'press-release') % no mixed percepts can be determined
            if ~isempty(trialKeyTimeUp(1:end-1)) && ~isempty(trialKeyTimeDown(2:end))
                
                sortedPercepts(p+1,trl).keyId = 0;
                % note, for replay use a different timeline:
                % sortedPercepts(p+1,trl).onset = trialKeyTimeUp(1:end-1);% - exp.trialStartTime(trl);
                sortedPercepts(p+1,trl).onset = trialKeyTimeDown(2:end);
                sortedPercepts(p+1,trl).duration = trialKeyTimeDown(2:end).*0;
                
                % add additional onsets and offsets of the mixed percepts
                trialKeyIdDown2 = [trialKeyIdDown; zeros(size( trialKeyTimeDown(2:end)))];
                trialKeyTimeDown2 = [trialKeyTimeDown;  trialKeyTimeDown(2:end)];
                trialKeyTimeUp2 = [trialKeyTimeUp(:); trialKeyTimeDown(2:end)];
            else
                trialKeyIdDown2 = trialKeyIdDown;
                trialKeyTimeDown2 = trialKeyTimeDown;
                trialKeyTimeUp2 = trialKeyTimeDown;
            end
            
        end
        
        
        % note: for replay use this timeline:
        % allPercepts{trl} = [[trialKeyTimeDown2 trialKeyTimeUp2] trialKeyIdDown2];
        allPercepts{trl} = [[trialKeyTimeDown2 trialKeyTimeUp2] - exp.trialStartTime(trl) trialKeyIdDown2];
        
        % sort percepts one final time in chronological order
        [~, order] = sort(allPercepts{trl}(:,1));
        allPercepts{trl}=allPercepts{trl}(order,:);
        
        if s.separateReturnTransitions
            % sort percepts chronologically:
            nonMixedPerceptIndices = find(allPercepts{trl}(:,3)~=0);
            mixedPerceptIndices = find(allPercepts{trl}(:,3)==0);
            nonMixedPercepts = allPercepts{trl}(nonMixedPerceptIndices, :);
            [~, order] = sort (nonMixedPercepts(:,1));
            nonMixedPercepts = nonMixedPercepts(order,:);
            
            perceptChanged = diff(nonMixedPercepts(:,3));
            allPercepts{trl}(mixedPerceptIndices(perceptChanged==0),3) = -1; %return transitions
        end
        
        
    else
        warning('No percepts in trial %u \n', trl);
        % for each percept type
        for p = 1:size(s.perceptKeys,1)
            
            sortedPercepts(p,trl).keyId = []; % theKey
            sortedPercepts(p,trl).onset = [];
            sortedPercepts(p,trl).duration = [];
            
        end % for each sortedPercepts type
        allPercepts{trl} = [];
        
    end % if there is at least one percept in this trial
    
end % for trl


% stats
for trl = 1:size(sortedPercepts,2)
    fprintf('KEYSTAT trial %u \n', trl)
    for p = 1:size(sortedPercepts,1)
      if ~isempty(sortedPercepts(p,trl).duration)
        fprintf('KEYSTAT percept %u key %u \n', p, sortedPercepts(p,trl).keyId);
        fprintf('KEYSTAT median duration: %f \n', median(sortedPercepts(p,trl).duration));
        fprintf('KEYSTAT mean   duration: %f \n', mean(sortedPercepts(p,trl).duration));
      end
    end
end


allPerceptsExperiment = cell2mat(allPercepts(~cellfun('isempty', allPercepts))');
if plotFlag && ~isempty(allPerceptsExperiment)
    figure;
    perceptColors = 'rgbk';
    for trl = 1:length(exp.trialStartTime)
        subplot(length(exp.trialStartTime), 1, trl);
        for p = 1:size(sortedPercepts,1)
            for i = 1:length(sortedPercepts(p,trl).onset)
                plot([sortedPercepts(p,trl).onset(i) sortedPercepts(p,trl).onset(i)+sortedPercepts(p,trl).duration(i)], [1 1].*sortedPercepts(p,trl).keyId, perceptColors(p), 'LineWidth', 3)
                hold on
                % note: for replay use this
                %                 xlim([exp.trialStartTime(trl) exp.trialEndTime(trl)]); %-exp.trialStartTime(trl));
                xlim([exp.trialStartTime(trl) exp.trialEndTime(trl)]-exp.trialStartTime(trl));%;
                ylim([0 200+1]);
                ylabel('sortedPercepts')
            end
        end
    end
    xlabel('time')
    
    allPerceptKeysExperiment = unique(allPerceptsExperiment(:,3));
    figure;
    for i = 1:length(allPerceptKeysExperiment);
        subplot(1, length(allPerceptKeysExperiment), i );
        hist(diff(allPerceptsExperiment(allPerceptsExperiment(:,3)==allPerceptKeysExperiment(i), [1 2]), 1, 2), 0:10 );
        xlim([0 10])
        xlabel('seconds');
        ylabel('frequency');
        title(sprintf('percept of key %u', allPerceptKeysExperiment(i)) );
    end
end

end

