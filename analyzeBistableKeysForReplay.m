function [sortedPercepts, allPercepts] = analyzeBistableKeysForReplay(exp, key, varargin)
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
s.perceptKeys               = [5 114; 6 115]; % beware of differences on windows and mac!
s.debug                     = 0;

plotFlag                    = 1;

% check optional arguments
for i = 1:length(varargin)
    
    if strcmp(varargin{i}, 'plotFlag')
        plotFlag = varargin{i+1};
    end
    
    if strcmp(varargin{i}, 'settings')
        s = varargin{i+1};
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

tmp = find(abs(diff(key.timeDown))<0.001);

if key.idDown(tmp)==key.idDown(tmp+1)

    warning('difference between two neighbouring keydowns is smaller than 0.1ms')

    key.idDown(tmp)=[];

    key.timeDown(tmp)=[];

end

 

tmp = find(abs(diff(key.timeUp))<0.001);

if key.idUp(tmp)==key.idUp(tmp+1)

    warning('difference between two neighbouring keyUps is smaller than 0.1ms')

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

if s.debug
    if isempty(key.timeDown)
        fprintf('no key presses were recorded \n')
        return
    elseif isempty(key.timeUp)
        fprintf('no key releases were recorded \n')
    else
        
        fprintf('total size of presses in the experiment %u \n', length(key.timeDown))
        fprintf('total size of releases in the experiment %u \n', length(key.timeUp))
        
        fprintf('first press at %4.2f s \n', key.timeDown(1)-exp.trialStartTime(1))
        fprintf('first release at %4.2f s \n', key.timeUp(1)-exp.trialStartTime(1))
        
        fprintf('last press at %4.2f s \n', key.timeDown(end)-exp.trialStartTime(1))
        fprintf('last release at %4.2f s \n', key.timeUp(end)-exp.trialStartTime(1))
    end
end

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
        
        if s.debug
            fprintf('=> in trial %u: \n', trl)
            fprintf('total size of presses %u \n', length(trialKeyTimeDown))
            fprintf('total size of releases %u \n', length(trialKeyTimeUp))
            
            fprintf('first press at %4.2f s \n', trialKeyTimeDown(1)-exp.trialStartTime(1))
            fprintf('first release at %4.2f s \n', trialKeyTimeUp(1)-exp.trialStartTime(1))
            
            fprintf('last press at %4.2f s \n', trialKeyTimeDown(end)-exp.trialStartTime(1))
            fprintf('last release at %4.2f s \n', trialKeyTimeUp(end)-exp.trialStartTime(1))
        end
        
        % clean the key responses:
        if ~isempty(trialKeyTimeDown) && ~isempty(trialKeyTimeUp)
            if trialKeyTimeDown(1)>trialKeyTimeUp(1) ||...
                    trialKeyIdDown(1)~=trialKeyIdUp(1) % if keyup preceeds keydown or keyUp and keyDown are not of the same percept
                if s.debug
                    fprintf('deleteing the first percept at trial start \n')
                end
                trialKeyTimeUp(1) = [];
                trialKeyIdUp(1) = [];
            end
            
            if trialKeyTimeDown(end)>trialKeyTimeUp(end) ||...
                    trialKeyIdDown(end)~=trialKeyIdUp(end) % if last keydown comes after last keyup
                if s.rejectLastPercept % delete the last keydown
                    if s.debug
                        fprintf('deleteing the last percept at trial end \n')
                    end
                    trialKeyTimeDown(end) = [];
                    trialKeyIdDown(end) = [];
                else  % add an artificial keyup at trial end
                    if s.debug
                        fprintf('adding artifical last percept at trial end \n')
                    end
                    trialKeyTimeUp(end+1) = exp.trialEndTime(trl);
                    trialKeyIdUp(end+1) = trialKeyIdDown(end);
                end
            end
        end % id ~isempty
        
        
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
            
        end % f press-release
        
        
        % for each percept type
        for p = 1:size(s.perceptKeys,1)
            
            theKey = intersect(s.perceptKeys(p,:), usedKeys);
            
            perceptStartTime = trialKeyTimeDown(trialKeyIdDown==theKey);
            perceptEndTime = trialKeyTimeUp(trialKeyIdUp==theKey);
            
            if s.debug
                fprintf('percept of the key %u\n', theKey)
                fprintf('size of presses %u \n', length(perceptStartTime))
                fprintf('size of releases %u \n', length(perceptEndTime))
            end
            
            sortedPercepts(p,trl).keyId = theKey;
            sortedPercepts(p,trl).onset = perceptStartTime - exp.trialStartTime(trl);
            sortedPercepts(p,trl).duration = perceptEndTime-perceptStartTime;
            
        end % for each sortedPercepts type
        
        % infer mixed perception duration from times when nothing was pressed (valid only for press-hold)
        if strcmp(s.responseType, 'press-hold')
            if ~isempty(trialKeyTimeUp(1:end-1)) && ~isempty(trialKeyTimeDown(2:end))
                sortedPercepts(p+1,trl).keyId = 0;
                sortedPercepts(p+1,trl).onset = trialKeyTimeUp(1:end-1);% - exp.trialStartTime(trl);
                sortedPercepts(p+1,trl).duration = trialKeyTimeDown(2:end) - trialKeyTimeUp(1:end-1);
                
                % add additional onsets and offsets of the mixed percepts
                trialKeyIdDown2 = [trialKeyIdDown(:); zeros(size(trialKeyTimeUp(1:end-1)))];
                trialKeyTimeDown2 = [trialKeyTimeDown(:); trialKeyTimeUp(1:end-1)];
                trialKeyTimeUp2 = [trialKeyTimeUp(:); trialKeyTimeDown(2:end)];
            else
                trialKeyIdDown2 = trialKeyIdDown;
                trialKeyTimeDown2 = trialKeyTimeDown;
                trialKeyTimeUp2 = trialKeyTimeUp;
            end
        end
        
%         allPercepts{trl} = [[trialKeyTimeDown2 trialKeyTimeUp2]-exp.trialStartTime(trl) trialKeyIdDown2];
        allPercepts{trl} = [[trialKeyTimeDown2 trialKeyTimeUp2] trialKeyIdDown2];
                      


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
    
    
            % temporary solution
        if isequal(size(allPercepts{trl}), [1 2]);
            allPercepts{trl} = [];
        end
    
end % for trl


% stats
for trl = 1:size(sortedPercepts,2)
    fprintf('KEYSTAT trial %u \n', trl)
    for p = 1:size(sortedPercepts,1)
        fprintf('KEYSTAT percept %u key %u \n', p, sortedPercepts(p,trl).keyId);
        fprintf('KEYSTAT median duration: %f \n', median(sortedPercepts(p,trl).duration));
        fprintf('KEYSTAT mean   duration: %f \n', mean(sortedPercepts(p,trl).duration));
    end
end


allPerceptsExperiment = cell2mat(allPercepts(~cellfun('isempty', allPercepts))');
if plotFlag && ~isempty(allPerceptsExperiment)
    figure;
    for trl = 1:length(exp.trialStartTime)
        subplot(length(exp.trialStartTime), 1, trl);
        for p = 1:size(sortedPercepts,1)
            for i = 1:length(sortedPercepts(p,trl).onset)
                try
                    plot([sortedPercepts(p,trl).onset(i) sortedPercepts(p,trl).onset(i)+sortedPercepts(p,trl).duration(i)], [p p], cols(p), 'LineWidth', 3)
                    xlim([exp.trialStartTime(trl) exp.trialEndTime(trl)]-exp.trialStartTime(1));
                    ylim([0 length(sortedPercepts)+1]);
                    ylabel('sortedPercepts')
                    hold on
                catch
                end
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

