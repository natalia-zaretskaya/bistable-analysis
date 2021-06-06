function out = analyzeReplay(data, replayTC, perceptChange, ifi, plotFlag, varargin)
% this function will analyze the following parameters of replay
% - mean/median reaction time (only for correct responses)
% - mean/median transparency level of each percept and mixed percept (only for correct responses)
% - nErrors (wrong button pressed)
% - nMissed (no button pressed)
%
% INPUTS:
% data -  shoudl be an n x 3 matrix with button presses. The columns are:
%   tonset
%   toffset
%   precept identity
%
% replaytc      - replay time courses from the experiment (stimulus.replay.)
% perceptChange - times at which a change occured from the experiment
% ifi           - interflip interval of the experiment
% plotFlag
% varargin{1} - transparency to percept mapping, e.g. [97 99]
%   first value corresponds to trasparency -1, second to 1
%
% If cell arrays are passed (e.g. 1cell = 1trial), each cell should contain
% matricies in the above format.
% If cell arrays are passed the parameters will also be computed per-cell
%
% Natalia Zaretskaya 2015.02.03
%
% see also: analyzeBistableKeys


if ~isempty(varargin)
  transparencyToKeyMapping = varargin{1};
else
  transparencyToKeyMapping = [];
end

if iscell(data) % also compute per-trial peramteres
    
    if plotFlag
        figure
        c = 1;
        for t = 1:length(data)
            uniquePercepts = unique(data{t}(:,3));
            uniquePercepts = uniquePercepts(uniquePercepts~=0);
            
            if isempty(transparencyToKeyMapping)
                transparencyToKeyMapping = [max(uniquePercepts(uniquePercepts>0)) min(uniquePercepts(uniquePercepts>0))];
            end
            
            percepts = data{t}(:,3);
            percepts(percepts==transparencyToKeyMapping(1)) = -1;
            percepts(percepts==transparencyToKeyMapping(2)) = 1;
            
            subplot(length(data), 1, c)
            plot((0:length(replayTC{t})-1).*ifi, replayTC{t}); hold on
            plot(perceptChange{t}(:,1), perceptChange{t}(:,2), 'g')
            plot(data{t}(:,1), percepts, 'r');
            c = c+1;
        end
        legend({'tc', 'expected press', 'actual press'})
        title('replay')
    end
    
    
    for t = 1:length(data)
        if ~isempty(data{t})
            trial = analyzeReplay(data{t}, replayTC{t}, perceptChange{t}, ifi, plotFlag); % cool! this function calls itself
            out.trial(t) = trial;
        else
            warning(sprintf('trial %u does not contain any data \n', t))
        end
    end
    
    
%     % combine data from all trials
       data = data(~cellfun('isempty', data));
       data = cell2mat(data(:));
       replayTC = cell2mat(replayTC);
       perceptChange = cell2mat(perceptChange(:));
end

% initialize variables
out.meanRT = [];
out.medianRT = [];
out.RT = [];
out.meanAlpha = []; % total viewing time of each percept
out.medianAlpha = [];
out.alpha = [];
out.nCorrect = [];
out.nWrong = [];
out.nMissed = [];
out.nPhysicalPercepts = [];

uniquePercepts = unique(data(:,3));


if ~isempty(varargin)
  transparencyToKeyMapping = varargin{1};
else
  transparencyToKeyMapping = [];
end

if isempty(transparencyToKeyMapping)
    transparencyToKeyMapping = [max(uniquePercepts(uniquePercepts>0)) min(uniquePercepts(uniquePercepts>0))];
end

% percept-to-tranparency mapping for now:
% TODO: find a universal approach for 3 percepts
allTransparencies = round(perceptChange(:,2));

allTransparencies(allTransparencies==-1) = transparencyToKeyMapping(1);
allTransparencies(allTransparencies==1) = transparencyToKeyMapping(2);

perceptChange(:,2) = allTransparencies;
rtInterval = [-0.2 2]; % in seconds
% forward analysis: stimulus -> response

out.perceptOrder = round(uniquePercepts);
for p = 1:length(uniquePercepts)
    
    
    thisPhysicalPerceptIndices = allTransparencies==uniquePercepts(p);
    thisPhysicalPerceptOnsets = perceptChange(thisPhysicalPerceptIndices,1);
    
    hits = 0;
    falseAlarms = 0;
    misses = 0;
    rt = [];
    alpha = [];
    
    for i = 2:length(thisPhysicalPerceptOnsets)-1 % 2:end-1 because the first is a mixed, and the last is a mixed
        
        % check whether there is any button press in this time interval
        response = intersect(find(data(:,1) > thisPhysicalPerceptOnsets(i)+rtInterval(1)), ...
            find(data(:,1) < thisPhysicalPerceptOnsets(i)+rtInterval(1)+rtInterval(2)));
        
        % compare response with percept
        if ~isempty(response)
            responseTime = data(response,1);
            responseId = data(response,3);
            
            correctResponse = find(responseId == uniquePercepts(p), 1);
            
            if isempty(correctResponse)
                falseAlarms = falseAlarms+1;
            else
                hits = hits+1;
                rt(end+1) = responseTime(correctResponse(1)) - thisPhysicalPerceptOnsets(i);
                
                ithTransparency = findNearestFrame(responseTime(correctResponse(1)), (0:1:length(replayTC)-1).*ifi);
                alpha(end+1) = replayTC(ithTransparency);
            end
            
        else % miss
            misses = misses+1;
            
        end
        
        
    end
    
    if ~isempty(rt)
      out.meanRT(p) = mean(rt);
      out.medianRT(p) = median(rt);
      out.RT{p} = rt;
    else
      out.meanRT(p) = NaN;
      out.medianRT(p) = NaN;
      out.RT{p} = NaN;
    end
    if ~isempty(alpha)
      out.meanAlpha(p) = mean(alpha);
      out.medianAlpha(p) = median(alpha);
      out.alpha{p} = alpha;
    else
      out.meanAlpha(p) = NaN;
      out.medianAlpha(p) = NaN;
      out.alpha{p} = NaN;
    end
    out.nCorrect(p) = hits;
    out.nWrong(p) = falseAlarms;
    out.nMissed(p) = misses;
    out.nPhysicalPercepts(p) = numel(thisPhysicalPerceptOnsets(2:end-1));
    
end % p


end


function i = findNearestFrame(thisTime, allTimes)

nearest = [];
i1 = find(allTimes<=thisTime, 1, 'last' );
if ~isempty(i1)
    nearest = [nearest; i1];
end

i2 = find(allTimes>=thisTime, 1 );
if ~isempty(i2)
    nearest = [nearest; i2];
end


distance = abs(allTimes(nearest)-thisTime);

i = nearest(distance==min(distance));

if ~isempty(i)
    i = i(end);
else
    
    error('could not find nearest frame; chech your timing')
    
end

end