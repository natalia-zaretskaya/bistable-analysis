function out = analyzeRivalry(data, varargin)
% this function will analyze the following parameters of binocular rivalry
% - mean/median/skewness/kurtosis percept duration
% - mean/median/skewness/curtosis mixed duration
% - reversal rate
% - rate of return transitions
% - mean/median duration fo each eye
% - time to onset of rivalry relative to trial start (added 13.08.2015)
%
% The input shoudl be either a n x 3(or n x 4/5 if eye/color information is present)
% tarray with percept data. The columns are:
%   tonset
%   toffset
%   identity
%   (and optionally eye)
%   (and optionally color)
%
% If cell array is passed (e.g. 1cell = 1trial), each cell should contain
% matricies in the above format.
% If cell array is passed the parameters will also be computed per-cell
%
% Natalia Zaretskaya 2015.02.03
% update 2016.12.16 added std and iqr to analyze the variation of the data
% update 2017.0329 added n to count occurances of each percept
% update 2020.0224 added color information for the analysis
% update 2020.0503 octave compatibility
% see also: analyzeBistableKeys


thresholdDurations = false;
thresholdValue = 0;
suppressText = false;
plotFlag = 0;

for i = 1:length(varargin)
    if strcmp(varargin{i}, 'durationThreshold')
        thresholdDurations = true;
        thresholdValue = varargin{i+1};
    elseif strcmp(varargin{i}, 'suppressText')
        suppressText = varargin{i+1};
    elseif strcmp(varargin{i}, 'plotFlag')
        plotFlag = varargin{i+1};
    end

end

% initialize variables
out.meanPercept = [];
out.stdPercept = [];
out.medianPercept = [];
out.iqrPercept = [];
out.skewnessPercept = [];
out.kurtosisPercept = [];
out.medianBothPercepts = [];
out.iqrBothPercepts = [];
out.medianBothMixed = [];
out.iqrBothMixed = [];
out.timePercept = []; % total viewing time of each percept
out.nPercept = []; % number of occurances of each percept
out.meanEye = [];
out.medianEye = [];
out.timeEye = [];
out.meanColor = [];
out.medianColor = [];
out.timeColor = [];
out.nTransitions = [];
out.nReversals = [];
out.nReturns = [];
out.timeToOnset = [];

if iscell(data) % also compute per-trial peramteres
    for t = 1:length(data)
        if ~isempty(data{t})
            out.trial(t) = analyzeRivalry(data{t}, 'durationThreshold', thresholdValue, 'suppressText', suppressText); % cool! this function calls itself
        else
            warning('trial %u does not contain any data \n', t);
        end
    end

    % combine data from all trials
    data = data(~cellfun('isempty', data));
    data = cell2mat(data(:));

end

if isempty(data)
    warning('data is empty; nothing to analyze')
    out.meanPercept = NaN;
    out.stdPercept = NaN;
    out.medianPercept = NaN;
    out.iqrPercept = NaN;
    out.skewnessPercept = NaN;
    out.kurtosisPercept = NaN;
    out.medianBothPercepts = NaN;
    out.iqrBothPercepts = NaN;
    out.medianBothMixed = NaN;
    out.iqrBothMixed = NaN;
    out.timePercept = NaN; % total viewing time of each percept
    out.nPercept = NaN; % number of occurances of each percept
    out.meanEye = NaN;
    out.medianEye = NaN;
    out.timeEye = NaN;
    out.meanColor = NaN;
    out.medianColor = NaN;
    out.timeColor = NaN;
    out.nTransitions = NaN;
    out.nReversals = NaN;
    out.nReturns = NaN;
    out.timeToOnset = NaN;
    out.trial = [];
    return
end


uniquePercepts = unique(data(:,3));
out.uniquePercepts = uniquePercepts;

out.timeToOnset = data(1,1); % time to rivalry onset = time until the first button press

% --- percept-sepcific data --- @
for p = 1:length(uniquePercepts)
    thisPerceptIndices = data(:,3)==uniquePercepts(p);
    thisPerceptDurations = data(thisPerceptIndices, 2) - data(thisPerceptIndices, 1);

    if thresholdDurations
        thisPerceptDurations = thisPerceptDurations(thisPerceptDurations>thresholdValue);
    end

if not(isempty(thisPerceptDurations))
    out.meanPercept(p) = mean(thisPerceptDurations);
    out.stdPercept(p) = std(thisPerceptDurations);
    out.medianPercept(p) = median(thisPerceptDurations);
    
    if length(thisPerceptDurations) == 1
      out.iqrPercept(p) = 0;
    else
      out.iqrPercept(p) = iqr(thisPerceptDurations);
    end
    out.skewnessPercept(p) = skewness(thisPerceptDurations);
    out.kurtosisPercept(p) = kurtosis(thisPerceptDurations);
    out.timePercept(p) = sum(thisPerceptDurations);
    out.nPercept(p) = numel(thisPerceptDurations);
end

end

% NOTE: median([a b]) ~= median([median(a) median(b)]) ~= mean([median(a) median(b)])
%   and not always: mean([a b]) ~= mean([mean(a) mean(b)])
thisPerceptIndices = data(:,3)>0;
thisPerceptDurations = data(thisPerceptIndices, 2) - data(thisPerceptIndices, 1);
if thresholdDurations
    thisPerceptDurations = thisPerceptDurations(thisPerceptDurations>thresholdValue);
end
out.medianBothPercepts = median(thisPerceptDurations);

if numel(thisPerceptDurations)>1
  out.iqrBothPercepts = iqr(thisPerceptDurations);
else
out.iqrBothPercepts = NaN;
end
out.meanBothPercepts = mean(thisPerceptDurations);
out.stdBothPercepts = std(thisPerceptDurations);

thisPerceptIndices = data(:,3)<=0;
thisPerceptDurations = data(thisPerceptIndices, 2) - data(thisPerceptIndices, 1);
if thresholdDurations
    thisPerceptDurations = thisPerceptDurations(thisPerceptDurations>thresholdValue);
end

if ~isempty(thisPerceptDurations)
  out.medianBothMixed = median(thisPerceptDurations);
else
  out.medianBothMixed =  NaN;
end

if numel(thisPerceptDurations)>1
out.iqrBothMixed = iqr(thisPerceptDurations);
else
out.iqrBothMixed = NaN;
end
out.meanBothMixed = mean(thisPerceptDurations);
out.stdBothMixed = std(thisPerceptDurations);

% --- Transition count--- @
out.nTransitions = numel(find(data(:,3)<=0)); % reversals (0) + returns (-1)
out.nReversals = numel(find(data(:,3)==0));
out.nReturns = numel(find(data(:,3)==-1));


% --- eye-specific ifnormation --- @
if size(data,2) >= 4

    thisEyeIndices1 = ((data(:,4)==1) + (data(:,3)==min(uniquePercepts(uniquePercepts>0))))==2; % image of key 5 -> left eye
    thisEyeIndices2 = ((data(:,4)==2) + (data(:,3)==max(uniquePercepts(uniquePercepts>0))))==2; % image of key 6 -> left eye
    thisEyeIndices = [thisEyeIndices1 + thisEyeIndices2] ==1;

    thisEyeDurations = data(thisEyeIndices, 2) - data(thisEyeIndices, 1);
    out.meanEye(1) = mean(thisEyeDurations);
    out.medianEye(1) = median(thisEyeDurations);
    out.timeEye(1) = sum(thisEyeDurations);

    thisEyeIndices1 = ((data(:,4)==2) + (data(:,3)==min(uniquePercepts(uniquePercepts>0))))==2; % image of key 5 -> right eye
    thisEyeIndices2 = ((data(:,4)==1) + (data(:,3)==max(uniquePercepts(uniquePercepts>0))))==2; % image of key 6 -> right eye
    thisEyeIndices = [thisEyeIndices1 + thisEyeIndices2] ==1;

    thisEyeDurations = data(thisEyeIndices, 2) - data(thisEyeIndices, 1);
    out.meanEye(2) = mean(thisEyeDurations);
    out.medianEye(2) = median(thisEyeDurations);
    out.timeEye(2) = sum(thisEyeDurations);

end


% --- color-specific ifnormation --- @
if size(data,2) == 5

% % %     for c = 1:2
% % %         thisColorIndices = data(:,5)==c;
% % %         thisColorDurations = data(thisColorIndices, 2) - data(thisColorIndices, 1);
% % %         out.meanColor(c) = mean(thisColorDurations);
% % %         out.medianColor(c) = median(thisColorDurations);
% % %         out.timeColor(c) = sum(thisColorDurations);
% % %         clear thisColorIndices thisColorDurations
% % %     end

if size(data,2) >= 4

    thisColorIndices1 = ((data(:,5)==1) + (data(:,3)==min(uniquePercepts(uniquePercepts>0))))==2; % image of key 5 -> red color
    thisColorIndices2 = ((data(:,5)==2) + (data(:,3)==max(uniquePercepts(uniquePercepts>0))))==2; % image of key 6 -> red color
    thisColorIndices = [thisColorIndices1 + thisColorIndices2] ==1;

    thisColorDurations = data(thisColorIndices, 2) - data(thisColorIndices, 1);
    out.meanColor(1) = mean(thisColorDurations);
    out.medianColor(1) = median(thisColorDurations);
    out.timeColor(1) = sum(thisColorDurations);

    thisColorIndices1 = ((data(:,5)==2) + (data(:,3)==min(uniquePercepts(uniquePercepts>0))))==2; % image of key 5 -> green color
    thisColorIndices2 = ((data(:,5)==1) + (data(:,3)==max(uniquePercepts(uniquePercepts>0))))==2;  % image of key 6 -> green color
    thisColorIndices = [thisColorIndices1 + thisColorIndices2] ==1;

    thisColorDurations = data(thisColorIndices, 2) - data(thisColorIndices, 1);
    out.meanColor(2) = mean(thisColorDurations);
    out.medianColor(2) = median(thisColorDurations);
    out.timeColor(2) = sum(thisColorDurations);

end




end


if ~suppressText
    fprintf('thresholding percept durations by %f \n', thresholdValue)
    if size(data,2) == 4
        fprintf('eye information is present! will also analyze eye dominance... \n');


        fprintf('eye dominance info:\n');
        fprintf('EYESTAT total viewing time LE: %f, RE: %f\n', out.timeEye(1), out.timeEye(2));
        fprintf('EYESTAT median dominance   LE: %f, RE: %f\n', out.medianEye(1), out.medianEye(2))
        fprintf('EYESTAT mean dominance     LE: %f, RE: %f\n', out.meanEye(1), out.meanEye(2))
    end
end

if plotFlag
    figure;
    subplot(3,3,1)
    bar([1 2], out.meanPercept(out.uniquePercepts>0)); title('mean percept')
    subplot(3,3,2)
    bar([1 2], out.medianPercept(out.uniquePercepts>0)); title('median percept')
    subplot(3,3,3)
    bar([1 2], out.timePercept(out.uniquePercepts>0)); title('total time percept')

    if size(data,2) >= 4
    subplot(3,3,4)
    bar([1 2], out.meanEye); title('mean eye')
    subplot(3,3,5)
    bar([1 2], out.medianEye); title('median eye')
    subplot(3,3,6)
    bar([1 2], out.timeEye); title('total time eye')

    subplot(3,3,7)
    bar([1 2], out.meanColor); title('mean color')
    subplot(3,3,8)
    bar([1 2], out.medianColor); title('median color')
    subplot(3,3,9)
    bar([1 2], out.timeColor); title('total time color')
    end
end

end
