function [cleanedPercepts] = cleanBistableKeys(originalPercepts, c, varargin)
%
% additional preselection of bistable responses
% (use e.g. after analyzeBistableKeys.m)
%
% INPUTS:
% ========
% PERCEPTS - matrix with onsets, offsets and ids of each percept in
% chronological order
%
% c - criteria sturcture with fields:
%   .minimalDuration - minimal allowed duration in time units
%   .maximalDuration - maximal allowed duration in time Units
%   .previousDifferent - only use percept onsets preceeded by a different
%   .previousMinimalDuration
%   .previousMaximalDuration
%   .nextDifferent
%   .nextMinimalDuration
%   .nextMaximalDuration
%
%
% OPTIONS:
% ========
%
% plotFlag - 0/1 whether to plot the data
%
%
% OUTPUTS:
% ========
% cleanPercept - same as percept, but containing only those originalPercepts that
% fullfill the c


% defaults:
if isempty(c)
    c.minimalDuration         = 0;
    c.maximalDuration         = Inf;
    c.previousDifferent       = 0;
    c.previousMinimalDuration = 0;
    c.previousMaximalDuration = Inf;
    c.nextDifferent           = 0;
    c.nextMinimalDuration     = 0;
    c.nextMaximalDuration     = Inf;
end
plotFlag = 0;

% check optional arguments
for i = 1:length(varargin)
    if strcmp(varargin{i}, 'plotFlag')
        plotFlag = varargin{i+1};
    end
    
end

cleanedPercepts = cell(size(originalPercepts));
% for each trial
for trl = 1:max(size(originalPercepts))
    
    d = originalPercepts{trl}(:,2)-originalPercepts{trl}(:,1);
    id = originalPercepts{trl}(:,3);
    
    for p = 1:length(id)
        
        % check each condition:
        if d(p) < c.minimalDuration; criteriaSatisfied = 0;  
        elseif d(p) > c.maximalDuration; criteriaSatisfied = 0;
            
        elseif p ~=1 && c.previousDifferent && id(p-1)==id(p); criteriaSatisfied = 0;
        elseif p ~=1 && d(p-1) < c.previousMinimalDuration; criteriaSatisfied = 0;
        elseif p ~=1 && d(p-1) > c.previousMaximalDuration; criteriaSatisfied = 0;
            
        elseif p ==1 && c.previousMinimalDuration~=0; criteriaSatisfied = 0;
        elseif p ==1 && c.previousMaximalDuration~=Inf; criteriaSatisfied = 0;
            
        elseif p ~=length(id) && c.nextDifferent && id(p+1)==id(p); criteriaSatisfied = 0;
        elseif p ~=length(id) && d(p+1) < c.nextMinimalDuration; criteriaSatisfied = 0;
        elseif p ~=length(id) && d(p+1) > c.nextMaximalDuration; criteriaSatisfied = 0;
            
        elseif p ==length(id) && c.nextMinimalDuration~=0;   criteriaSatisfied = 0;
        elseif p ==length(id) && c.nextMaximalDuration~=Inf; criteriaSatisfied = 0;
            
        else
            criteriaSatisfied = 1;            
        end
        
        
        if criteriaSatisfied
            cleanedPercepts{trl} = [cleanedPercepts{trl}; originalPercepts{trl}(p,:)];
        end
        
        
    end
    
end % for trl

% stats
for trl = 1:length(cleanedPercepts)
    fprintf('CLEANKEYSTAT trial %u \n', trl)
    
    if ~isempty(cleanedPercepts{trl})
        perceptTypes = unique(cleanedPercepts{trl}(:,3));
        
        for p = 1:length(perceptTypes)
            
            idx = cleanedPercepts{trl}(:,3) == perceptTypes(p);

            fprintf('CLEANKEYSTAT percept %u key %u \n', p, perceptTypes(p));
            fprintf('CLEANKEYSTAT median duration: %f \n', ...
                median(cleanedPercepts{trl}(idx,2)-cleanedPercepts{trl}(idx,1)));
            fprintf('CLEANKEYSTAT mean   duration: %f \n', ...
                mean(cleanedPercepts{trl}(idx,2)-cleanedPercepts{trl}(idx,1)));
        end
        
    end
    
end

if plotFlag
    
    figure;
    
    for trl = 1:length(cleanedPercepts)
        
        subplot(length(cleanedPercepts), 1, trl);
        
        perceptTypes = unique(originalPercepts{trl}(:,3));
        
        for p = 1:length(perceptTypes)
            
            % plot original data
            idx = find(originalPercepts{trl}(:,3) == perceptTypes(p));
            for i = 1:length(idx)
                plot([originalPercepts{trl}(idx(i),1) originalPercepts{trl}(idx(i),2)], [p p], 'LineWidth', 3, 'Color', 'b')
                hold on
            end
            ylim([0 length(perceptTypes)+1]);
            
            % plot clean data
            if ~isempty(cleanedPercepts{trl})
                
                idx = find(cleanedPercepts{trl}(:,3) == perceptTypes(p));
                
                for i = 1:length(idx)
                    plot([cleanedPercepts{trl}(idx(i),1) cleanedPercepts{trl}(idx(i),2)], [p p]+0.25, 'LineWidth', 3, 'Color', 'r')
                end
                
                hold on
                ylim([0 length(perceptTypes)+1]);
                
            end
            
        end
        
    end
    
end


end


