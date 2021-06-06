function [ exp, key , info ] = createSurrogateRivalry( data, returnOnNth, design, perceptKeys, nrepeats )
% Creates surrogate rivalry time course from real known distributions
% of  in the format readable by analyzeBistableKeys.m
%
% Inputs:
%
%   data - cell array 3x1 containgin vectors durations
%     data{1} - durations of percept 1
%     data{2} - durations of percept 2
%     data{3} - durations of mixed percepts
%
% returnOnNth - vector or scalar - each nth traisiton is a return transiton
%   (i.e. p1-mixed-p1); If vector, as sequence is randomized
%
% design - structure with fields
%    .trialDuration - in seconds
%    .trialSequence - vector with integer trial types (zero is treated as baseline)
%
% perceptKeys - 1x2 [keyID1 keyID2] integers - pretend that subject pressed
% these keys
%
% nrepeats - will repeat this number of times, with each iteration
%   delivering unique percept durations (i.e. drawing without repetitions)
%
% Natalia Zaretskaya 05.02.2015



for n = 1:nrepeats

    key(n).idDown = [];
    key(n).idUp = [];
    key(n).timeDown = [];
    key(n).timeUp = [];
    
    exp(n).trialStartTime = [];
    exp(n).trialEndTime = [];

    
    % start creating a sequence:
    currentTime = 0;
    
    for t = 1:length(design.trialDuration)
        
        exp(n).trialStartTime(t) = currentTime;
        
        if design.trialSequence(t)% in not a baseline period
            
            % randomize on every trial
            returnOnNth = returnOnNth(randperm(length(returnOnNth)));
            
            % continue iterating between percep and mix
            p = randperm(2); % start with a random percept
            c = 1; % iteration count
            
            while currentTime<=exp(n).trialStartTime(t)+design.trialDuration(t)
                
                % start with mixed,do not record button presses for mixed
                i = randi(length(data{3}));
                currentTime = currentTime+data{3}(i);
                
                % delete the chosen mixed percept
                data{3}(i) = [];
                
                % choose duration
                i = randi(length(data{p(1)}));
                
                % record current time
                key(n).idDown = [key(n).idDown; perceptKeys(p(1))];
                key(n).timeDown = [key(n).timeDown; currentTime];
                
                % add duration
                currentTime = currentTime+data{p(1)}(i);
                
                % record current time
                key(n).idUp = [key(n).idUp; perceptKeys(p(1))];
                key(n).timeUp = [key(n).timeUp; currentTime];
                
                % delete the chosen percept
                data{p(1)}(i) = [];
                
                % exchange percepts, but not on the 5th iteration
                if mod(c, returnOnNth(1))
                    p = [p(2) p(1)];
                else 
                    returnOnNth = [returnOnNth(2:end) returnOnNth(1)];
                    c = 0;
                end
                
                c = c+1;
            end
            
            
        else % if baseline
            
        end % if not baseline trial
        
        currentTime = exp(n).trialStartTime(t) + design.trialDuration(t); % recet current time
        
        exp(n).trialEndTime(t) = exp(n).trialStartTime(t) + design.trialDuration(t);
        
    end % trial
      
end % for n

% also save the information about the used dataset:
info.mean = [mean(data{1}) mean(data{2}) mean(data{3})];
info.median = [median(data{1}) median(data{2}) median(data{3})];
info.returnOnNth = returnOnNth;

end

