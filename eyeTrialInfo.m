function data = eyeTrialInfo(data, eyeTrials, trialSequence)
% adds percept-eye combination information to rivalry data
%
% INPUTS:
% data          - cell array, second output of analyzeBistableKeys.m
% eyeTrials     - 2x1 cell array were eyeTrials{1} contains trial ids with
% stimulus1 -> left eye, and stimulus2 ->right eye, and eyeTrials{2} trial
% ids where stimulus is swaped between the eyes
% trialSequence - array with trial ids as they occured in the experiment
%
% OUTPUTS:
% data          - same as input data, but with a forth column containing
% eye info for each percept (1 - stimulus1->left eye, 2 - stimulus1->right eye
%
% Natalia Zaretskaya 2016.07.12


tmp = cell2mat(data');
percepts = unique (tmp(:,3));
percepts = percepts(percepts>0); % real percepts (i.e. not mixtures)
if numel(percepts)==1
    warning('eye dominance could not be determined because only one percept was reported');
    return
else
    for t = 1:length(data)
        if ~isempty(data{t})
            
            eyeVector = zeros(size(data{t},1),1);
            
            if ismember(trialSequence(t), eyeTrials{1})
                
                thisPerceptIndices = data{t}(:,3)== percepts(1);
                eyeVector(thisPerceptIndices) = 1;
                thisPerceptIndices = data{t}(:,3)== percepts(2);
                eyeVector(thisPerceptIndices) = 2;
                
            elseif ismember(trialSequence(t), eyeTrials{2})
                thisPerceptIndices = data{t}(:,3)== percepts(1);
                eyeVector(thisPerceptIndices) = 2;
                thisPerceptIndices = data{t}(:,3)== percepts(2);
                eyeVector(thisPerceptIndices) = 1;
                
            else
                warning('undefined eye information in this trial')
            end
            
            data{t} = [data{t} eyeVector];
        end
    end % t
    
end



end