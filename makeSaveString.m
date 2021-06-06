function saveString = makeSaveString(experimentType, runString)
    
if isempty(runString)
    saveString = 'tmp.mat';
else
    thedatestr = datestr(now); thedatestr(isspace(thedatestr)) = '_';  thedatestr(thedatestr==':') = '_';
    saveString = [experimentType, '_', runString, '_', thedatestr, '.mat'];
end


if ~isdir('log')
    mkdir('log')
end
saveString = fullfile('log', saveString);

fprintf('=>Saving data under %s\n', saveString);

end



