
tVar = {};
for k = 1:numel(EEG.event)
    tVar{k} = EEG.event(k).type;
end



EventType ={'111', '112', '113', '114', '121', '122', '123', '124', '131', '132', '133', '134'};
idx = {};
for ThisEvent = 1:length(EventType)
idx{ThisEvent} = strfind(tVar,EventType{ThisEvent});
Tally = sum(cell2mat(idx{ThisEvent}));
disp([EventType{ThisEvent} ' was noted this many times:' num2str(Tally)]);
end

