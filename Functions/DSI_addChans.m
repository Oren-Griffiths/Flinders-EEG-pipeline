
function EEG = DSI_addChans(EEG)
% takes the imported chan labels from default .edf file for DSI-24 and
% relabels them in a way that supports application of channel locations.

if ~isfield(EEG, 'chanlocs')
    disp('no channel labels imported. Import error');
    return
else % channel location made it along
    if numel(EEG.chanlocs) == 25
        % add channel labels to raw data structure.
        EEG.chanlocs(1).labels = 'P3';
        EEG.chanlocs(2).labels = 'C3';
        EEG.chanlocs(3).labels = 'F3';
        EEG.chanlocs(4).labels = 'Fz';        
        EEG.chanlocs(5).labels = 'F4';
        EEG.chanlocs(6).labels = 'C4';
        EEG.chanlocs(7).labels = 'P4';
        EEG.chanlocs(8).labels = 'Cz';
        EEG.chanlocs(9).labels = 'CM';
        EEG.chanlocs(10).labels = 'A1';
        EEG.chanlocs(11).labels = 'Fp1';
        EEG.chanlocs(12).labels = 'Fp2';
        EEG.chanlocs(13).labels = 'T3';
        EEG.chanlocs(14).labels = 'T5';
        EEG.chanlocs(15).labels = 'O1';
        EEG.chanlocs(16).labels = 'O2';
        EEG.chanlocs(17).labels = 'X3';
        EEG.chanlocs(18).labels = 'X2';
        EEG.chanlocs(19).labels = 'F7';
        EEG.chanlocs(20).labels = 'F8';
        EEG.chanlocs(21).labels = 'X1';
        EEG.chanlocs(22).labels = 'A2';
        EEG.chanlocs(23).labels = 'T6';
        EEG.chanlocs(24).labels = 'T4';
        EEG.chanlocs(25).labels = 'Trigger';
        
        % reorganize the structure so that scalp channels come first.
        ChannelReorder = ...
            [1 2 3 4 5 6 7 8 11 12 13 14 15 16 19 20 23 24 10 22 9 21 18 17 25 ];
        % switch to tables for manipulation and then back. 
        T = struct2table(EEG.chanlocs);
        Tnew = T(ChannelReorder,:);
        EEG.chanlocs = table2struct(Tnew);
        % and change the data too!
        EEG.data = EEG.data(ChannelReorder,:);
        % manipulation done.

        
    else % wrong number of channels. something has gone awry.
        disp('error in number of channels for DSI-24. Check raw file.');
        return
    end
end



end