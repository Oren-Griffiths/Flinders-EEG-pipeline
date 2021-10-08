
% quick script to visualize which channels are problematic for SVT artefact
% rejection. 

% This has all been added to X6 script now.  

eeglab;

load('testEEG.mat');
data = testEEG.data;
chanlocs = testEEG.chanlocs;
SVT = 100;
Subject_Path = [fileparts(pwd) '\102\'];

DrawARfigs(data,SVT,chanlocs, Subject_Path);

function DrawARfigs(data,SVT,chanlocs, Subject_Path)

samples = [1:size(data,2)];
threshold_pos = SVT*ones(1,size(data,2));
second_threshold_pos = 2*SVT*ones(1,size(data,2));
threshold_neg = -1*SVT*ones(1,size(data,2));
second_threshold_neg = -2*SVT*ones(1,size(data,2));
%
if size(data,1) < 60
    NoOfChans = 32;
else
    NoOfChans = 64;
end

NoOfChans = 4;

figure;
for ThisChan = 1:NoOfChans
    ARfails(ThisChan) = 0; % initialize a counter per channel.
    for ThisEpoch = 1:size(data,3)
       subplot(round(sqrt(NoOfChans))+1, round(sqrt(NoOfChans))+1,ThisChan);
        hold on
        if max(abs(data(ThisChan,:,ThisEpoch))) > SVT
            tempfig = plot(samples, data(ThisChan,:,ThisEpoch), '-r'); % red
            ARfails(ThisChan) = ARfails(ThisChan) + 1;
        else
            plot(samples, data(ThisChan,:,ThisEpoch), '-b'); % blue
        end
        % add threshold indicators
        plot(samples, threshold_pos,'-k', samples, threshold_neg, '-k', ...
            samples, second_threshold_pos,'--k', samples, second_threshold_neg, '--k');
        hold off
        drawnow;
    end
    set(gca,'FontSize',8);

    % calculate proportion failed.
    title(chanlocs(ThisChan).labels);
    ylim([-3*SVT 3*SVT]);
    xlim([1 size(data,2)]);
    disp(num2str(ThisChan));
end

out_filename = [Subject_Path 'Figures\ARFigure_SVT',num2str(SVT) , '_practice.pdf'];
% make an enormous figure 25 inches by 25 inches.
set(gcf,'PaperPositionMode','manual','PaperUnits','Inches','PaperSize',[25 25], 'PaperPosition', [0 0 25 25] );
saveas(gcf,out_filename,'pdf');
close(gcf);

figure;
plot(ARfails/size(data,3));

% put the correct labels on x-axis
JustLabels = extractfield(chanlocs(1:NoOfChans),'labels');
xlim([1 NoOfChans]);
xticks([1:NoOfChans]);
xticklabels(JustLabels);
ylabel('PropnEpochsRejected');
% save that figure
out_filename = [Subject_Path 'Figures\ARFigure_Rejection',num2str(SVT) , '_practice.pdf'];
saveas(gcf,out_filename,'pdf');

end