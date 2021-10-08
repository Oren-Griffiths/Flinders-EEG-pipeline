
% a simple function to plot the output of AR calculations in EEGlab. 
% % indicative calls ot the function are set out below.
% eeglab;
% load('testEEG.mat');
% data = testEEG.data;
% chanlocs = testEEG.chanlocs;
% SVT = 100;
% Subject_Path = [fileparts(pwd) '\102\'];

% write this line in the command line. 
% DrawARfigs(data,SVT,chanlocs, Subject_Path);

function DrawARfigs(data,SVT,chanlocs, Subject_Path, imageType)
% imageType determines whether to draw png or pdf (pdf larger and will fail
% in parallel mode as doesn't have enough memory per core). 
samples = [1:size(data,2)];
% samples = ((samples-1)/256) - 0.2; % this is how you adjust for baseline.
threshold_pos = SVT*ones(1,size(data,2));
second_threshold_pos = 2*SVT*ones(1,size(data,2));
threshold_neg = -1*SVT*ones(1,size(data,2));
second_threshold_neg = -2*SVT*ones(1,size(data,2));
% parameters needed for looping
if size(data,1) < 60
    NoOfChans = 32;
else
    NoOfChans = 64;
end
NoOfEpochs = size(data,3);

% start drawing the figure.
figure;

for ThisChan = 1:NoOfChans    
    % nominate somewhere to draw.
    subplot(round(sqrt(NoOfChans))+1, round(sqrt(NoOfChans))+1,ThisChan);
    
    % collect data per channel. 
    one_channel = abs(squeeze(data(ThisChan,:,:)))'; % gives epochs x samples.
    % and only positive values (so can compare to one threshold).
    
    maxvals = max(one_channel,[],2); 
    % find max val per epoch, in (n x 1) list of  epochs
    [~, idx] = sort(maxvals); % sort list of maxvals.
    
    % reorder one_channel by sorted size
    for ThisEpoch = 1:NoOfEpochs
    one_channel_sorted(ThisEpoch,:) = one_channel(idx(ThisEpoch),:);
    end
    
    % find nth epoch
    boundaries = [1 round(0.1*NoOfEpochs) , round(0.2*NoOfEpochs), ...
        round(0.3*NoOfEpochs), round(0.4*NoOfEpochs), ...
        round(0.5*NoOfEpochs),  round(0.6*NoOfEpochs), ...
         round(0.7*NoOfEpochs),  round(0.8*NoOfEpochs), ...
          round(0.9*NoOfEpochs) NoOfEpochs];
    
    labels = {'1%', '10%', '20%', '30%', '40%', '50%', '60%', '70%', ...
        '80%', '90%', '100%'};
    
    % in case any elements are zero because <10 epochs, then limit to
    % non-zero integers
    boundaries = boundaries(boundaries > 0);
    labels = labels(boundaries > 0);
    
    if ~isempty(boundaries) % skip this facet if no relevant epochs.
        % list of epochs that correspond to each decile in terms of max value.
        deciles = one_channel_sorted(boundaries,:);
        
        hold on
        for ThisDecile = size(deciles,1):-1:1 % start of line drawing loop
            if max(deciles(ThisDecile,:) > SVT)
                % red for those that break threshold.
                line(samples, deciles(ThisDecile,:)-100*ThisDecile, 'Color',[1 0 0]);
                text(samples(end),0-100*ThisDecile,...
                    labels{ThisDecile}, 'FontSize',8, 'Color', 'r');
            else
                % black to blue for those that keep under threshold.
                line(samples, deciles(ThisDecile,:)-100*ThisDecile, 'Color',[0 0 ThisDecile/size(deciles,1)]);
            end
        end % of line drawing loop.
        
        hold off
    end
    % and add some labels.
    title(chanlocs(ThisChan).labels);
    % ylim([-3*SVT 3*SVT]);
    xlim([1 size(data,2)]);
    set(gca,'FontSize',8);
    disp(num2str(ThisChan));
end

if strcmp(imageType, 'pdf')
    % if explicitly requested a high quality vector (pdf) image then do
    % that.
    out_filename = [Subject_Path 'Figures\X6_',num2str(SVT), '_ARbyDecile.pdf'];
else % just give a shitty low-res png
    out_filename = [Subject_Path 'Figures\X6_',num2str(SVT), '_ARbyDecile.png'];
end
% make an enormous figure 10 inches by 10 inches.
% ideally should be 25 inches square, but often can be too slow. 
if strcmp(imageType, 'pdf')
    set(gcf,'PaperPositionMode','manual','PaperUnits','Inches','PaperSize',[25 25], 'PaperPosition', [0 0 25 25] );
    tic;
    saveas(gcf,out_filename,'pdf');
    close(gcf);
    display(['Time taken to draw highres PDF is ' num2str(toc)]);
else % just save a quick png
    set(gcf,'PaperPositionMode','manual','PaperUnits','Inches','PaperSize',[25 25], 'PaperPosition', [0 0 25 25] );
    tic;
    saveas(gcf,out_filename,'png');
    close(gcf);
    display(['Time taken to draw lowres PNG is ' num2str(toc)]);
end


% calculate proportion of epochs that break threshold (two loops)
for ThisChan = 1:NoOfChans
    tempdata = squeeze(abs(data(ThisChan,:,:)));
    maxes = max(tempdata, [], 1);
    idx = (maxes>SVT);
    ARfails(ThisChan) = mean(idx);
end

% draw a figure with proportion of epochs failed
figure;
plot(ARfails);

% put the correct labels on x-axis
ChanLabels = extractfield(chanlocs(1:NoOfChans),'labels');

for ThisChan = 1:NoOfChans
text(ThisChan,ARfails(ThisChan),...
 ChanLabels{ThisChan}, 'FontSize',18, 'Color', 'k');
end

xlim([1 NoOfChans]);
xticks([1:NoOfChans]);
xticklabels(ChanLabels);
ylabel('PropnEpochsRejected');
set(gca,'FontSize',18);
% save that figure
set(gcf,'PaperPositionMode','manual','PaperUnits','Inches','PaperSize',[25 25], 'PaperPosition', [0 0 25 25] );
out_filename = [Subject_Path 'Figures\X6_',num2str(SVT) , '_ARperChannel.pdf'];
saveas(gcf,out_filename,'pdf');
close(gcf);

end