%% script to check motion segmentation performance

clear; close all;

max_NumHypoPerFrame = 100;
FrameGap = 1;
model_type = lower('homography');

AlphaRange = 5:20;

%% Load Seq Information
temp = load('../../Data/SeqList.mat');
SeqList = temp.SeqList;

seq_range = 1:length(SeqList);

valid_alpha = [];
avg_error = [];
med_error = [];

for Alpha = AlphaRange
    error = [];
    for iter = 1:100
        %% Load Segmentation Result
        result_path = fullfile('../../Results/MoSeg/',model_type, '/', num2str(iter));
        
        result_filepath = fullfile(result_path,sprintf('Error_RandSamp_nhpf-%d_alpha-%g.mat',...
            max_NumHypoPerFrame,Alpha));
        
        if ~exist(result_filepath,'file')
            continue;
        else
            temp = load(result_filepath);
        end
        
        error = [error ; temp.error];      
    end
    if(~isempty(error))
        avg_error = [avg_error mean(mean(error,1))];
        med_error = [med_error median(median(error,1))];
        valid_alpha = [valid_alpha Alpha];
    end
end

% avg_error = mean(error,2);
% med_error = median(error,2);

figure;
colors = colormap(hsv(3));
plot(valid_alpha,100*avg_error,'color',colors(1,:),'marker','x'); hold on;
plot(valid_alpha,100*med_error,'color',colors(2,:),'marker','x'); hold on;

xlabel('alpha');
ylabel('error (%)');

grid on;

legend({'Mean Error','Median Error'});

title('Fundamental Matrix Motion Segmentation Performance');
