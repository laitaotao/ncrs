%% script to check motion segmentation performance

clear; close all;

max_NumHypoPerFrame = 100;
FrameGap = 1;
model_type = lower('fundamental');

AlphaRange = 11;%5:15;

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
        avg_error = mean(error,1)*100;
        med_error = median(error,1)*100;
        valid_alpha = [valid_alpha Alpha];
    end
end

save('F_results.mat','avg_error','med_error');
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
