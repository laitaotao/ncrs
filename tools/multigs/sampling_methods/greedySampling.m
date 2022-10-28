%-------------------------------------------------------------------------
% Accelerated Hypothesis Generation for Multi-Structure Robust Fitting
%-------------------------------------------------------------------------
% The demo code in this package implements the guided-sampling method for
% multi-structure robust fitting proposed in:
%
% T.-J. Chin, J. Yu and D. Suter
% Accelerated Hypothesis Generation for Multi-Structure Robust Fitting
% In Proc. European Conf. on Computer Vision, Crete, Greece, 2010.
%
% T.-J. Chin, J. Yu and D. Suter
% Accelerated Hypothesis Generation for Multi-Structure Data via Preference Analysis
% To appear in IEEE Trans. on Pattern Analysis and Machine Intelligence.
%
% Copyright (c) 2010 Tat-Jun Chin and Jin Yu
% School of Computer Science, The University of Adelaide, South Australia
% http://www.cs.adelaide.edu.au/~{tjchin,jinyu}
%
% The program is free for non-commercial academic use. Any commercial use
% is strictly prohibited without the authors' consent. Please acknowledge
% the authors by citing the above paper in any academic publications that
% have made use of this package or part of it.
%
% If you encounter any problems or questions please email to
% tjchin@cs.adelaide.edu.au.
%
% This program makes use of Peter Kovesi and Andrew Zisserman's MATLAB
% functions for multi-view geometry
% (http://www.csse.uwa.edu.au/~pk/Research/MatlabFns/
%  http://www.robots.ox.ac.uk/~vgg/hzbook/code/).


function [ par res inx tim ] = greedySampling(lim,data,numHypo,model_type)
% input:
% lim (1x1) = Maximum CPU seconds allowed.
% data (dxn) = Input data of dimensionality d.
% M (1x1) = Maximum number of hypotheses to be generated.
% model_type (string) = Type of model to be estimated.
%
% output:
% par (dxM) = Parameters of the putative models.
% res (nxM) = Residuals as measured to the putative models.
% inx (pxM) = Indices of p-subsets
% tim (1xM) = CPU time for generating each model.



%---------------------------
% Model specific parameters.
%---------------------------
[ fitfn resfn degenfn psize numpar ] = getModelParam(model_type);


%------------------------
% Check other parameters.
%------------------------
degenmax = 10;  % Max number of inner loop to avoid degeneracy.

%-----------------
% Prepare storage.
%-----------------
n = size(data,2);
par = zeros(numpar,numHypo);
res = zeros(n,numHypo);
inx = zeros(psize,numHypo);
tim = zeros(1,numHypo);

% if(n/100 < 300)
%     k = 10;
% else
%     k = 20;
% end
k = 10;

Threshold = 3.5;
eResiduals = 0.0005;

USE_GRAOUSE = 0;
%-----------------
% Random sampling.
%-----------------
fprintf('Greedy sampling for %.2f seconds...',lim);
t0 = cputime;

dataOrg = data;
y1 = dataOrg(1:3, :);
y2 = dataOrg(4:6, :);

dat_img_1 = normalise2dpts(y1);
dat_img_2 = normalise2dpts(y2);
data = [ dat_img_1; dat_img_2 ];

[LPMCorrectIndex] = LPM(dataOrg);%1:size(w_Ndata,2);%

tim = zeros(1,1);
hypo_count=1;
X = data;
%tic;
X_inliers = X(1:2, LPMCorrectIndex);
sigmaExp = 0.5;
[nearPtsTab] = calcNearPtsTab(X_inliers, 'exp', sigmaExp);
remainingPoints = 1:length(LPMCorrectIndex);
Converged = 0;

X_inliers = X(:, LPMCorrectIndex);
while(Converged==0)
    
    %     [Xs,Is] = datasample(X, floor(SampFrac*N), 2, 'Replace', false, 'Weights', W);
    %     Ws = W(Is);
    if (length(remainingPoints)<= k)
        remainingPoints = 1:length(LPMCorrectIndex);
    end
    Xout = X_inliers(:,remainingPoints);
    %run HMSS on selcted sample
    nearPtsT = nearPtsTab(remainingPoints,remainingPoints);
    choice = randperm(length(remainingPoints),1);
    nearPointsCdf = nearPtsT(choice,:);% + nearPtsTab2(ind,:))/2;
    nearPointsCdf(choice) = 0;
    nearPointsCdf = nearPointsCdf / sum(nearPointsCdf);
    nearPointsCdf = cumsum(nearPointsCdf);
    rndSub = rand(psize,1);
    [dum, pinxSub] = histc(rndSub,[0 nearPointsCdf]);
    
    [theta_f, sigma_f, ins, ~,  pinx] = FLKOSfitArbitraryModel(Xout, k, model_type, Threshold, pinxSub, eResiduals);
    identPoints = remainingPoints(ins);
    if( length(identPoints) >psize )
        remainingPoints = setdiff(remainingPoints ,identPoints);
    end
    %get the residual for all data points
    [ ht st ] =feval(resfn, theta_f, X);
    
    par(:,hypo_count) = st;
    res(:,hypo_count) = ht;
    inx(:,hypo_count) = pinx;
     
    if (USE_GRAOUSE ==0)
        if hypo_count==numHypo
            Converged = 1;
        end
    end
    hypo_count = hypo_count+ 1;
end
par = par(:,1:(hypo_count-1));
res = res(:,1:(hypo_count-1));
inx = inx(:,1:(hypo_count-1));
fprintf('done \n');
end
