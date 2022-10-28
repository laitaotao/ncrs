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


function [ par res inx tim K] = mAGS(lim,data,M,model_type,score)
% input:
% lim (1x1) = Maximum CPU seconds allowed.
% data (dxn) = Input data of dimensionality d.
% M (1x1) = Maximum number of hypotheses to be generated.
% blksize (1x1) = Block size of Multi-GS.
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
degenmax = 100;  % Max number of inner loop to avoid degeneracy.

% if (mod(M,blksiz)~=0)
%     %error('Bad block size!');
% end

%-----------------
% Prepare storage.
%-----------------
n = size(data,2);
par = zeros(numpar,M);
res = zeros(n,M);
K = zeros(1,M);
resinx2 = zeros(n,M);
res2 = zeros(n,M);
inx = zeros(psize,M);
tim = zeros(1,M);

%-----------------
% Guided sampling.
%-----------------
%fprintf('Multi-GS sampling for %.2f seconds...',lim);

c = 1;
d_scale = 0.1;
blksiz = 100;
degenCount = 0;
time_for_sort = 0;
%tic;
lambda1   = 0.8;
numNeigh1 = 6;
pts = data';

runNum = 1000;

tic;
for m=1:M
    
    degencnt = 0;
    isdegen = 1;
    ff = 0;
    while (isdegen==1)&&(degencnt<=degenmax)
        % Increment degeneracy count.
        degencnt = degencnt + 1;
        
        if m<= runNum
            % Uniform sampling for the first block.
            [ pinx ] = randsample(n,psize);
            kk = 0;
        else
            % Weighted sampling
            
            itorCount = 0;
            seedinx = 0;
            if(rand < 0.4)  %0.5 2.64  0.6 2.68  0.4 2.60  0.45 2.69 0.3 2.77 (itorCount< 10)
                            %0.5 2.62   0.4 2.78   (seedinx == 0)
                while(itorCount< 20)%(itorCount< 10) %1 2  (seedinx == 0) //0.5 2.60  //1000 0.5 2.50 // 0.4 2.51  //0.6 2.63   //0.45 2.68
                    % //0.4 400 2.61  0.4 800  2.485
                    choice = randsample(n, 1);%randsample(n,1,true,qua_score);%
                    distance1 = (pts(:,1) - pts(choice,1)).^2 + (pts(:,2) - pts(choice,2)).^2;
                    [~, index_distance1] = sort(distance1);
                    distance2 = (pts(:,4) - pts(choice,4)).^2 + (pts(:,5) - pts(choice,5)).^2;
                    [~, index_distance2] = sort(distance2);
                    [p2, C] = point_cosF(index_distance1(1:numNeigh1+3), index_distance2(1:numNeigh1+3), lambda1, numNeigh1);
                    if(p2 == 1)
                        seedinx = choice;
                        break;
                    end
                    itorCount = itorCount + 1;
                end
            end
            if(seedinx == 0)
                seedinx = randsample(n, 1);
            end
            
            
            [ pinx sig_inx ff kk] = weightedSampling(n,psize,resinx2(:,1:win),m_now, d_scale,score, seedinx);
        end
        
        if(length(pinx)<psize)
            continue;
        end
        
        psub = data(:,pinx);
        
        if isempty(strfind(model_type,'fundamental'))
            % Fit the model, and check for degeneracy.
            if(length(pinx) == psize)
                isdegen = feval(degenfn,psub);
            end
        else
            % Check for degeneracy.
            [ isdegen F ] = feval(degenfn,psub);
        end
        % need delete
        if degencnt>=2
            degenCount = degenCount + 1;
        end
    end
    
    if isempty(strfind(model_type,'fundamental'))
        % Fit the model on the p-subset
        st = feval(fitfn,psub);
        % Compute residuals.
        ds = feval(resfn,st,data);
    else
        % Compute residuals.
        [ ds st ] = feval(resfn,F,data);
    end
    
    % Store.
    par(:,m) = st;
    res(:,m) = ds;
    K(m) = kk;
    inx(:,m) = pinx;
    m_vector = repmat(m,n,1);
    res2(:,c) = ds;
    resinx2(:,c) = m_vector;
    
    
    tim(1,m) = toc;
    tic;
    if(m~=1)
        tim(1,m) = tim(1,m) + tim(1,m-1);
    end
    
    if tim(1,m)>=lim
        par = par(:,1:m);
        res = res(:,1:m);
        K = K(1:m);
        inx = inx(:,1:m);
        tim = tim(:,1:m);
        break;
    end
    
    % Update indices of sorted hypotheses.
    if(m>=runNum)&&(mod(m,blksiz)==0)%&&(m>=1000)
        
        % Intersection width.
        win = round(d_scale * m);
        m_now = m;
        [ foo resinx] = sort(res2(:,1:c),2);
        if(win + blksiz > m)
            c = m;
        else
            c = win + blksiz;
        end
        
        res2 = foo(:,1:c);
        for p=1:n
            resinx2(p,1:c) = resinx2(p,resinx(p,1:c));
        end
    end
    c = c+1;
end
end


%-------------------
% Weighted sampling.
%-------------------
function [pinx sig_particles_index2 ff kk] = weightedSampling(n, psize, resinx, m_now, d_scale, score, seedinx)
% n (1x1)      = Size of data.
% psize (1x1)  = Size of p-subset.
% resinx (nxm) = Indices of sorted hypotheses for each datum.
% win (1x1)    = Intersection width.

% Storage.
pinx = zeros(psize,1);

% First index is found by uniform sampling.
%seedinx = randsample(n,1);

% seedLabel = label(seedinx);
% if(seedLabel~=0)
%     debug = 1;
% end
pinx(1,:) = seedinx;
ff=0;
now_inx = round(m_now*d_scale);
new_w_temp = computeIntersectionForIT(resinx(seedinx,1:now_inx)',resinx(:,1:now_inx)',[m_now, d_scale]);%computeweightModify(resinx(seedinx,:)',resinx',m_now, d_scale);%
[II2, EE2]=Entropy_Thresholding(new_w_temp, 2);
sig_particles_index2 = find(II2>EE2);
if isempty(sig_particles_index2)
    new_w = 0;
    kk = 0;
else
    ww1 = (new_w_temp(sig_particles_index2)/sum(new_w_temp(sig_particles_index2)));
    ww2 = score(sig_particles_index2);
    ww3 = ww2/sum(ww2);
    new_w = ww1.*ww3;
    kk = length(sig_particles_index2);
    %     new_w = ww1;
end
w_sum = sum(new_w);
if (w_sum > 0)&&(length(sig_particles_index2)> psize)
    %ppppp = randsample(length(new_w),psize,true,new_w);
    zz = cumsum(new_w/sum(new_w));
    [ppppp]= guidedrandsample2(zz,psize,new_w);
else
    ppppp = randsample(n,psize);
    zz = [];
    sig_particles_index2 = [];
end

for i=1:psize
    if (w_sum > 0)&&(length(sig_particles_index2)> psize)
        pinx(i,:) = sig_particles_index2(ppppp(i));
    else
        pinx(i,:) =ppppp(i);
    end
end
end