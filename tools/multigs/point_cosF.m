function [p,C]=point_cosF(neighborX,neighborY,lambda,K)
%% *************************************
% Multi-scale LPM.
% 3 scale: K1= K-4£» K2=K-2£» K3= K.
% C=sum(ci/Ki).
% lambda is ranged from  0 to 1
%% **************************************
[~, L] = size(neighborX);
C = 0;
if K<4
    Error(' K must not less than 4£¬ so that the multi-scale not less than [2 4 6]')
end
Km = K+2 : -2 : K-2;%:K-4;
M  = length(Km);

for KK = Km
    neighborX = neighborX(2:KK+1, :);
    neighborY = neighborY(2:KK+1, :);
    neighborIndex = [neighborX; neighborY];
    index = sort(neighborIndex);
    temp1 = diff(index);
    temp2 = (temp1 == zeros(size(temp1, 1), size(temp1, 2)));
    ni = sum(temp2); %  the number of common elements in the K-NN
    %% cost calculation
    %****
    c1 = KK-ni;
    %*****
    C =  C+ (c1)/KK;%
end
p = (C./M) <= lambda.*ones(1,L);%C/M

