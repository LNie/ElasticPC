function [mpc,mpc_cube] = ElasticPC(fData,alpha,varargin)

% The elastic PC algorithm for calculating minimum partial correlations (MPC)
% Writen by Lei NIE (nieleimail@gmail.com)
% 7 Feb. 2015

% Input:
% fData is a matrix of fMRI data; each column is the time series of a ROI
% alpha is the given significant level
% varargin{1}: previous alpha
% varargin{2}: previous mpc cube
% varargin{3}: previous mpc

% Output:
% mpc is a matrix contains corrected z-score of minimum partial correlations
% mpc_cube is a 3-D matrix that contains mpc for each order


[nSamples,nNodes] = size(fData);
full_Covariance = cov(fData);
mpc_cube = Inf*ones(nNodes,nNodes,nNodes-1);

if nargin == 2
    pre_alpha = 0;
    full_Correlation = corr(fData);
    full_Correlation = 0.5.*log((1+full_Correlation)./(1-full_Correlation));
    full_Correlation = sqrt(nSamples - 3)*abs(full_Correlation);
    full_Correlation(eye(nNodes,nNodes)==1) = 0;
    pre_mpc_cube = Inf*ones(nNodes,nNodes,nNodes-1);
    pre_mpc_cube(:,:,1) = full_Correlation;
    mpc = Inf*ones(nNodes,nNodes);
else
    pre_alpha = varargin{1};
    mpc = varargin{2};
    pre_mpc_cube = varargin{3};
end

mpc_cube(:,:,1) = pre_mpc_cube(:,:,1);
pre_cutoff = norminv(1 - 0.5*pre_alpha);
cutoff = norminv(1 - 0.5*alpha);
mpc_square = pre_mpc_cube(:,:,1);
pre_mpc_square = pre_mpc_cube(:,:,1);

for ord = 1:(nNodes-2)
    pre_mpc_square = min(pre_mpc_square,pre_mpc_cube(:,:,ord));
    skeleton = ones(nNodes,nNodes);
    skeleton(mpc_square<=cutoff) = 0;
    pre_skeleton = ones(nNodes,nNodes);
    pre_skeleton(pre_mpc_square<=pre_cutoff) = 0;
    
    mpc_square = min(pre_mpc_square,mpc_square);
    [X,Y] = find(skeleton); 
    for i = 1:length(Y)
        x = X(i);
        y = Y(i);
        tmp_vec = skeleton(:,y);
        tmp_vec(x) = 0;
        if sum(tmp_vec) >= ord
            pre_tmp_vec = pre_skeleton(:,y);
            if pre_tmp_vec(x) == 1
                pre_tmp_vec(x) = 0;
                com_vec = tmp_vec.*pre_tmp_vec;
                new_vec = tmp_vec.*(1-com_vec);
                if ord ==1 
                    allSets = find(new_vec);
                else
                    com_nebs = find(com_vec)';
                    new_nebs = find(new_vec)';
                    allSets = multi_combs([new_nebs com_nebs],length(new_nebs),ord);
                end
            else
                if ord ==1
                    allSets = find(tmp_vec);
                else
                    nebs = find(tmp_vec)';
                    allSets = multi_combs(nebs,length(nebs),ord);
                end
            end
            nSets = size(allSets,1);
            tmp_b = zeros(ord+2,2);
            tmp_b(1,1) = 1;
            tmp_b(2,2) = 1;
            for j = 1:nSets
                tmp_Covariance = full_Covariance([x,y,allSets(j,:)],[x,y,allSets(j,:)]);
                tmp_val = tmp_Covariance\tmp_b;
                pcor = tmp_val(1,2)/sqrt(tmp_val(1,1)*tmp_val(2,2));
                pcor = 0.5*log((1+pcor)/(1-pcor));
                score = sqrt(nSamples - ord - 3)*abs(pcor);                     
                if score < mpc_square(x,y)
                    mpc_square(x,y) = score;
                    mpc_square(y,x) = score;
                end
            end
        end
    end
    mpc_cube(:,:,ord+1) = mpc_square;
end
mpc = min(min(mpc_cube,[],3),mpc);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = multi_combs(v,k,m)
% v is a row vector that presents the element set
% take m elements from v
% at least one element of the first k has to be chosen

n = length(v);
if k == 0
    P = [];
elseif n == m
    P = v;
elseif m == 1
    P = v(1:k).';
else
    P = [];
    if m < n && m > 1
        for i = 1:min(k,(n-m+1))
            Q = multi_combs(v(i+1:n),n-i,m-1);
            P = [P; [v(ones(size(Q,1),1),i) Q]]; %#ok<AGROW>
        end
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
