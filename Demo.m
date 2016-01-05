% This is a demo for the elastic PC algorithm
% Please download and unzip the test data from http://www.fmrib.ox.ac.uk/analysis/netsim/sims.tar.gz
results = cell(28,50,5);
for i= 1:28
    disp(i);
    load(['sim',num2str(i),'.mat']);
    for j = 1:Nsubjects
        ts1=ts(((j-1)*Ntimepoints+1):(j*Ntimepoints),:);
        ts1=ts1-repmat(mean(ts1),Ntimepoints,1);
        results{i,j,1} = corr(ts1);  % full correlation
        results{i,j,2} = partialcorr(ts1); % partial correlation

% 		Elastic PC starts from 0.05 to 0.25
        [results{i,j,3},mpc_cube] = ElasticPC(ts1,0.05);
		[results{i,j,4},mpc_cube] = ElasticPC(ts1,0.15,0.05,results{i,j,3},mpc_cube);
		[results{i,j,5},mpc_cube] = ElasticPC(ts1,0.25,0.15,results{i,j,4},mpc_cube);
		
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the average c-sensitivity for each method
n_methods = 5;
cs = zeros(28,50,n_methods);
for i=1:28
    load(['sim',num2str(i),'.mat']);
    for j = 1:Nsubjects
        gtnet=squeeze(net(j,:,:));
        gtstructure = tril(abs(gtnet)+abs(gtnet)',-1);
        gtstructure(gtstructure~=0)=1;

        for k = 1:n_methods
            fc = abs(results{i,j,k});
            fc = fc.*(1-eye(size(fc))); % remove diagonal
            
%           Calculate the c-sensitivity
            TP = fc(gtstructure==1);
            FP = fc(tril((1-gtstructure),-1)==1);
            threshold = prctile(abs(FP),95);
            cs(i,j,k) = length(find(abs(TP)>threshold))/sum(gtstructure(:));    
        end
    end
end

% Average
ave_cs = zeros(28,n_methods);
for i=1:28
    for k = 1:n_methods
        ave_cs(i,k) = mean(cs(i,:,k));
    end
end
    




