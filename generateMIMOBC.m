function generateMIMOBC(M,N,K,num_H) % K- users, M - transmit antenna, N- receive antenna
tic
P=10;% set the power
rng('default');
X=zeros(M*N*K,num_H);% creat initial X and Y
Y=zeros(M*N*K,num_H);
for loop = 1:num_H % training samples
    capacity=0;
    Covar=zeros(M,N,K);% initialization
    Downlink=zeros(M,N,K);%initialization
    % generate channel data
    H = 1/sqrt(2) * (randn(M,N,K) + i * randn(M,N,K));
    % shape(H): (2,2,3)
    temp_H = reshape(H,M*N*K,1);% reshape the matrix to array (12,1)
    % use IWF method to get covariance matrix
    [capacity, Covar] = iterative_waterfill(H,P,50);
    % use transfer covariance matrix to downlink matrix
    Downlink = MAC_to_BC(H, Covar);
    %  reshape the matrix to array (12,1)
    temp_Downlink = reshape(Downlink,M*N*K,1);
    % store one sample, data X, Y
    X(:,loop)=temp_H;
    Y(:,loop)=temp_Downlink;
    if mod(loop,5000)==0
        fprintf('.');
    end
end
% save the dataset
save(sprintf('MIMOBC%d_%d.mat',num_H, K),'X','Y');
fprintf('Generate Done! \n');
toc
end