% This code need Matlab R2018a
% load the DNN model which is trained by python.
modelfile = 'model5layers1r001.json'; % load DNN structure
weights = 'model5layers1r001.h5'; % load weights and bias in DNN model
net = importKerasNetwork(modelfile,'WeightFile',weights, ...
      'OutputLayerType','regression');
  
% layers = importKerasLayers(modelfile,'ImportWeights',true, ...
%       'WeightFile',weights,'OutputLayerType','regression')
sample=5000; % set test samples 
T1=0;
T2=0;
T3=0;
rate_nn_sum=[];
rate_iwf_sum=[];
K=3;
M=2;
N=2;
for loop=1:sample
    temp_H2d=zeros(24,1);
    % generate channel matrix H 
    H = 1/sqrt(2) * (randn(M,N,K) + 1i * randn(M,N,K));
    temp_H = reshape(H,M*N*K,1); % reshape H
    for l=1:M*N*K
        temp_H2d(l,:)=real(temp_H(l)); % transform complex number to real 
        temp_H2d(12+l,:)=imag(temp_H(l));% [R1 R2 R3 I1 I2 I3] format
    end
    % regenerate test dataset
    temp_H4d=temp_H2d';
    tic
    p_nn = predict(net,temp_H4d); % feed the data into DNN model
    T1=T1+toc; % caluculate DNN time
    p_nnM=zeros(M*N*K,1);
    for l=1:M*N*K     
        p_nnM(l)=p_nn(l)+1i*p_nn(12+l);
    end
    p_nnM=reshape(p_nnM,M,N,K);  % transform DNN vector output to covariance matrix 
    rate_nn=0;
    % begin to calculate DNN sum rate capacity
    for index=1:K
        if index==1
            rate_nn = rate_nn + real(log2(det(eye(M) +  H(:,:,index)' *p_nnM(:,:,index)* H(:,:,index))));
        else
        p_nnM(:,:,index)=p_nnM(:,:,index)+p_nnM(:,:,index-1);
        rate_nn = rate_nn + real(log2(det(eye(M) +  H(:,:,index)' *p_nnM(:,:,index)* H(:,:,index))/det(eye(M)+ H(:,:,index)' *p_nnM(:,:,index-1)* H(:,:,index))));
        end
    end
    rate_nn_sum=[rate_nn_sum rate_nn];% save sum rate value of each test sample 
    
    tic
    [capacity,Covar] = iterative_waterfill(H,10,50);
    Downlink = MAC_to_BC(H, Covar);% obtain covariance matrix using IWF algorithm
    T2=T2+toc; % calculate IWF time 
    rate_iwf=0;
    % calculate IWF sum rate capacity 
    for index=1:K
        if index==1
            rate_iwf = rate_iwf + real(log2(det(eye(M) +  H(:,:,index)' *Downlink(:,:,index)* H(:,:,index))));
        else
        Downlink(:,:,index)=Downlink(:,:,index)+Downlink(:,:,index-1);
        rate_iwf = rate_iwf + real(log2(det(eye(M) +  H(:,:,index)' *Downlink(:,:,index)* H(:,:,index))/det(eye(M)+ H(:,:,index)' *Downlink(:,:,index-1)* H(:,:,index))));
        end
    end
    rate_iwf_sum=[rate_iwf_sum rate_iwf];
%    % save IWF sum rate capacity of each test sample

end

fprintf('Testing performance: %.2f%% sum-rate in %.2f%% time\n',sum(rate_nn_sum./rate_iwf_sum)/sample*100, T1/T2*100);
fprintf('DNN sum rate: %.2f s\n  IWF sum rate: %.2f s\n ',sum(rate_nn_sum)/sample, sum(rate_iwf_sum)/sample);
fprintf('DNN time: %.2f s\n IWF time: %.2f s\n ', T1, T2);

% %%%%%% PART E: Plot Testing Performance %%%%%
disp('####### PART E: Plot Testing Performance #######');
figure(1)
cdfplot(rate_nn_sum)
hold on;
cdfplot(rate_iwf_sum)
hold on;
legend('DNN','IWF');
xlabel('rate');
ylabel('cumulative probability');
% savefig(sprintf('DNN_CDF2d_%d_%d_%d_%d',K,10000,ENCODE1,neurons));

figure(2)
[counts,centers]=hist(rate_nn_sum);
bar(centers,counts,'FaceColor',[0,1,1]);

hold on;
[counts2,centers2]=hist(rate_iwf_sum);
bar(centers2,counts2,'FaceColor',[0,1,0]);
hold on;
legend('DNN','IWF');
xlabel('rate');
ylabel('samples');