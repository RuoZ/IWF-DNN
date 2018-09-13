function changeComplex(X,Y,M,N,K,num_H) 
% use this function to transfer complex number to real part and imaginary
% part
% final formate: ABAB(real1,imag1,real2,imag2...) 
tic
rng('default');
XX=zeros(2*M*N*K,num_H);
YY=zeros(M*N*K,num_H);
for loop = 1:M*N*K % training samples
    XX(loop*2-1,:)=real(X(loop,:));
    XX(loop*2,:)=imag(X(loop,:));
    YY(loop*2-1,:)=real(Y(loop,:));
    YY(loop*2,:)=imag(Y(loop,:));
%     Y(:,loop)=temp_Downlink;
end
save(sprintf('MIMOComplex%d_%d.mat',num_H, K),'XX','YY');
fprintf('Generate Done! \n');
toc
end