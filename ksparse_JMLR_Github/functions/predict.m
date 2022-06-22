function [Ytest] = Chambolle_Predict(Xtest,W,mu)
%PROXIMAL_DOUGLAS_RACHFORD_PREDICT Summary of this function goes here
%   Detailed explanation goes here
k = size(mu,1);
Ytest = zeros(size(Xtest,1),1);

for i = 1:size(Xtest,1)
    distmu = zeros(1,k);
    XWi = Xtest(i,:)*W;
    for j = 1:k
        distmu(j) = norm(XWi-mu(j,:),1);
    end
    [~,Ytest(i)] = min(distmu);
end
end