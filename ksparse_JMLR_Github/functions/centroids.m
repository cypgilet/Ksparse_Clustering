function [mu] = centroids(XW,Y,k)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

mu = zeros(k,size(XW,2));
for i = 1:k
    for j = 1:size(XW,2)
        Ci = XW(Y==i,:);
        mu(i,j) = mean(Ci(:,j));
    end
end

end

