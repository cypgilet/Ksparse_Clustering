function [ACC_glob, tab_acc] = compute_accuracy(idxR,idx,k)
%----- INPUT
% idxR : real labels
% idx  : estimated labels
% k    : number of class
%----- OUTPUT
% ACC_glob : global accuracy
% tab_acc  : accuracy per class


% Global accuracy :
y = sum(idxR==idx);
ACC_glob = y/length(idxR);

% Accuracy per class :
tab_acc = zeros(1,k);
for j =1:k
    ind=find(idxR==j);
    tab_acc(1,j) = sum(idx(ind)==j)/length(ind);
end

end