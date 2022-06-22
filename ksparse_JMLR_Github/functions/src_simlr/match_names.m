function [labels] = match_names(idx1,idx2,k)
% Match estimated class names to real class names
%----- INPUT
% idxR : real labels
% idx  : estimated labels
% k    : number of class
%----- OUTPUT
% label : estimated labels with matching class names


P = perms(1:k);
Tacc = zeros(size(P,1),2);


for i = 1:size(P,1)
    echange = P(i,:);
    idxT = zeros(size(idx2));
    for j=1:k
        idxT(find(idx2==j)) = echange(j);
    end
    [ACC_glob, tab_acc] = compute_accuracy(idx1,idxT,k);
    Tacc(i,1) = ACC_glob;
    Tacc(i,2) = mean(tab_acc);
end


ind = find(Tacc(:,2)==max(Tacc(:,2)));
if length(ind)>1
    T = Tacc(ind,:);
    ind = find(T(:,1)==max(T(:,1)));
end
echange = P(ind(1),:);


idxN = zeros(size(idx2));
for j=1:k
    idxN(find(idx2==j)) = echange(j);
end

labels = idxN;

end