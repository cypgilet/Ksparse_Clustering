function [labels,T,mu] = spcl_normalized_sym(X,k,sigma,nb_kms)
n = size(X,1);

W = squareform(pdist(X,'squaredeuclidean'));
W = exp(-((W)/(2*sigma^2)));
for i = 1:n
    W(i,i) = 0;
end

D = zeros(n,n);
for  i = 1:n
    D(i,i) = sum(W(i,:));
end
D2 = diag(1./(sqrt(diag(D))));
L = eye(n,n)-D2*W*D2;

%[vect_p,val_p] = eigs(L,k,'sm');
[vect_p,val_p]=eig(L);
[~,ind] = sort(diag(val_p));
vect_p = vect_p(:,ind);

U = vect_p(:,1:k);
T = zeros(n,k);
for i = 1:n
    T(i,:) = U(i,:)/norm(U(i,:),2); 
end

[Y,mu] = kmeans(T,k,'Replicates',nb_kms);

labels = Y;

end

