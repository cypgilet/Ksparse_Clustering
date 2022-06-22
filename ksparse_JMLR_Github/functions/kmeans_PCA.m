function [score,Ykms] = kmeans_PCA(X,k,nb_components,nb_kms)
[~,score] = pca(X,'NumComponents',nb_components);
Ykms = kmeans(score,k,'Replicates',nb_kms);
end

