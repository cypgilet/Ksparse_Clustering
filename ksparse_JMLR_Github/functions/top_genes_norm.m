function [ topGenes ] = top_genes_norm(w,featuresnames)

%------- TopGenes :
normGene = zeros(size(w,1),1);
for i = 1:size(w,1)
    normGene(i,1) = norm(w(i,:));
end
[normSort,ind] = sort(normGene,'descend');
ind = ind(1:sum(normGene~=0));
normSort = normSort(1:sum(normGene~=0));
if size(featuresnames,1) < size(featuresnames,2)
    featuresnames = featuresnames';
end
topGenes = table(featuresnames(ind), ind, normSort);

end

