function [ topGenes ] = top_genes_diff_simu(w,featuresnames)

%------- TopGenes :
normGene = zeros(size(w,1),1);
weightDiffexpr = zeros(size(w,1),1);
indgeneNDE = [];
for i = 1:size(w,1)
    normGene(i,1) = norm(w(i,:));
    weightDiffexpr(i,1) = max(w(i,:))-min(w(i,:));
    if contains(featuresnames(i),'NDE')
        i2remove = indgeneNDE;
        indgeneNDE = [i2remove, i];
    end
end
weightDiffexpr(indgeneNDE) = [];
featuresnames(indgeneNDE) = [];
[diffSort,ind] = sort(weightDiffexpr,'descend');
if size(featuresnames,1) < size(featuresnames,2)
    featuresnames = featuresnames';
end
topGenes = table(featuresnames(ind), ind, diffSort);

end

