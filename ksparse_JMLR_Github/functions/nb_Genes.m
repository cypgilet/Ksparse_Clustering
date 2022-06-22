function [nbG, indGene_w] = nb_Genes(w)
% Return the number of selected genes from the matrix w

d = size(w,1);
ind_genes = zeros(d,1);

for i = 1:d
    if norm(w(i,:))>0
        ind_genes(i,1) = 1;
    end
end

indGene_w = find(ind_genes==1);
nbG = sum(ind_genes);

end

