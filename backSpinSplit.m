function [dataout_sorted,genes_order,cells_order,genes_gr_level,cells_gr_level,cells_gr_level_sc,genes_bor_level,cells_bor_level] = backSpinSplit(data,numLevels,issort)



genes_bor_level = cell(numLevels,1);
cells_bor_level = cell(numLevels,1);
[N,M] = size(data);
genes_order = [1:N]';
cells_order = [1:M]';
if ~isempty(issort)
    if issort==0
        [~,genesorder,cellsorder] = sort_and_noplot(data,1,1,20,0);%sorting by spin first
        genes_order = genes_order(genesorder);
        cells_order = cells_order(cellsorder);
    end
end

genes_gr_level = ones(N,numLevels+1);
cells_gr_level = ones(M,numLevels+1);
cells_gr_level_sc = zeros(M,numLevels+1);
for i=1:numLevels
    k=0;
    for j=1:length(unique(cells_gr_level(:,i)));%  2^(i-1)
        g_settmp = find(genes_gr_level(:,i)==j);
        c_settmp = find(cells_gr_level(:,i)==j);
        datatmp = data(genes_order(g_settmp),cells_order(c_settmp));
        if length(g_settmp)>1 & length(c_settmp)>1
            [sorted_data_resort1,genes_resort1,cells_resort1,gr1,gr2,genesgr1,genesgr2,sc1,sc2] = ...
                divide_to_2and_resort(datatmp,i==numLevels);
            if ~isempty(gr2)
                genes_order(g_settmp) = genes_order(g_settmp(genes_resort1));
                cells_order(c_settmp) = cells_order(c_settmp(cells_resort1));
                genes_gr_level(g_settmp(genesgr1),i+1) = k+1;
                genes_gr_level(g_settmp(genesgr2),i+1) = k+2;
                cells_gr_level(c_settmp(gr1),i+1) = k+1;
                cells_gr_level(c_settmp(gr2),i+1) = k+2;
                cells_gr_level_sc(c_settmp(gr1),i+1) = sc1;
                cells_gr_level_sc(c_settmp(gr2),i+1) = sc2;
                k = k+2;
            else
                genes_gr_level(g_settmp,i+1) = k+1;
                cells_gr_level(c_settmp,i+1) = k+1;
                cells_gr_level_sc(c_settmp,i+1) = cells_gr_level_sc(c_settmp,i);
                k = k+1;
            end
        else
            genes_gr_level(g_settmp,i+1) = k+1;
            cells_gr_level(c_settmp,i+1) = k+1;
            cells_gr_level_sc(c_settmp,i+1) = cells_gr_level_sc(c_settmp,i);
            k = k+1;
        end
    end
    genes_bor_level{i} = find(diff(genes_gr_level(:,i+1))>0)+1;
    cells_bor_level{i} = find(diff(cells_gr_level(:,i+1))>0)+1;
end


dataout_sorted = data(genes_order,cells_order);













