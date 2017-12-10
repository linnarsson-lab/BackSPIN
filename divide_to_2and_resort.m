function [sorted_data_resort1,genes_resort1,cells_resort1,gr1,gr2,genesgr1,genesgr2,sc1,sc2] = divide_to_2and_resort(sorted_data,flag_lastlev)

% calc correlation matrix for cells
Rcells = corr_mat(sorted_data);
Rgenes = corr_mat(sorted_data');
% look for the optimal breaking point
[N,~] = size(Rcells);
sc = zeros(N,1);
for i=1:N-1
    i
    %     tmp1 = Rcells(1:i,1:i);
    %     tmp2 = Rcells(i+1:N,i+1:N);
    %     sc(i) = (sum(tmp1(:))+sum(tmp2(:)))/ (i^2+(N-i)^2) ;
    if i==1
        tmp1 = Rcells(1:i,1:i);
        tmp1 = sum(tmp1(:));
        tmp2 = Rcells(i+1:N,i+1:N);
        tmp2 = sum(tmp2(:));
        sc(i) = (tmp1+tmp2)/ (i^2+(N-i)^2) ;
    else
        tmp1 = tmp1 + sum(Rcells(i,1:i)) + sum(Rcells(1:i,i));
        tmp2 = tmp2 - sum(Rcells(i,i:end)) - sum(Rcells(i:end,i));
        sc(i) = (tmp1+tmp2)/ (i^2+(N-i)^2) ;
    end
end


[~,bp1] = max(sc);
sc1 = Rcells(1:bp1,1:bp1);
sc1 = triu(sc1);
sc1 = mean(sc1(sc1(:)~=0));
sc2 = Rcells(bp1+1:N,bp1+1:N);
sc2 = triu(sc2);
sc2 = mean(sc2(sc2(:)~=0));
aveall = triu(Rcells);
aveall = mean(aveall(aveall(:)~=0));

max([sc1,sc2])/aveall
if max([sc1,sc2])/aveall>1.0;%;1.15
    
    
    % divide to two groups
    gr1 = 1:bp1;
    gr2 = bp1+1:N;
    % and assign the genes into the two groups
    mean_gr1 = mean(sorted_data(:,gr1),2);
    mean_gr2 = mean(sorted_data(:,gr2),2);
    [center_gr1,flip_flag1] = min([calc_loccenter(sorted_data(:,gr1),2),calc_loccenter(fliplr(sorted_data(:,gr1)),2)],[],2);
    [center_gr2,flip_flag2] = max([calc_loccenter(sorted_data(:,gr2),2),calc_loccenter(fliplr(sorted_data(:,gr2)),2)],[],2);
    sorted_data_tmp = sorted_data;
    sorted_data_tmp(flip_flag1==2,gr1) = fliplr(sorted_data(flip_flag1==2,gr1));
    sorted_data_tmp(flip_flag2==2,gr2) = fliplr(sorted_data(flip_flag2==2,gr2));
    loc_center = calc_loccenter(sorted_data_tmp,2);
%     loc_center = (center_gr1.*sum(sorted_data_to_min(:,gr1),2) + (center_gr2+bp1).*sum(sorted_data_to_min(:,gr2),2))./sum(sorted_data_to_min,2);
    
%     loc_center = calc_loccenter(sorted_data);
%     imax = zeros(size(loc_center));
    imax = ones(size(loc_center));%changed mar27 2017
    imax(loc_center<=bp1) = 1;
    imax(loc_center>bp1) = 2;
    % % % % % % % % % % % % % % % % % % % % % % % % % % % %     modify
%     if length(gr1)>1 | length(gr2)>1
%         loc_var_gr1 = calc_locvar(sorted_data(:,gr1));
%         loc_var_gr2 = calc_locvar(sorted_data(:,gr2));
%         [~,imax] = min([loc_var_gr2./mean_gr1,loc_var_gr1./mean_gr1],[],2);
%     else
%         [~,imax] = max([mean_gr1,mean_gr2],[],2);
%     end
    % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    
%         d = abs(mean_gr1-mean_gr2);
%         if length(d)>20
%             for i=1:length(d)
%                 if d(i)<0.2
%                     in = Rgenes(i,:)>prctile(Rgenes(i,:),100-100*20/length(d));
%                     mean_gr1(i) = mean(sum(sorted_data(in,gr1)),2);
%                     mean_gr2(i) = mean(sum(sorted_data(in,gr2)),2);
%                 end
%             end
%         end
%         [~,imax] = max([mean_gr1,mean_gr2],[],2);
    
    
    genesgr1 = find(imax==1);
    genesgr2 = find(imax==2);
    if isempty(genesgr1)
        [~,in] = max(mean_gr1);
        genesgr1 = in;
        genesgr2 = setdiff(genesgr2,in);
    elseif isempty(genesgr2)
        [~,in] = max(mean_gr2);
        genesgr2 = in;
        genesgr1 = setdiff(genesgr1,in);
    end
    
    % resort each group
    datagr1 = sorted_data(genesgr1,gr1);
    datagr1 = datagr1 - repmat(mean(datagr1,2),1,length(gr1));
    if min(size(datagr1))>1
        if flag_lastlev~=1
            [~,genesorder1,cellorder1] = sort_and_noplot(datagr1,1,3,2,0);% sort_and_noplot(datagr1,1,1,8,0)
        else
            [~,genesorder1,cellorder1] = sort_and_noplot(datagr1,1,1,2,0);% sort_and_noplot(datagr1,1,1,8,0)
        end
    elseif length(genesgr1)==1
        genesorder1 = 1;
        [~,cellorder1] = sort(datagr1(1,:));
    elseif length(gr1)==1
        cellorder1 = 1;
        [~,genesorder1] = sort(datagr1(:,1));
    end
    
    datagr2 = sorted_data(genesgr2,gr2);
    datagr2 = datagr2 - repmat(mean(datagr2,2),1,length(gr2));
    if min(size(datagr2))>1
        if flag_lastlev~=1
            [~,genesorder2,cellorder2] = sort_and_noplot(datagr2,1,3,2,0);% sort_and_noplot(datagr2,1,1,8,0);
        else
            [~,genesorder2,cellorder2] = sort_and_noplot(datagr2,1,1,2,0);% sort_and_noplot(datagr2,1,1,8,0);
        end
    elseif length(genesgr2)==1
        genesorder2 = 1;
        [~,cellorder2] = sort(datagr2(1,:));
    elseif length(gr2)==1
        cellorder2 = 1;
        [~,genesorder2] = sort(datagr2(:,1));
    end
    
    
    % contcatenate cells and genes order
    genes_resort1 = [genesgr1(genesorder1); genesgr2(genesorder2) ];
    cells_resort1 = [gr1(cellorder1)';gr2(cellorder2)'];
    sorted_data_resort1 = sorted_data(genes_resort1,cells_resort1);
    
    genesgr1 = 1:length(genesgr1);
    genesgr2 = length(genesgr1)+1:length(sorted_data(:,1));
    
else
    
    sorted_data_resort1 = sorted_data;
    genes_resort1 = 1:length(sorted_data(:,1));
    cells_resort1 = 1:length(sorted_data(1,:));
    gr1 = 1:length(sorted_data(1,:));
    gr2 = [];
    genesgr1 = 1:length(sorted_data(:,1));
    genesgr2 = [];
    sc1 = Rcells;
    sc1 = triu(sc1);
    sc1 = mean(sc1(sc(:)~=0));
    sc2 = 0;
end





