function [sorted_data,genesorder,cellsorder,Z_cells,Z_genes] = sort_and_noplot(data,corr_dist_flag,sorttype_flag,sort_cyc,flag_tree)

if (nargin <3)
    error('At least two input arguments required. USAGE: em(x,k,min_alpha)');
elseif (nargin == 3)
    sort_cyc = 8;
elseif (nargin == 4)
    flag_tree = 0;
elseif (nargin >5)
    error('Too many input arguments. USAGE: em(x,k,min_alpha)');
end


%  sorttype_flag=1 sort by both genes and samples, =2 only genes, =3 only
%  samples

[N,M] = size(data);
if corr_dist_flag==1 | isempty(corr_dist_flag)
    if sorttype_flag==1 | sorttype_flag==3
        corrmat_s = corr_mat(data);%NxN
        genesorder = 1:N;
    end
    if sorttype_flag==1 | sorttype_flag==2
        corrmat_g = corr_mat(data');%MxM
    end
elseif corr_dist_flag==2
    if sorttype_flag==1 | sorttype_flag==3
        genesorder = 1:N;
        corrmat_s = 1-distance_mat(data);%NxN
    end
    if sorttype_flag==1 | sorttype_flag==2
        corrmat_g = 1-distance_mat(data');%MxM
    end
end

% sort_score_i1 = inf;
if sorttype_flag==1 | sorttype_flag==3
    wid_vec = [round(0.4*M):-round(0.01*M):1,1];
%     XI = randperm(M);
%     [~,mx] = max(corrmat_s);
    [~,XI] = sort(sum(corrmat_s));
    for i=1:length(wid_vec)
        i
        if i<=length(wid_vec)
            tmpmat = 1-corrmat_s(XI,:);
            tmpmat = tmpmat(:,XI);
            [sorted_ind,sort_score] = sort_mat_by_neighborhood(tmpmat, wid_vec(i), sort_cyc,2);
            XI=XI(sorted_ind);
        else
            tmpmat = 1-corrmat_s(XI,:);
            tmpmat = tmpmat(:,XI);
            [sorted_ind,sort_score] = sort_mat_by_neighborhood(tmpmat, wid_vec(i), sort_cyc,2);
            XI=XI(sorted_ind);
        end
        %         if sum(diff(sort_score)<0)<sort_score_i1
        %             XI=XI(sorted_ind);
        %             sort_score_i1 = sum(diff(sort_score)<0);
        %         end
    end
    cellsorder = XI;
else
    cellsorder = [1:M];
end

% sort_score_i1 = inf;
if sorttype_flag==1 | sorttype_flag==2
    wid_vec = [round(0.4*N):-round(0.01*N):1,1];
%     XI = randperm(N);
%     [~,mx] = max(corrmat_g);
    [~,XI] = sort(sum(corrmat_g));
    for i=1:length(wid_vec)
        i
        if i<=length(wid_vec)
            tmpmat = 1-corrmat_g(XI,:);
            tmpmat = tmpmat(:,XI);
            [sorted_ind] = sort_mat_by_neighborhood(tmpmat, wid_vec(i), sort_cyc,2);
            XI=XI(sorted_ind);
        else
            tmpmat = 1-corrmat_g(XI,:);
            tmpmat = tmpmat(:,XI);
            [sorted_ind] = sort_mat_by_neighborhood(tmpmat, wid_vec(i), sort_cyc,2);
            XI=XI(sorted_ind);
        end
        %         if sum(diff(sort_score)<0)<sort_score_i1
        %             XI=XI(sorted_ind);
        %             sort_score_i1 = sum(diff(sort_score)<0);
        %         end
    end
    genesorder = XI;
else
    genesorder = [1:N];
end

sorted_data = data(genesorder,cellsorder);

if flag_tree==1
    Z_cells = tree_for_spin_clust(1-corrmat_s(cellsorder,cellsorder));
    Z_genes = tree_for_spin_clust(1-corrmat_g(genesorder,genesorder));
end



