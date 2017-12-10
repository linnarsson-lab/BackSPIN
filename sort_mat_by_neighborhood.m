function [sorted_ind,sort_score] = sort_mat_by_neighborhood(dist_matrix, wid,N,w_pow)

M = length(dist_matrix(:,1));

XI=1:M;
XIbest = XI;
fitscor = sum(diag(dist_matrix(XI,XI),1));
for i=1:N
    tmpmat = dist_matrix(XI,:);
    tmpmat = tmpmat(:,XI);
    %     [sorted_ind] = side2side_sort(dist_matrix);
    [sorted_ind, sort_score, mismatch] = neighborhood_sort(tmpmat, wid,w_pow);
%     XI=XI(sorted_ind);
    %     if sum(mismatch(:))<sort_score_i1
    %         XI=XI(sorted_ind);
    %         sort_score_i1 = sum(diff(sort_score)<0);
    %     end
%     sum(diag(dist_matrix(XI,XI),1))
    
    s1 = sum(diag(dist_matrix(XI(sorted_ind),XI(sorted_ind)),1));
    XI=XI(sorted_ind);
%     if s1<fitscor
%         XIbest = XI;
%         fitscor =s1
%     end
end

sorted_ind = XI;%XIbest;

