% create Z a tree matrix in matlab form. D should be a distance matrix. if
% D is a correlation matrix send 1-D so that diagonal=0. Z is M-1x3 matrix
% where M is the size of D. D(:,1:2) contain the two clusters merged in
% each step and the Z(:,3) contain the distance.
function Z = tree_for_spin_clust(D)

Z = zeros(length(D)-1,3);
diag1 = diag(D,1);
clus_in = 1:length(D);
D1 = D;
for i=1:length(D)-1
    i
    [val,im] = min(diag1);
    Z(i,:) = [clus_in(im),clus_in(im+1),val];
%     if im~=1
%         tmp1 = mean(D1(im:im+1,:));
%         tmp2 = mean(D1(im-1:im,:));
%         tmp1 = max(D1(im:im+1,:));
%         tmp2 = max(D1(im-1:im,:));
%         tmp1 = geomean(D1(im:im+1,:));
%         tmp2 = geomean(D1(im-1:im,:));
%     else
%         tmp1 =1;
%         tmp2 = mean(D1(im:im+1,:));
%         tmp2 = max(D1(im:im+1,:));
%         tmp2 = geomean(D1(im:im+1,:));
%     end
%     D1(im,:) = [tmp1(1:im-1),1,tmp2(im+1:end)];
%     D1(:,im) = D1(im,:);    
%     D1(im,:) = mean(D1(im:im+1,:));
%     D1(:,im) = D1(im,:)';
    D1(:,im+1) = [];
    D1(im+1,:) = [];
    clus_in(im) = i + length(D);
    clus_in(im+1) =[];
    diag1 = diag(D1,1);
end

    
    









