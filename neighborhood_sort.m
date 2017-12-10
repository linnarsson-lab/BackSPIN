%%%%%%%%%%%%%% Neighborhood sort %%%%%%%%%%%%%%%%%%%
function [sorted_ind, sort_score, mismatch] = neighborhood_sort(dist_matrix, wid, w_pow, weights_mat, flag)
% [sorted_ind, sort_score, mismatch] = neighborhood_sort(dist_matrix, wid, perc, flag)
% 05/09/02 - A function that accepts a distance matrix to be sorted
% Uses matrix form of weights: A normalized gaussian running along the main
% diagonal.
% dist_matrix = matrix to be sorted
% wid = the variance of the Gaussian that determines the neighborhood (default is 1)
% perc = percentile of values above which moves are allowed (default is 0)
% flag = 1 draws the weight matrix
% sorted_ind = the sorted indices
% sort_score = the scores
% mismatch = the smoothed matrix.


% Prepare the weights, must be same length as input matrix.
mat_size = length(dist_matrix);
if (nargin < 2)
    wid = 1;
else if (nargin < 4)
		[i,j] = meshgrid(1:mat_size,1:mat_size);
		weights_mat = abs(exp(-((i-j).^w_pow)/wid/mat_size)) - 1e-8;%elimnate values that are too small and make copmutation slow
        weights_mat(weights_mat<0) = 0;
		weights_mat=weights_mat./repmat(sum(weights_mat,1),mat_size,1);
		weights_mat=weights_mat./repmat(sum(weights_mat,2),1,mat_size);
		weights_mat = (weights_mat + weights_mat')/2;
    end
end

if (wid <= 0)
    errordlg('Non positive width');
    uiwait;
    return
end

% Score obtained by multiplying the distance matrix times the weight function

mismatch = dist_matrix * weights_mat;


[val,mn] = min(mismatch,[],2);
main_diag = diag(mismatch,0);
sort_score = 1:mat_size;
mx=max(val);
% f=find(main_diag > prctile(main_diag,perc));
% sort_score(f)=mn(f) - 0.1*sign((mat_size/2-mn(f))).*val(f)/(mx);
sort_score=mn - 0.1*sign((mat_size/2-mn)).*val/(mx);

% output options
if ((nargin>4) & (flag==1))
    figure; 
    imagesc(weights_mat);
end

% tic
% sort_score = assignment(mismatch);
% toc
% trace(mismatch)

% Sorting the matrix
[~, sorted_ind] = sort(sort_score);

% set(findobj('tag','data_scores_2'),'string',num2str(trace(dist_matrix(sorted_ind,sorted_ind)*weights_mat)));