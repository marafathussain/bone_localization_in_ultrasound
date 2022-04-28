%%
function [best_group_index,eigValues] = cluster(bone,K)
y = length(bone);
bone = [bone; (1:y)];
%figure,plot(bone(1,:), bone(2,:),'r+'), title('Original Data Points');
A = zeros(y);

% Estimating Sigma
sigma = zeros(1,y);
for m = 1:y
    mi = m - K; if mi < 1; mi = 1;end
    ma = m + K; if ma > y; ma = y;end
    
    dist = zeros(1,ma-mi+1);
    z = 1;
    for g = mi:ma
        dist(z) = sqrt((bone(1,m) - bone(1,g))^2 + (bone(2,m) - bone(2,g))^2);
        z = z + 1;
    end
    dist = sort(dist);
    dist = dist(dist > 0);
%     if length(dist) < 7
%         p = 1;
%     end
    sigma(m) = dist(K);
    clear dist;
end
        
for i = 1:y
    for j = 1:y        
        A(i,j) = exp(-((bone(1,i) - bone(1,j))^2 + (bone(2,i) - bone(2,j))^2)/(1*sigma(i)*sigma(j))); 
    end
end  

degs = sum(A, 2);
degs(degs == 0) = eps;
D = sparse(1:size(A, 1), 1:size(A, 2), degs);

%L1 = D - A;  %%% I used following line
L1 = A;
Di = spdiags(1./(degs.^0.5), 0, size(D, 1), size(D, 2));
NL = Di * L1 * Di;    

[~,eigValues] = eig(NL);
eigValues = sum(eigValues,2);

group_num = [3 4 5 6 7];
[~,best_group_index,~,~] = cluster_rotate(A,group_num,0,2);
