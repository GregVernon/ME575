function f = bezier(nodes,basisFun)
degree = size(nodes,1)-1;
nPts = size(basisFun,1);

f = zeros(nPts,1);
for ii = 1:degree+1
    f = f + basisFun(:,ii) .* nodes(ii,:);
end
end