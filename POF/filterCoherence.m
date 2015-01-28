function op = filterCoherence(op0, coh, show)
% keep only some coherence order of a product operator
% input
%   op0 : the operator to filter
%   coh:  the coherence orders to filter
%   show: (optional), if it is not empty, then show decomposition result
% output
%   op  : the matrix of the product operator after filtration

% spinSys = setSpinSys(3);
% op0 = spinOperator('x',3) + spinOperator('xz',3);
% coh = [1 -1];
nargin
if nargin < 3 
    show = [];
end
allComp = spinDecomposition(op0, 'pm', show);
op = zeros(size(op0));
for k = 1:length(allComp)    
    ord = sum((allComp(k).opNotation == 'p') - (allComp(k).opNotation == 'm'));
    if any(ord == coh)
        op = op + allComp(k).coeff * allComp(k).opMat;
    end
end
spinDecomposition(op, 'pm', show);
end
