function [allComp str] = spinDecomposition(op0, decompMode, show)
% calculate all the compoents of a product operator
% input
%   op0 : the operator to decomposite
%   decompMode: (optional), value can be 'xyz' or 'pm'
%               'xyz': in the mode of Ix, Iy and Iz
%               'pm':  in the mode of I+, I- and Iz
%   show: (optional), if it is not empty, then show decomposition result
% output
%   allComp: all the components in a structure array
%       .opNoation  : the notation of the basic operator
%       .opMat      : the matrix of the basic operator
%       .norm       : the normalization coefficient of the basic operator
%       .coeff      : the amount of the basic operator
%op0 = spinOperator('zz',spinSys) + 0.5*spinOperator('xx',spinSys);
nSpin = round(log(size(op0,1)) / log(2));
spinSys = setSpinSys(nSpin, zeros(1,nSpin), zeros(1,nSpin), zeros(nSpin));
allComp = struct('opNotation',cell(1,4^nSpin),'opMat',[],'norm',[],'coeff',[]);
o1 = 'x'; o2 = 'y';
if nargin>=2 && strcmp(decompMode,'pm');
    o1 = 'p'; o2 = 'm';
end
for k = 0:4^nSpin-1
    digCode = dec2base(k,4,nSpin);
    opCode = strrep(strrep(strrep(strrep(digCode,'0','e'),'1',o1),'2',o2),'3','z');
    [op opNotation]=spinOperator(opCode, spinSys);
    allComp(k+1).opNotation = opNotation;
    allComp(k+1).opMat = op;
    allComp(k+1).norm = sum(op(:).^2);
    allComp(k+1).coeff = sum(sum(op0.*op)) / allComp(k+1).norm;
end
allCoeff = [allComp.coeff];
[coeff ind] = sort(abs(allCoeff),'descend');
nComp = sum(abs(coeff)>eps*10);
allComp = allComp(ind(1:nComp));

if nargin >= 3 && ~isempty(show)
    str = [];
    for k = 1:nComp
        str = [str, num2str(allComp(k).coeff, '%6.4f'), ' ', allComp(k).opNotation, ' + '];
    end
    if length(str) > 2 
        str(end-1:end) = [];
    end
    disp(str);
end
end