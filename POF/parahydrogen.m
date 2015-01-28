% A pulse is applied to a spin state of hydrogen
% This script plot the product operator dependence the flip angle.
nSpin = 2;
Tp = zeros(4); Tm = zeros(4); T0 = zeros(4); S0 = zeros(4);
Tp(1,1) = 1; Tm(4,4) = 1; 
T0(2,3) = 1; T0(3,2) = 1; 
S0(2,3) = -1; S0(3,2) = 1;
T = Tp + T0 + Tm;
spinDecomposition(Tp,'xyz','show');
spinDecomposition(Tm,'xyz','show');
spinDecomposition(T0,'xyz','show');
spinDecomposition(S0,'xyz','show');

op1 = pulse(Tm, [1,2],pi/2,0);
spinDecomposition(op1,'xyz','show');
flipAngles = linspace(0,pi,100);
% this is the initial state
% op0 = T0, S0*1i, Tm, Tp, 
op0 = Tp * 1;

res = zeros(length(flipAngles),4^nSpin);
notations = cell(1,4^nSpin);
kNotation = 1;
for k = 1:length(flipAngles)
    alpha = flipAngles(k);
    op = pulse(op0,[1 2], alpha,0);
    comp = spinDecomposition(op,'xyz','');
    for n = 1:length(comp)
        tmp = find(strcmp(notations, comp(n).opNotation), 1);
        if ~isempty(tmp)
            res(k,tmp) = comp(n).coeff;
        else
            notations{kNotation} = comp(n).opNotation;
            res(k, kNotation) = comp(n).coeff;
            kNotation = kNotation + 1;
        end
    end
end
validTerm = any(res,1);
notations = notations(validTerm);
res = res(:,validTerm);
[notations,sortIndex] = sort(notations);
res = res(:,sortIndex);

% combine the same trend
%{
nTerms = sum(validTerm);
uniqueTerm = true(1,nTerms);
for k = 2:nTerms
    tmp = res(:,1:k-1) - repmat(res(:,k),1,k-1);
    kUnique = find(sum(abs(tmp)) < length(ts) * 0.001, 1);
    if ~isempty(kUnique)
        uniqueTerm(k) = false;
        %notations{kUnique} = [notations{kUnique},'&',notations{k}];
    end
end
notations = notations(uniqueTerm);
for k = 1:length(notations)
    tmp = notations{k};
    notations{k} = tmp(tmp~='I');
end
res = res(:,uniqueTerm);
%}

%{
flagInterest = false(1, length(notations));
for k = 1:length(notations)
    if any(strcmp(termInterest,notations{k}))
        flagInterest(k) = true;
    end
end
res = res(:,flagInterest);
notations = notations(flagInterest);
%}
res = real(res);
plt(flipAngles,res,'traceID',notations);