%   I & S : proton from parahydrogen,chemical shift can be the same or different
%% define spin system
nSpin = 2;
csI = -1; csS = 1; 
JIS = 10; 
ts = linspace(0,0.1,50);
B = 10;
coupling = 'weak';
termInterest = {'2I1xI2x','2I1yI2y','2I1zI2z','1I1x','1I1y','1I1z','1I2x','1I2y','1I2z'};

csList = [csI, csS];
JMat = [0,  JIS;
      JIS,  0;];

spinSys = setSpinSys(nSpin, csList, zeros(1,nSpin), JMat);

op0 = spinOperator('ex',spinSys) - 0*spinOperator('ze',spinSys); 
spinDecomposition(op0,'xyz','');
spinDecomposition(op0,'pm','');

res = zeros(length(ts),4^nSpin);
notations = cell(1,4^nSpin);
kNotation = 1;
for k = 1:length(ts)
    t = ts(k);
    op = evolve(op0, spinSys, t, B, coupling);
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
res = real(res(:,sortIndex));
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

%%{
flagInterest = true(1, length(notations));
for k = 1:length(notations)
    if any(strcmp(termInterest,notations{k}))
        flagInterest(k) = true;
    end
end
res = res(:,flagInterest);
notations = notations(flagInterest);

plt(ts,real(res),'traceID',notations);