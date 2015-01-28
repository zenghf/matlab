% spin system: ABCD, 
%   R & T : proton from ligand, chemical shift can be the same or different
% ref: J. Chem. Phys. 131 194505(2009) 
% A theoretical basis for spontaneous polarization transfer in non-hydrogenative parahydrogen-induced polarization.
%% define spin system
nSpin = 2;
csR = 10; 
csT = -1;
JRT = 7;
ts = linspace(0,1,1000);
B = 100e-4;
coupling = 'strong';
%{
termInterest = {'2I1xI2x','2I1yI2y','2I1zI2z','1I1x','1I1y','1I1z','1I2x','1I2y','1I2z',...
                '2I3xI4x','2I3yI4y','2I3zI4z','1I3x','1I3y','1I3z','1I4x','1I4y','1I4z'};
%}
csList = [csR, csT];
JMat = [ 0,      JRT;
        JRT,    0];

spinSys = setSpinSys(nSpin, csList, zeros(1,nSpin), JMat);

coeffs = rand(1,5);%[1 1 1 1 1];
a0 = coeffs(1); a = coeffs(2); b = coeffs(3); c = coeffs(4); d = coeffs(5);
% ref page 194505-4

op0 = a0 * spinOperator('zz',spinSys)/2 +...
      a  *(spinOperator('xx',spinSys) + spinOperator('yy',spinSys)) +...
      b  *(spinOperator('yx',spinSys) - spinOperator('xy',spinSys)) +...
      c  * spinOperator('ez',spinSys) +...
      d  * spinOperator('ze',spinSys);
spinDecomposition(op0,'xyz','show');
spinDecomposition(op0,'pm','show');



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

plt(ts,res,'traceID',notations);