function spinSys = setSpinSys(spinList, csList, csaList, JMat)
% create a spin system
% input
%   spinList    : cell matrix of the spin names
%                 integer n: create a spin system contains n 1H
%   csList      : chemical shift list in ppm, default is 0
%   csaList     : chemical shift anisotropy list in ppm, default is 0
%   JMat        : J matrix in Hz, default is 0
% output
%   spinSys
%       .spinList, .csList, .csaList, .JMat, .nSpin, .spinNumberList, .spinNameList
if nargin == 0
    spinList = {'1H','1H'};
end
if nargin >= 1
    if isa(spinList, 'double')
        tmp  = spinList; spinList = cell(1, tmp);
        for k = 1:tmp
            spinList{k} = '1H';
        end
    end
end
nSpin = length(spinList);
if nargin<2 || isempty(csList)
    csList = zeros(1,nSpin);
    %disp('set all chemical shift to be 0 ppm');
end
if nargin<3 || isempty(csaList)
    csaList = zeros(1,nSpin);
    %disp('set all chemical shift anisotropy to be 0 ppm');
end
if nargin<4 || isempty(JMat)
    JMat = zeros(nSpin);
    %disp('set all J coupling to be 0 Hz');
end
spinNumberList = arrayfun(@(x) gmr(x,'spin'), spinList);
spinNameList = cell(1,nSpin);

singleLetterNameList = 'ISRTPQ';
spinNameList{1} = 'I1';
spinNameMat = zeros(nSpin, length(singleLetterNameList));
for k = 1:nSpin
    cDigit = sum(strcmp(spinList{k},spinList(1:k-1))) + 1;
    nLetter = sum(any(spinNameMat>0));
    if cDigit == 1
        cLetter = nLetter + 1;
    else
        cLetter = find(strcmp(spinList{k}, spinList(spinNameMat(1,1:nLetter))),1);
    end
    spinNameMat(cDigit, cLetter) = k;
    spinNameList{k} = [singleLetterNameList(cLetter),int2str(cDigit)];
end

parReturn = {'spinList','csList','csaList','JMat','nSpin','spinNumberList','spinNameList'};
for k = parReturn
    spinSys.(k{:}) = eval(k{:});
end
validateSpinSys(spinSys);