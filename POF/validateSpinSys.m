function validateSpinSys(spinSys)
% validate the spin system
% check the chemical shift list, chemical shift anisotropy list
% and J matrix

nSpin = spinSys.nSpin;
csList = spinSys.csList;
csaList = spinSys.csaList;
JMat = spinSys.JMat;

if length(csList) ~= nSpin
    error('chemical shift list dimension does not match');
end
if any(imag(csList))
    error('chemical shift must be real')
end


if length(csaList) ~= nSpin
    error('chemical shift anisotropy list dimension does not match');
end
if any(imag(csaList))
    error('chemical shift anisotropy must be real')
end

if length(size(JMat)) ~= 2 || any(size(JMat) ~=nSpin)
    error('J matrix dimension does not match');
end
if any(imag(JMat(:)))
    error('J matrix must be real')
end

if any(any(JMat - JMat'))
    error('J matrix must be symmetric')
end
if any(diag(JMat))
    error('diagonal of J matrix must be zero ')
end
