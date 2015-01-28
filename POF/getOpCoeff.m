function coeff = getOpCoeff(op, spinSys, opCode)
% decomposition of operator op, calculate the coefficient of opCode 
% input
%   op : the operator to decomposite
%   spinSys : spin system
%   opCode : string notation of the component either in 'xyz' mode or 'pm' mode
% output
%   coeff      : the amount of the basic operator
% coeff = getOpCoeff(op, spinSys, 'zex'); 
%      get the coeffecient of 2I1zI3x of 3-spin system.
opMat = spinOperator(opCode, spinSys);
coeff = sum(sum(opMat .* op)) / sum(opMat(:).^2);
end