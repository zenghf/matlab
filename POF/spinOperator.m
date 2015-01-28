function [op varargout]=spinOperator(varargin)
% generate the matrix for the product operater
% input: the last parameter can be total number of spins
%        of a spin system defined by setSpinSys
%       the operator for each spin can be
%           'xyzpme', represent Ix, Iy, Iz, I+, I-, E
%      e.g.: spinOperator('xyp',4)
%            spinOperator(1,'p',3,'m',4)
% output 
%  op: the matrix of the product operator
%  the 2nd output: Notation of the product operator, 
%      e.g.: 4I1xI2yI3p, 2I1pI3m
narginchk(2,inf);
if isa(varargin{end}, 'struct')
    spinSys = varargin{end};
    nSpin = spinSys.nSpin;
else
    nSpin = varargin{end};
    spinSys = setSpinSys(nSpin);
end
orders = repmat('e',1,nSpin);
if isa(varargin{1},'double')    
    for k = 1:(nargin-1)/2
        orders(varargin{2*k-1}) = varargin{2*k};
    end
elseif isa(varargin{1},'char')
    orders(1:length(varargin{1})) = varargin{1};
end

Ix = [0 1/2;1/2 0];
Iy = [0 1i/2;-1i/2 0];
Iz = [1/2 0;0 -1/2];
Ip = [0 0;1 0]/sqrt(2);
Im = [0 1;0 0]/sqrt(2);
Ie = eye(2);

op = 1;
for k = 1:nSpin
    op = kron(op, eval(['I',orders(k)]));
end

opNotation = repmat(' ',1,3*nSpin);
n = -1;
for k = 1:nSpin            
    if orders(k) ~= 'e'
        opNotation(3*k-2:3*k) = [spinSys.spinNameList{k},orders(k)];
        n = n + 1;
    end
end
opNotation = [int2str(2^n),strrep(opNotation,' ', '')];
op = op * 2^n;
if nargout == 2
    varargout{1} = opNotation;
end