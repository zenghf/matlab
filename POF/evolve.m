function op = evolve(op0, spinSys, t, B, coupling)
% evolution of product operator under chemical shift and J coupling
% input
%   op0     : the initial product operator
%   spinSys : the spin system with defined chemical shift list and J
%             coupling matrix
%   t       : the time of the evolution
%   B       : the static magnetic field strength in T, default is 0
%   coupling: 'strong' or 'week', default is 'weak' coupling
% output
%   op      : the product operator after evolution
narginchk(3,5);
if nargin < 4
    B = 0;
end
if nargin < 5
    coupling = 'weak';
end

% spinSys = setSpinSys(2, [1 1],[0 0],[0 5;5 0]);
% op0 = spinOperator('xx',spinSys);
% t = 1;
% B = 1;
% coupling = 'strong';

nSpin = spinSys.nSpin;
validateSpinSys(spinSys);

HJ = 0;
strong = strcmp(coupling,'strong');
for i = 1:nSpin
    for j = 1:nSpin
        HJ = HJ + (i~=j) * spinSys.JMat(i,j)* pi/2 * ...
            (        spinOperator(i,'z',j,'z',spinSys)+...
            strong * spinOperator(i,'x',j,'x',spinSys)+...
            strong * spinOperator(i,'y',j,'y',spinSys)  );
    end
end
HCS = 0;
for k = 1:nSpin
    HCS = HCS - spinSys.csList(k) * gmr('1H',B) * 2e-6 * pi * ...
          spinOperator(k,'z',spinSys);
end

H = HJ + HCS;
rotMat = expm(-1i*H*t);
op = rotMat * op0 * rotMat';
