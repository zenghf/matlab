% spin system: ABCD, 
%   R & T : proton from ligand, chemical shift can be the same or different
% ref: J. Chem. Phys. 131 194505(2009) 
% A theoretical basis for spontaneous polarization transfer in non-hydrogenative parahydrogen-induced polarization.
%% define spin system
nSpin = 2;
csR = 0; 
csT = 0;
JRT = 87;
ts = linspace(0,1,10);
B = 100e-4;
coupling = 'weak';
%{
termInterest = {'2I1xI2x','2I1yI2y','2I1zI2z','1I1x','1I1y','1I1z','1I2x','1I2y','1I2z',...
                '2I3xI4x','2I3yI4y','2I3zI4z','1I3x','1I3y','1I3z','1I4x','1I4y','1I4z'};
%}
csList = [csR, csT];
JMat = [ 0,      JRT;
        JRT,    0];

spinSys = setSpinSys(nSpin, csList, zeros(1,nSpin), JMat);

op0 = spinOperator('ez',spinSys) +...
      spinOperator('xe',spinSys);
spinDecomposition(op0,'xyz','show');
spinDecomposition(op0,'pm','show');

H = 1000 * spinOperator('xe',spinSys) + spinOperator('zz',spinSys) * pi * JRT;
spinDecomposition(H,'xyz','show');
spinDecomposition(H,'pm','show');
t = 0.05;
rotMat = expm(-1i*H*t);
op = rotMat * op0 * rotMat';
spinDecomposition(op, 'xyz', 'show');