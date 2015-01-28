% spin system: ABCD, 
%   R & T : proton from ligand, chemical shift can be the same or different
% ref: J. Chem. Phys. 131 194505(2009) 
% A theoretical basis for spontaneous polarization transfer in non-hydrogenative parahydrogen-induced polarization.
%% define spin system
nSpin = 2;
csR = 0; 
csT = 0;
JRT = 0;
ts = linspace(0,1,10);
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

coeffs = [-1 -1/2 0 0 0];
a0 = coeffs(1); a = coeffs(2); b = coeffs(3); c = coeffs(4); d = coeffs(5);
% ref page 194505-4

op0 = a0 * spinOperator('zz',spinSys)/2 +...
      a  *(spinOperator('xx',spinSys) + spinOperator('yy',spinSys)) +...
      b  *(spinOperator('yx',spinSys) - spinOperator('xy',spinSys)) +...
      c  * spinOperator('ez',spinSys) +...
      d  * spinOperator('ze',spinSys);
spinDecomposition(op0,'xyz','show');
spinDecomposition(op0,'pm','show');

op1 = pulse(op0,[1 ],pi/6,0);
spinDecomposition(op1,'xyz','show');
spinDecomposition(op1,'pm','show');