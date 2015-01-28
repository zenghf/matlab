%% define spin system
nSpin = 3;
csN = 0;  % I: 15N; S: 1H 
csH1 = 0;
csH2 = 0;
J = 87;
JNH1 = J;
JNH2 = J;
ts = linspace(0,1,10);
B = 750 / 42.576;
coupling = 'weak';

csList = [csN, csH1, csH2];
JMat = [ 0,      JNH1,   JNH2;
        JNH1,    0,      0;
        JNH2,    0,      0];

spinSys = setSpinSys(nSpin, csList, zeros(1,nSpin), JMat);
NX = spinOperator('xee', spinSys);
NY = spinOperator('yee', spinSys);
HX = spinOperator('exe', spinSys) + spinOperator('eex', spinSys);
HY = spinOperator('eye', spinSys) + spinOperator('eey', spinSys);

op0 = spinOperator('zee', spinSys);% + spinOperator('eze', spinSys) + spinOperator('eez', spinSys);
%tau = 0.25 / J;

% spins are : [15N, 1H1, 1H2]

% [-- tau -- pi pulse on all spins -- tau --]
evolveJ = @(op, tau) evolve(pulse(evolve(op, spinSys, tau, B, coupling), [1 2 3], pi, 0), spinSys, tau, B, coupling); 

show = [];

tauList = linspace(0, 0.25, 100);
I3x = zeros(size(tauList));

for k = 1:length(tauList)
    tau = tauList(k) / J;

    op1 = pulse(op0, [1], pi/2, pi/2);
    op2 = evolveJ(op1, tau);
    op3 = pulse(op2, [1 2 3], pi/2, 0);
    op4 = evolveJ(op3, tau);
    
%     spinDecomposition(op0,'xyz',show);
%     spinDecomposition(op1,'xyz',show);
%     spinDecomposition(op2,'xyz',show);
%     spinDecomposition(op3,'xyz',show);    
%     ac = spinDecomposition(op4,'xyz',show);
    I3x(k) = getOpCoeff(op4, spinSys, 'eex');
end
plt(tauList, real(I3x));
