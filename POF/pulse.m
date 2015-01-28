function op = pulse(op0, spins, alpha, phi)
% calculate the product operator after applying a pulse
% input
%   op0     : the initial operator 
%   spins   : the spins to appy pulse
%   alpha   : flip angle / rad, 90 degree pulse : pi/2
%   phi     : phase / rad, x: 0; y: pi/2; -x: pi

% output
%   op      : the matrix of the product operator after the pulse



%   psi     : offset / rad = offest/Hz * pulse time/s * 2 * pi

% op0 = spinOperator('z',1);
% spins = [1];
% alpha = pi*1/2; phi = pi*1/2; psi = pi*1/4; 


nSpin = round(log(size(op0,1)) / log(2));
tmp = 0;psi=0;
for k = 1:length(spins)
    n = spins(k);
    tmp = tmp + 1i * (alpha*cos(psi)*...
        (spinOperator(n,'x',nSpin)*cos(phi) + spinOperator(n,'y',nSpin)*sin(phi))...
        +sin(psi)*spinOperator(n,'z',nSpin));
end
rotMat = expm(tmp);
op = rotMat * op0 * rotMat';



%a = spinDecomposition(op, 'xyz', 'show');