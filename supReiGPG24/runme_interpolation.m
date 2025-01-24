%% RUNME_INTERPOLATION

% Script file to check numerically results in the interpolation paper.

%
% This file is part of the archive Code and Results for Numerical 
% Experiments in "..."
% Copyright (c) 2024 Sean Reiter, Steffen W. R. Werner, and ...
% All rights reserved.
% License: BSD 2-Clause license (see COPYING)
%

clc;
clear all;
close all;

% Get and set all paths
[rootpath, filename, ~] = fileparts(mfilename('fullpath'));
loadname            = [rootpath filesep() ...
    'data' filesep() filename];
savename            = [rootpath filesep() ...
    'results' filesep() filename];

% Add paths to drivers and data
addpath([rootpath, '/drivers'])
addpath([rootpath, '/data'])

% Write .log file, put in `out' folder
if exist([savename '.log'], 'file') == 2
    delete([savename '.log']);
end
outname = [savename '.log']';

diary(outname)
diary on; 

fprintf(1, ['SCRIPT: ' upper(filename) '\n']);
fprintf(1, ['========' repmat('=', 1, length(filename)) '\n']);
fprintf(1, '\n');


%% 
% Build toy problem.

% Full-order model.
n = 8;  p = 2;  m = 2;     % Specify state, input, output dimensions
A = - diag(10*rand(n, 1));
B = 5*rand(n, m);
C = 5*rand(p, n);
E = eye(n, n);

% Allocate quadratic outputs; build both representations.
M        = zeros(p, n^2);
MOutputs = zeros(n, n, p); % Contains p individual outputs
for i = 1:p
    tmpMatZ            = 5*rand(n, n);
    MOutputs(:, :, i) = tmpMatZ + tmpMatZ';
    M(i, :)           = reshape(MOutputs(:, :, i), 1, []);
end

% Order of reduction.
r = 4;  Er = eye(r, r);

% Generate random reduced model (complex conjugate poles).
poles       = -1*rand(r/2, 1) + 1i*rand(r/2, 1);
poles       = [poles; conj(poles)]; 
rightEigVec = rand(r, r/2) + 1i*rand(r, r/2);
rightEigVec = [rightEigVec, conj(rightEigVec)];
tmp         = rightEigVec\(eye(r, r));  
leftEigVec  = tmp.';                                % Transpose, not conjugate transpose.
% Ar          = diag(poles);
Ar          = rightEigVec*diag(poles)*leftEigVec.'; % Not in diagonalized form.
Br          = rand(r, m);
Cr          = rand(p, r);
Mr          = zeros(p, r^2);
MrOutputs   = zeros(r, r, p);
for i = 1:p
    tmpMatZ            = rand(r, r);
    MrOutputs(:, :, i) = tmpMatZ + tmpMatZ';
    Mr(i, :)           = reshape(MrOutputs(:, :, i), 1, []);
end

% Interpolation points and tangent directions from poles and residues.
rightDir    = Br.'*leftEigVec;   % Right tangent directions; store as m x 1 columns
leftDirLin  = Cr*rightEigVec;    % Left tangent directions (linear output); store as p x 1 columns
leftDirQuad = zeros(p, 1, r, r); % Left tangent directions (quadratic output); store as p x 1 columns
for i = 1:r
    for j = 1:r
        leftDirQuad(:, :, i, j) = Mr*kron(rightEigVec(:, i), rightEigVec(:, j));
    end
end

%%
% Numerical checks on Kronecker product and vectorization identities.

Im = eye(m, m);
% Check Lemma: H2(s_1, s_2)*kron(u, Im) =  H2(s_2, s_1)*kron(Im, u).
% (Asrises in construction of rows of W.)
u = rand(n, 1);
for i = 1:r
    for l = 1:r
        lhs1 = M*kron(((-poles(l)*E - A)\B), ((-poles(i)*E - A)\B))*kron(rightDir(:, l), Im);
        lhs2 = M*kron(((-poles(i)*E - A)\B), ((-poles(l)*E - A)\B))*kron(Im, rightDir(:, l));
            fprintf(1, 'Relative error in symmetry of H2 for pair (%d, %d): %d\n', i, l, ...
                norm((lhs1 - lhs2)/lhs1, 2))
            fprintf(1, '---------------------------------------------------------------------------------\n')
    end
end

% Check Lemma: \frac{\partial}{\partial s_1}H2(s, z)*kron(u, Im) =  \frac{\partial}{\partial s_2}H2(z, s)*kron(Im, u).
% (Asrises in construction of rows of W.)
u = rand(n, 1);
for i = 1:r
    for l = 1:r
        % Note: This formulation is fine because E = I.
        lhs1 = M*kron(-(-poles(l)*E - A)\((-poles(l)*E - A)\B), ((-poles(i)*E - A)\B))*kron(rightDir(:, l), Im);
        lhs2 = M*kron(((-poles(i)*E - A)\B), -(-poles(l)*E - A)\((-poles(l)*E - A)\B))*kron(Im, rightDir(:, l));
            fprintf(1, 'Relative error in symmetry of H2s derivatives for pair (%d, %d): %d\n', i, l, ...
                norm((lhs1 - lhs2)/lhs1, 2))
            fprintf(1, '---------------------------------------------------------------------------------\n')
    end
end


% % Check vectorization and Kronecker product identity: 
% %   vec(B*X*A).' = vec(X).'*kron(A, B)
% tmp1 = zeros(p, n);
% ell  = 1;
% i    = 2;
% for k = 1:p 
%     tmp1(k, :) = rightDir(:, ell)'*B.'*((-poles(ell)*E.'-A.')\(MOutputs(:, :, k)*((-conj(poles(i))*E-A)\E)));
% end
% tmp2 = M*kron(((-conj(poles(i))*E-A)\E), ((-(poles(ell))*E-A)\(B*conj(rightDir(:, ell)))));
% fprintf(1, 'Relative error in vectorization and Kronecker product identity: vec(B*X*A).T = vec(X).T*kron(A, B): %d\n', ...
%     norm((tmp1 - tmp2)/tmp1, 2))
% fprintf(1, '--------------------------------------------------------------------------------------------------\n')


%%
% Check formulae for interpolation bases.

% 1. From matrix equations.
%    Solve in the basis of Ar's eigenvectors; Ar already diagonalized
%       Hit with leftEigVec (T) from right -> Br'*leftEigVec = rightDir
%       Note: (Ar.'*leftEigVec) = leftEigVec*diag(poles); we want to solve
%       for, e.g., XMatEq*leftEigVec, so just plug in diag(poles) for second
%       argument.
PrMatEq = lyap(Ar, diag(poles), Br*rightDir); % Ar*(PrMatEq*leftEigVec) + (PrMatEq*leftEigVec)*diag(poles) + Br*(Br*leftEigVec).' = 0
XMatEq  = lyap(A,  diag(poles), B*rightDir);  % A*(XMatEq*leftEigVec)   + (X*leftEigVec)*diag(poles)       + B*(Br*leftEigVec).'  = 0


%       Hit with rightEigVec (S) from right -> Cr*S = leftDirLin
rhsQ = Cr.'*leftDirLin; 
rhsZ = -C.'*leftDirLin; 

% Below, PrMatEq = (PrMatEq*leftEigVec) and (X*leftEigVec)
% That c.o.b is implicitly hidden, so need to pull out its inverse
for i = 1:p
    rhsQ  = rhsQ + 2*MrOutputs(:, :, i)*(PrMatEq)*(rightEigVec.'*MrOutputs(:, :, i)*rightEigVec);
    rhsZ  = rhsZ - 2*MOutputs(:, :, i)*(XMatEq)*(rightEigVec.'*MrOutputs(:, :, i)*rightEigVec);
end
%       Note: (Ar*rightEigVec) = rightEigVec*diag(poles); we want to 
%       solve for, e.g., ZMatEq*rightEigVec, so just plug in diag(poles) 
%       for second argument.

% Ar.'*(QrMatEq*RightEigVec) + (QrMatEq*RightEigVec)*diag(poles) + rhsQ = 0
QrMatEq = lyap(Ar.', diag(poles), rhsQ);  
% A.'*(ZMatEq*RightEigVec)   + (ZMatEq*RightEigVec)*diag(poles)  - rhsZ = 0
ZMatEq  = lyap(A.',  diag(poles), rhsZ);  

% 2. From rational Kyrlov subspace formulation.
XKrylov = zeros(n, r);  PrKrylov = zeros(r, r);
ZKrylov = zeros(n, r);  QrKrylov = zeros(r, r);
for i = 1:r
    XKrylov(:, i)  = (-poles(i)*E - A)\(B*rightDir(:, i));
    PrKrylov(:, i) = (-poles(i)*Er - Ar)\(Br*rightDir(:, i));
end

for i = 1:r
    tmpSumZ  = C.'*leftDirLin(:, i);
    tmpSumQr = Cr.'*leftDirLin(:, i);
    for l = 1:r
        tmpMatZ  = zeros(n, p);
        tmpMatQr = zeros(r, p);
        for j = 1:p
            tmpMatZ(:, j)  = MOutputs(:, :, j)*XKrylov(:, l);
            tmpMatQr(:, j) = MrOutputs(:, :, j)*PrKrylov(:, l);
        end
        tmpSumZ  = tmpSumZ + 2*tmpMatZ*leftDirQuad(:, :, i, l);
        tmpSumQr = tmpSumQr + 2*tmpMatQr*leftDirQuad(:, :, i, l);
    end
    ZKrylov(:, i)  = -(-poles(i)*E-A.')\tmpSumZ;
    QrKrylov(:, i) = (-poles(i)*Er-Ar.')\tmpSumQr;
end

fprintf(1, 'Relative error in equivalent formulations of X : %d\n', ...
    norm(XMatEq  - XKrylov, 'fro')/norm(XMatEq, 'fro'))
fprintf(1, '--------------------------------------------------------------\n')
fprintf(1, 'Relative error in equivalent formulations of Pr: %d\n', ...
    norm(PrMatEq - PrKrylov, 'fro')/norm(PrMatEq, 'fro'))
fprintf(1, '--------------------------------------------------------------\n')
fprintf(1, 'Relative error in equivalent formulations of Z : %d\n', ...
    norm(ZMatEq  - ZKrylov, 'fro')/norm(ZMatEq, 'fro'))
fprintf(1, '--------------------------------------------------------------\n')
fprintf(1, 'Relative error in equivalent formulations of Qr: %d\n', ...
    norm(QrMatEq - QrKrylov, 'fro')/norm(QrMatEq, 'fro'))
fprintf(1, '--------------------------------------------------------------\n')

% 3. Alternative representation for rows of Z and Qr. 
%    Note: This is less efficient due to the way the solves have to be
%          done, but for small problems like this it is fine.
%    We expect ZKrylov.' = ZKrylovTrans

ZKrylovTrans  = zeros(r, n); 
QrKrylovTrans = zeros(r, r); 
for i = 1:r
    tmpSolveZ  = ((-poles(i)*E - A)\E);
    tmpSolveQr = ((-poles(i)*Er - Ar)\Er);
    tmpSumZ    = -leftDirLin(:, i).'*C*(tmpSolveZ);
    tmpSumQr   = leftDirLin(:, i).'*Cr*(tmpSolveQr);
    for l = 1:r
        tmpSumZ  = tmpSumZ - leftDirQuad(:, :, i, l).'*M*kron(tmpSolveZ, ((-poles(l)*E - A)\(B*rightDir(:, l)))) ...
                    - leftDirQuad(:, :, l, i).'*M*kron(((-poles(l)*E - A)\(B*rightDir(:, l))), tmpSolveZ);
        tmpSumQr = tmpSumQr + leftDirQuad(:, :, i, l).'*Mr*kron(tmpSolveQr, ((-poles(l)*Er - Ar)\(Br*rightDir(:, l)))) ...
                    + leftDirQuad(:, :, l, i).'*Mr*kron(((-poles(l)*Er - Ar)\(Br*rightDir(:, l))), tmpSolveQr);
    end
    ZKrylovTrans(i, :)  = tmpSumZ;
    QrKrylovTrans(i, :) = tmpSumQr;
end

fprintf(1, 'Relative error in Krylov-based formulations of Z  and Z.T : %d\n', ...
    norm(ZKrylovTrans.'  - ZKrylov, 'fro')/norm(ZKrylov, 'fro'))
fprintf(1, '--------------------------------------------------------------\n')
fprintf(1, 'Relative error in Krylov-based formulations of Qr and Qr.T: %d\n', ...
    norm(QrKrylovTrans.' - QrKrylov, 'fro')/norm(QrKrylov, 'fro'))
fprintf(1, '--------------------------------------------------------------\n')


%%
% Do the matrices enforce the expected interpolation conditions?
% Bases are computed from the previous section.

W = ((ZKrylov.'*XKrylov)\ZKrylov.');
V = XKrylov;
Ar = W*A*V; 
Br = W*B;    
Cr = C*V;
Mr = M*(kron(V, V));
Ir = eye(r, r);

fprintf(1, 'TESTING INTERPOLATION CONDITIONS (Krylov based construction).\n')
fprintf(1, '-----------------------------------------------------------------\n')
% Define transfer functions
HLin     = @(s)      C* ((s*E - A)\B);
HLinRed  = @(s)      Cr*((s*Ir- Ar)\Br);
HQuad    = @(s1, s2) M* kron(((s1*E - A)\B),    ((s2*E - A)\B));
HQuadRed = @(s1, s2) Mr*kron(((s1*Ir - Ar)\Br), ((s2*Ir - Ar)\Br));

fprintf(1, '1. Linear, right tangential interpolation conditions.\n')
fprintf(1, '-----------------------------------------------------------------\n')
for k = 1:r
    fprintf(1, 'Relative error in condition %d: %d\n', k, ...
        norm(HLin(-poles(k))*rightDir(:, k) - HLinRed(-poles(k))*rightDir(:, k), 2)/norm(HLin(-poles(k))*rightDir(:, k), 2))
    fprintf(1, '-----------------------------------------------------------------\n')
end

fprintf(1, '2. Quadratic, right tangential interpolation conditions.\n')
 fprintf(1, '-----------------------------------------------------------------\n')
for j = 1:r
    for k = 1:r
        fprintf(1, 'Relative error in condition (%d, %d): %d\n', j, k, ...
            norm(HQuad(-poles(j), -poles(k))*kron(rightDir(:, j), rightDir(:, k)) ...
            - HQuadRed(-poles(j), -poles(k))*kron(rightDir(:, j), rightDir(:, k)), 2)/ ...
            norm(HQuad(-poles(j), -poles(k))*kron(rightDir(:, j), rightDir(:, k)), 2))
        fprintf(1, '-----------------------------------------------------------------\n')
    end
end

Im = eye(m, m);
fprintf(1, '3. Mixed, left tangential interpolation conditions.\n')
fprintf(1, '-----------------------------------------------------------------\n')
for k = 1:r
    % Quadratic component of mixed term
    quadTermFO = zeros(1, m);
    quadTermRO = zeros(1, m);
    for j = 1:r
        quadTermFO = quadTermFO + leftDirQuad(:, :, k, j).'*HQuad(-poles(k), -poles(j))*kron(Im, rightDir(:, j)) + ...
            leftDirQuad(:, :, j, k).'*HQuad(-poles(j), -poles(k))*kron(rightDir(:, j), Im);
        quadTermRO = quadTermRO + leftDirQuad(:, :, k, j).'*HQuadRed(-poles(k), -poles(j))*kron(Im, rightDir(:, j)) + ...
            leftDirQuad(:, :, j, k).'*HQuadRed(-poles(j), -poles(k))*kron(rightDir(:, j), Im);
    end
    % Linear component of mixed term
    mixedTermFO = quadTermFO + leftDirLin(:, k).'*HLin(-poles(k));
    mixedTermRO = quadTermRO + leftDirLin(:, k).'*HLinRed(-poles(k));
    fprintf(1, 'Relative error in condition %d: %d\n', k, norm(mixedTermFO - mixedTermRO, 2)/norm(mixedTermFO, 2))
    fprintf(1, '-----------------------------------------------------------------\n')
end

HLinDeriv         = @(s)      -C* ((s*E - A)\((s*E - A)\B));
HLinRedDeriv      = @(s)      -Cr*((s*Ir- Ar)\((s*Ir- Ar)\Br));
HQuadPartials1    = @(s1, s2) -M* kron((s1*E - A)\((s1*E - A)\B),    (s2*E - A)\B);
HQuadRedPartials1 = @(s1, s2) -Mr*kron((s1*Ir- Ar)\((s1*Ir- Ar)\Br), (s2*Ir - Ar)\Br);
HQuadPartials2    = @(s1, s2) -M* kron((s1*E - A)\B,                 (s2*E - A)\((s2*E - A)\B));
HQuadRedPartials2 = @(s1, s2) -Mr*kron(((s1*Ir - Ar)\Br),            ((s2*Ir - Ar)\((s2*Ir - Ar)\Br)));

fprintf(1, '4. Mixed, bi-tangential Hermite interpolation conditions.\n')
fprintf(1, '-----------------------------------------------------------------\n')
for k = 1:r
    % Quadratic component of mixed term
    quadTermFO = 0;
    quadTermRO = 0;
    for j = 1:r
        quadTermFO = quadTermFO + leftDirQuad(:, :, k, j).'*HQuadPartials1(-poles(k), -poles(j))*kron(rightDir(:, k), rightDir(:, j)) + ...
            leftDirQuad(:, :, j, k).'*HQuadPartials2(-poles(j), -poles(k))*kron(rightDir(:, j), rightDir(:, k));
        quadTermRO = quadTermRO + leftDirQuad(:, :, k, j).'*HQuadRedPartials1(-poles(k), -poles(j))*kron(rightDir(:, k), rightDir(:, j)) + ...
            leftDirQuad(:, :, j, k).'*HQuadRedPartials2(-poles(j), -poles(k))*kron(rightDir(:, j), rightDir(:, k));
    end
    % Linear component of mixed term
    mixedTermFO = quadTermFO + leftDirLin(:, k).'*HLinDeriv(-poles(k))*rightDir(:, k);
    mixedTermRO = quadTermRO + leftDirLin(:, k).'*HLinRedDeriv(-poles(k))*rightDir(:, k);
    fprintf(1, 'Relative error in condition %d: %d\n', k, norm(mixedTermFO - mixedTermRO, 2)/norm(mixedTermFO, 2))
    fprintf(1, '-----------------------------------------------------------------\n')
end


%%
% Solve without the change of basis. 
XMatEq  = lyap(A,  Ar.', B*Br.');  

rhsZ = -C.'*Cr; 
for i = 1:p
    rhsZ  = rhsZ - 2*MOutputs(:, :, i)*XMatEq*MrOutputs(:, :, i);
end
% A.'*(ZMatEq*RightEigVec)   + (ZMatEq*RightEigVec)*diag(poles)  - rhsZ = 0
ZMatEq  = lyap(A.', Ar, rhsZ);  

W = ((ZMatEq.'*XMatEq)\ZMatEq.');
V = XMatEq;
Ar = W*A*V; 
Br = W*B;    
Cr = C*V;
Mr = M*(kron(V, V));
Ir = eye(r, r);

fprintf(1, 'TESTING INTERPOLATION CONDITIONS (Sylvester based construction).\n')
fprintf(1, '-----------------------------------------------------------------\n')
% Define transfer functions
HLin     = @(s)      C* ((s*E - A)\B);
HLinRed  = @(s)      Cr*((s*Ir- Ar)\Br);
HQuad    = @(s1, s2) M* kron(((s1*E - A)\B),    ((s2*E - A)\B));
HQuadRed = @(s1, s2) Mr*kron(((s1*Ir - Ar)\Br), ((s2*Ir - Ar)\Br));

fprintf(1, '1. Linear, right tangential interpolation conditions.\n')
fprintf(1, '-----------------------------------------------------------------\n')
for k = 1:r
    fprintf(1, 'Relative error in condition %d: %d\n', k, ...
        norm(HLin(-poles(k))*rightDir(:, k) - HLinRed(-poles(k))*rightDir(:, k), 2)/norm(HLin(-poles(k))*rightDir(:, k), 2))
    fprintf(1, '-----------------------------------------------------------------\n')
end

fprintf(1, '2. Quadratic, right tangential interpolation conditions.\n')
 fprintf(1, '-----------------------------------------------------------------\n')
for j = 1:r
    for k = 1:r
        fprintf(1, 'Relative error in condition (%d, %d): %d\n', j, k, ...
            norm(HQuad(-poles(j), -poles(k))*kron(rightDir(:, j), rightDir(:, k)) ...
            - HQuadRed(-poles(j), -poles(k))*kron(rightDir(:, j), rightDir(:, k)), 2)/ ...
            norm(HQuad(-poles(j), -poles(k))*kron(rightDir(:, j), rightDir(:, k)), 2))
        fprintf(1, '-----------------------------------------------------------------\n')
    end
end

Im = eye(m, m);
fprintf(1, '3. Mixed, left tangential interpolation conditions.\n')
fprintf(1, '-----------------------------------------------------------------\n')
for k = 1:r
    % Quadratic component of mixed term
    quadTermFO = zeros(1, m);
    quadTermRO = zeros(1, m);
    for j = 1:r
        quadTermFO = quadTermFO + leftDirQuad(:, :, k, j).'*HQuad(-poles(k), -poles(j))*kron(Im, rightDir(:, j)) + ...
            leftDirQuad(:, :, j, k).'*HQuad(-poles(j), -poles(k))*kron(rightDir(:, j), Im);
        quadTermRO = quadTermRO + leftDirQuad(:, :, k, j).'*HQuadRed(-poles(k), -poles(j))*kron(Im, rightDir(:, j)) + ...
            leftDirQuad(:, :, j, k).'*HQuadRed(-poles(j), -poles(k))*kron(rightDir(:, j), Im);
    end
    % Linear component of mixed term
    mixedTermFO = quadTermFO + leftDirLin(:, k).'*HLin(-poles(k));
    mixedTermRO = quadTermRO + leftDirLin(:, k).'*HLinRed(-poles(k));
    fprintf(1, 'Relative error in condition %d: %d\n', k, norm(mixedTermFO - mixedTermRO, 2)/norm(mixedTermFO, 2))
    fprintf(1, '-----------------------------------------------------------------\n')
end

HLinDeriv         = @(s)      -C* ((s*E - A)\((s*E - A)\B));
HLinRedDeriv      = @(s)      -Cr*((s*Ir- Ar)\((s*Ir- Ar)\Br));
HQuadPartials1    = @(s1, s2) -M* kron((s1*E - A)\((s1*E - A)\B),    (s2*E - A)\B);
HQuadRedPartials1 = @(s1, s2) -Mr*kron((s1*Ir- Ar)\((s1*Ir- Ar)\Br), (s2*Ir - Ar)\Br);
HQuadPartials2    = @(s1, s2) -M* kron((s1*E - A)\B,                 (s2*E - A)\((s2*E - A)\B));
HQuadRedPartials2 = @(s1, s2) -Mr*kron(((s1*Ir - Ar)\Br),            ((s2*Ir - Ar)\((s2*Ir - Ar)\Br)));

fprintf(1, '4. Mixed, bi-tangential Hermite interpolation conditions.\n')
fprintf(1, '-----------------------------------------------------------------\n')
for k = 1:r
    % Quadratic component of mixed term
    quadTermFO = 0;
    quadTermRO = 0;
    for j = 1:r
        quadTermFO = quadTermFO + leftDirQuad(:, :, k, j).'*HQuadPartials1(-poles(k), -poles(j))*kron(rightDir(:, k), rightDir(:, j)) + ...
            leftDirQuad(:, :, j, k).'*HQuadPartials2(-poles(j), -poles(k))*kron(rightDir(:, j), rightDir(:, k));
        quadTermRO = quadTermRO + leftDirQuad(:, :, k, j).'*HQuadRedPartials1(-poles(k), -poles(j))*kron(rightDir(:, k), rightDir(:, j)) + ...
            leftDirQuad(:, :, j, k).'*HQuadRedPartials2(-poles(j), -poles(k))*kron(rightDir(:, j), rightDir(:, k));
    end
    % Linear component of mixed term
    mixedTermFO = quadTermFO + leftDirLin(:, k).'*HLinDeriv(-poles(k))*rightDir(:, k);
    mixedTermRO = quadTermRO + leftDirLin(:, k).'*HLinRedDeriv(-poles(k))*rightDir(:, k);
    fprintf(1, 'Relative error in condition %d: %d\n', k, norm(mixedTermFO - mixedTermRO, 2)/norm(mixedTermFO, 2))
    fprintf(1, '-----------------------------------------------------------------\n')
end


%% Finished script.
fprintf(1, 'FINISHED SCRIPT.\n');
fprintf(1, '================\n');
fprintf(1, '\n');

diary off