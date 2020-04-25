function Hmatrix=makeHmatrix(nfo)
% Calculates the hamiltonian.
depth=nfo.depth;                    % static lattice depth
k=nfo.k;                    % quasimomentum
theta=nfo.theta;            % modulation phase
alpha=nfo.alpha;            % depth of modulation
numStates=nfo.numStates;    % number of basis states

%%%%%%%%%%%%%% crystal momentum operator %%%%%%%%%%%%%%%%
C=0;
for ii=2:2:numStates
   C=[C; ii; 0]; 
end
C(end)=[];
D=zeros(numStates,1);
M1=gallery('tridiag',C,D,-C);

%%%%%%%%%%%%% kinetic energy operator %%%%%%%%%%%%%%
M2=0;
for ii=2:2:(numStates-1)
    M2=[M2 ii^2 ii^2];
end
M2=sparse(-diag(M2));

%%%%%%%%%%%%%%%%%% lattice operator %%%%%%%%%%%%%%%%%%%%
rowsL=[3:1:numStates];
colsL=[1:1:numStates-2];
valsL=[sqrt(2) ones(1,numStates-3)];
M3=sparse(rowsL,colsL,valsL,numStates,numStates)+...
   sparse(colsL,rowsL,valsL,numStates,numStates);

%%%%%%%%%%%%%% Put it all together %%%%%%%%%%%%%%%%%
Hmatrix=(k^2*speye(numStates)-2*1i*k*M1-M2)...
    -depth/4*M3*(1+alpha*sin(theta));
end

