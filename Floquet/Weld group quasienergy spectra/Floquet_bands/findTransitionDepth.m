function freq=findTransitionDepth(depth,n1,n2,k)
%FINDLATTICEDEPTH Summary of this function goes here
%   Detailed explanation goes here
depths=linspace(1,20,100);
paramsBANDS;

switch nargin
    case 1
        n1=1;
        n2=3;
        k=0;
end   

fr=25.18;
H1=makeHmatrix(k,depth);
[~,b1]=eig(full(H1));
lambda1=b1*ones(numStates,1); 
eng=abs(lambda1(n1)-lambda1(n2));

freq=eng*fr;
end
