function [hF,funcs] = computeTransitionOverlaps(A,maxBand)
%COMPUTETRANSITIONOVERLAPS Summary of this function goes here
%   Detailed explanation goes here

if nargin~=2
   A=5;
   maxBand=5;
end

numStates=35;
kpoints=501;
kvec=linspace(-1,1,kpoints);
% maxBand=5;

%% Make Hamiltonian

% First Derivative
C=0;
for ii=2:2:numStates
   C=[C; ii; 0]; 
end
C(end)=[];
D=zeros(numStates,1);
M1=gallery('tridiag',C,D,-C);

% Second derivative matrix
M2=0;
for ii=2:2:(numStates-1)
    M2=[M2 ii^2 ii^2];
end
M2=sparse(-diag(M2));

% Cosine Matrix
rowsL=[3:1:numStates];
colsL=[1:1:numStates-2];
valsL=[sqrt(2) ones(1,numStates-3)];
M3=sparse(rowsL,colsL,valsL,numStates,numStates)+...
   sparse(colsL,rowsL,valsL,numStates,numStates);

% Make the hamiltonian function
Hmat=@(k) (k^2*speye(numStates)-2*1i*k*M1-M2)...
    -A/4*M3;

%% Calculate eigenstates and overlaps

Ht=-M3/4;
% Ht=M3;
Cmat=zeros(numStates,numStates,length(kvec));
for kk=1:kpoints
    [v,b]=eig(full(Hmat(kvec(kk))));   
%     eng=b*ones(size(v,1),1);                    % Enegies    
%     lambdaFull=eng-depth/2;                     % Offset    

    Cmat(:,:,kk)=abs(ctranspose(v)*Ht*v);
end
%% Plot it

hF=figure(1736);
set(hF,'color','w','units','pixels','resize','off');
clf
hF.Position=[380 75 400 900];

lst=combnk(1:maxBand,2);
lst=flip(lst);


ps=cell(1,maxBand-1);

pds=cell(1,maxBand-1);
for nn=1:size(lst,1)  
    n1=lst(nn,1);
    n2=lst(nn,2);    
    
    ps{n1}=[ps{n1} n2];
    
    subplot(maxBand-1,1,maxBand-n1);      
    co=get(gca,'colororder');
    
    yy=reshape(Cmat(n1,n2,:),[kpoints 1]);
    pds{n1}{end+1}=plot(kvec,yy,'linewidth',2,'color',co(n2,:));    
    hold on
end

for nn=1:(maxBand-1)
    subplot(maxBand-1,1,maxBand-nn);
    ylim([0 .35]);
    set(gca,'FontSize',14,'box','on',...
        'linewidth',1);
    str=num2str(nn);
    ylabel(str);
    set(gca,'YTickLabel',[]);
    
    set(get(gca,'YLabel'),'color',co(nn,:));
    
    if nn==(maxBand-1)
       title('\langle i | H_{AM} | j \rangle'); 
    end

   if nn==1
       xlabel('Quasimomentum');
   else
       set(gca,'XTickLabel',[]);
   end
    hold on   
   p=get(gca,'Position');
   l1=legend(strsplit(num2str(ps{nn})),...
       'orientation','horizontal','location','northwest');
   l1.Position(1)=p(1);
   l1.Position(2)=p(2)+p(4)-l1.Position(4);
   l1.AutoUpdate='off';
end

hold all

%% Function and Special plots

pos={};

funcs=struct;
funcs.Add=@addTransition;
funcs.Clr=@clearTransitions;

    function clearTransitions
        for n=1:length(pos)
           delete(pos{n}); 
        end
    end

    function addTransition(nI,nF,k)
        X=pds{nI}{nF-nI}.XData;
        ind=find(k<X,1);
        
        x=pds{nI}{nF-nI}.XData(ind);
        y=pds{nI}{nF-nI}.YData(ind);
        
        figure(hF)
        subplot(maxBand-1,1,maxBand-nI);
        pos{end+1}=plot([-1 1]*x,[1 1]*y,...
            'kx','markersize',8,'linewidth',2);
        
    end


end



