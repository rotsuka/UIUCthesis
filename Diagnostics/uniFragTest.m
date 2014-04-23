% uniFragTest: create a deterministic fragility model with multiple
% fragility curves and simulate probabilitic draws from the distribution
%
% NOTE: NOT A FUNCTION

tic

clear all;close all;clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inputs

M=8 % magnitude
R=[10] % distance (km)
numBr=length(R); % number of bridges

% Table 6 from paper. MSC steel bridge. 1=insignificant, 2=moderate,
% 3=complete (from Nielson)

% MSC Steel
med=[.18 .31 .5];
dsp=[.55 .55 .55];

% % MSC Slab
% l=[.17 .45 1.73];
% b=[.7 .7 .7];

% % MSC Conc
% l=[.15 .52 1.03];
% b=[.70 .70 .70];

% % MSSS Conc
% l=[.20 .57 1.17];
% b=[.65 .65 .65];

% % MSSS Steel
% l=[.24 .44 .82];
% b=[.5 .5 .5];

delta=.001;
SF=1/delta; % scale factor for rounding

% Vector of Sa inputs to build curve
PGAVec=0:delta:2;

% 100% traffic capacity
qmax=2000;

% Number of draws (ensembles)
noEn=100;

F=0;
S_SR=[0 0 0 0];
S_HR=[0 0 0 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% METHODS

% Uses Campbell 1997 attenuation law to get the attenuation. Sa can be a
% vector
PGA=getPGA(M,R,F,S_SR,S_HR)

% Fragility curves 
p=getFragCurves(PGAVec,med,dsp);

% Takes the Sa from the attenuation law and gives the corresponding
% probability of exceedance
Pex=getPex(PGA,p,SF,PGAVec);

% Take the probability of exceedance and compute the probability of damage
% states as a vector
ProbDS=getProbVec(Pex)

% Get the range of the four limit states
DSlims=getDSRange(ProbDS)

% Simulate damage states
for n=1:noEn
    
    for i=1:numBr
        
        % Draw from the DS distribution
        draw=rand;
        
        % Determine the damage state
        [ratio,DState]=compCap(DSlims(i,:),draw);
        
        % Compute new traffic capacity
        cap=ratio*qmax;
        
        disp(['Bridge ' num2str(i) ', Draw ' num2str(n) ': ' num2str(draw) ', DS: ' DState ', Capacity: ' num2str(cap)]);
        
    end
    
    disp('...');
    
end

% Plot fragility curves
plot(PGAVec,p(1,:),PGAVec,p(2,:),'--',PGAVec,p(3,:),':')
% set(gcf,'defaulttextinterpreter','latex');
title('Fragility Curves');
xlabel('A_H (g)');
ylabel('Probability of exceedence');
legend('Slight','Moderate','Complete');
set(gca,'fontsize',9)
set(findall(gcf,'type','text'),'fontsize',9);



toc