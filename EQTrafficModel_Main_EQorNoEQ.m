% EQTrafficModel_Main_EQorNoEQ: the main file to execute. This will run a
% scenario with an earthquake or without an earthquake.

tic

close all; clear all; clc

% Warnings:
warning('Triangular fundamental diagram is assumed');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % INPUT PARAMETERS (USER INPUTS)
% 
% % Nodes
% 
% % NOTE: mapNodes does not need to be copied so passing the map back through
% % the function should be fine, but pass back valueNodes just in case
% numNodes=5;
% [valueNodes,keyNodes]=inpNodes(numNodes);
% 
% % Links
% 
% % Pass back the cell of links instead of the actual map so multiple map
% % object with difrerent handles can be created
% numLinks=6; % internal/boundary combined total
% [valueLinks,keyLinks]=inpLinks(numLinks);
% 
% % Sensors
% 
% %NOTE: mapSensors does not need to be copied so passing the map back through
% % the function should be fine, but pass back valueSensors just in case
% numSensors=2;
% [valueSensors,keySensors]=inpSensors(numSensors);
% 
% % EQ 
% newEQ=inpEQ;
% 
% % Tell which links are actually bridges
% bridgeLinks=inpBridge;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PRELOADED ROAD GEOMETRY

% load Geom_4link5node_2Sens_SensNoise5_MD.mat
load Geom_4link5node_2Sens_SensNoise5_HI.mat
% load Geom_4link5node_2Sens_SensNoise5_TO.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1 if there is EQ input in the evolution model
isEQInp=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SIMULATION AND FILTER VARIABLES

tmax=2400; % time of simulation [seconds]
DT=5; % time step [seconds]

% Obtain time vector
measTime=0:DT:tmax;

% Time at which EQ occurs (make sure it is a factor of tmax and a multiple of DT)
teq=600;
keq=(teq/DT)+1; % Kalman time step at which EQ occurs

% EnKF Variables
noEn=100; % number of ensemble members
tolerance=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BRIDGE VARIABLES (Nielson 2007) 

% Vector of PGA inputs to build fragility curve. These variables describe
% how accurately we want to map computed PGA values to our fragility curves
delta=.001; % The degree of accuracy
SF=1/delta; % Scale factor for rounding purposes
PGAVec=0:delta:3; % Choose a reasonable range

% % MSC Concrete
% med=[.15 .52 1.03];
% dsp=[.70 .70 .70];

% MSC Steel
med=[.18 .31 .5];
dsp=[.55 .55 .55];

% % MSSS Concrete
% med=[.20 .57 1.17];
% dsp=[.65 .65 .65];
 
% % MSSS Steel
% med=[.24 .44 .82];
% dsp=[.5 .5 .5];

% Parameters for attenuation (Refer to Campbell 1997)
F=0; % Numerical value for type of fault, scalar
S_SR=[0]; % Soil param, can be a vector if multiple bridges
S_HR=[0]; % Soil param, can be a vector if multiple bridges

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIAL AND BOUNDARY CONDITIONS, AND NOISE MODELS 

% Inital conditions. Currently the order is [internal links;left BC;right
% BC]

% Of accurate/noise solution (TOTAL NOT PER LANE)
rho0Acc=[40;35;35;35;60;20];

% Of approx (Kalman) solution (TOTAL NOT PER LANE)
% MD
% rho0App=[20;20;30;30;40;20]; % Test 1

% HI
rho0App=[20;20;30;30;40;20]; % Test 1

% TO
% rho0App=[20;20;30;30;40;20]; % Test 1

% Boundary conditions for acc and app models. This is dependent on the
% order of rho0App and rho0Acc
startAcc=rho0Acc(end-1); % left BC
endAcc=rho0Acc(end); % right BC
startApp=rho0App(end-1); % left BC
endApp=rho0App(end); % right BC

% Boundary condition noise (density percentages)
startSig=25;
endSig=25;
noiseBC=[startSig endSig];

% Sending and receiving function uncertainty (units of flow) for noise
% model
mu_s=150; % sending function
mu_q=200; % qmax value
mu_r=300; % receiving function
var_s=100^2; % variance on sending
var_q=150^2; % varience on qmax 
var_r=200^2; % variance on receiving

% Lognormal params. Do not alter this code as it provides the conversion
% between the normal distribution parameters and lognormal distribution.
err_S=log(mu_s^2/sqrt(mu_s^2+var_s)); % lambda
err_Q=log(mu_q^2/sqrt(mu_q^2+var_q)); % lambda
err_R=log(mu_r^2/sqrt(mu_r^2+var_r)); % lambda
Q_S=log(1+var_s/mu_s^2); % xi^2
Q_Q=log(1+var_q/mu_q^2); % xi^2
Q_R=log(1+var_r/mu_r^2); % xi^2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CREATE MULTIPLE MAPLINKS/MAPNODES FOR EXACT/ACCURATE ETC

% For true solution (DT,rho0Acc)
mapLinksNoise=containers.Map(keyLinks,valueLinks);

% For Kalman solution (DT,rho0App). A second map object is necessary so that
% the open loop model can be run.

% % For these solutions slightly perturb values of true parameters.
% nvalueLinks=alterParams(valueLinks);
% mapLinksApp=containers.Map(keyLinks,nvalueLinks); % filter, EQ
% mapLinksApp2=containers.Map(keyLinks,nvalueLinks); % open loop, EQ

mapLinksApp=containers.Map(keyLinks,valueLinks);
mapLinksApp2=containers.Map(keyLinks,valueLinks);

% mapNodes
mapNodes=containers.Map(keyNodes,valueNodes);

% mapSensors (ONLY TRAFFIC DENSITY SENSORS)
mapSensors=containers.Map(keySensors,valueSensors);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% COMPUTE NECESSARY QUANTITIES FOR EACH MAPLINKS CONTAINER

% Compute length of each link and overwrite previous map obj.
mapLinksNoise=getLengthCoord(mapLinksNoise,mapNodes,numLinks);
mapLinksApp=getLengthCoord(mapLinksApp,mapNodes,numLinks);
mapLinksApp2=getLengthCoord(mapLinksApp2,mapNodes,numLinks);

% Compute dx,number of cells, and start/end cell for link using CFL
% condition and overwrite previous map obj.
mapLinksNoise=getCellDisc(DT,mapLinksNoise,keyLinks,numLinks);
mapLinksApp=getCellDisc(DT,mapLinksApp,keyLinks,numLinks);
mapLinksApp2=getCellDisc(DT,mapLinksApp2,keyLinks,numLinks);

% The total number of cells takes into account the CTM model and the two
% ghost cells.
totCells=compTotCells(mapLinksNoise,numLinks);

% Add the initial conditions and overwrite previous map obj
mapLinksNoise=addInitialCond(mapLinksNoise,keyLinks,rho0Acc);
mapLinksApp=addInitialCond(mapLinksApp,keyLinks,rho0App);
mapLinksApp2=addInitialCond(mapLinksApp2,keyLinks,rho0App);

% Determine which links are bridges and overwrite previous map obj
mapLinksNoise=getBridges(mapLinksNoise,numLinks,bridgeLinks);
mapLinksApp=getBridges(mapLinksApp,numLinks,bridgeLinks);
mapLinksApp2=getBridges(mapLinksApp2,numLinks,bridgeLinks);

% Attach links onto nodes and overwrite previous map obj
mapNodes=attachLinks(mapNodes,numNodes,mapLinksApp,numLinks);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MEASUREMENT MATRIX, H

% Rearrange sensor order and give a global cell position. keySensors here
% gives the keys back in the desired order.
[mapSensors keySensors sortInd]=getSensorPos(mapSensors,mapLinksApp,...
    numSensors);

% Attach sensors onto respective links. No need to attach sensors to
% mapLinksApp2.
mapLinksApp=attachSensors(mapLinksApp,mapSensors,numSensors,keySensors);

% Get measurement matrix, H and its dimensions, kk and n
[H kk n]=getHSens(mapSensors,keySensors,totCells,numSensors);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RUN CTM (forward progression) for acc/noise solutions.

% Initialize. NOTE: Parameters of the true model are given a "Noise" or "N"
% ending in their name.
xMatNoise=zeros(length(measTime),totCells);

% Initial conditions. Implementation depends on the mapSensors being crated
% in order of cell discretization
x0Noise=getStates0(mapLinksNoise,totCells,numLinks);

% For true model
xMatNoise(1,:)=x0Noise';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For true model, determine qmax due to damage from the earthquake. 'N'
% here indicates those variables related with the true model.

% Compute bridge location w.r.t. epicenter
mapLinksNoise=getDistanceEQ(mapLinksNoise,numLinks,newEQ);

% Construct vector of distances
distVecN=getDistVec(mapLinksNoise,keyLinks,bridgeLinks,numLinks);

% Get peak ground accelerations
PGAN=getPGA(newEQ.mag,distVecN,F,S_SR,S_HR);

% Create deterministic fragility curves for fragility model. 
p=getFragCurves(PGAVec,med,dsp);

% Takes the Sa from the attenuation law and gives the corresponding
% probability of exceedance
PexN=getPex(PGAN,p,SF,PGAVec);

% Take the probability of exceedance and compute the probability of damage
% states as a vector
ProbDSN=getProbVec(PexN);

% Get the range of the four limit states
DSlimsN=getDSRange(ProbDSN);

% Forward CTM model
for k=1:length(measTime)-1
    
    % Assign the parameter based on the probability table. We only want to
    % assign ONCE at the instance the EQ occurs
    if k==keq
        
        mapLinksNoise=getMapParamsN(mapLinksNoise,numLinks,keyLinks,DSlimsN);
        
    end
    
    % Density propogation via CTM
    xMatNoise(k+1,1:totCells)=forwardPredict(DT,xMatNoise(k,:),mapLinksNoise,...
        mapNodes,totCells,startAcc,endAcc,0,err_S,err_Q,err_R,Q_S,Q_Q,Q_R,noiseBC);
    
end

% Uncomment if you want to just run deterministic model
% densPlotter_one(xMatNoise,mapLinksNoise,numLinks,measTime);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% KALMAN FILTER PREPARATION

% Create a measurement noise matrix.
noisyMat=getNoisyMat(mapSensors,numSensors,keySensors,measTime);

% Get the measurement vector by adding the noisyMat to the true solution 
y=getyMeas(H,xMatNoise,measTime,numSensors,noisyMat);

% Kalman time
K=length(measTime);

% Initialize 3D arrays and 4D arrays that contain density values and
% covariance errors for all ensembles and all time steps (m=prior,
% p=posterior)
[xk_m xk_p Pk_m Pk_p]=initialize3DArray(n,K);
[xi_m xi_p]=initialize4DArray(n,noEn,K);

% qmaxArray to keep track of qmax of every ensemble for every link
qmaxArray=zeros(numLinks,1,noEn,K);

% Averaged measurement noise
v_bar=zeros(numSensors,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GET INITIAL STATES AND INITIAL ERROR COV OF STATES

% Initial states
xk_p(1:totCells,1,1)=getStates0(mapLinksApp,totCells,numLinks);

x0App=xk_p(:,1,1);

% Initial error covariance matrix
Pk_p(:,:,1)=initialCov(x0Noise,x0App,tolerance,totCells);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% THE SUBSCRIPT 2 IS FOR THE FORWARD PROPOGATION OF ENSEMBLES WITH NO
% SENSOR DATA (OPEN LOOP)

% Initialize (no need to create add'l model noise or meas noise arrays)
[xk_m2 xk_p2 Pk_m2 Pk_p2]=initialize3DArray(n,K);

xk_p2(:,1,1)=xk_p(:,1,1);

[xi_m2 xi_p2]=initialize4DArray(n,noEn,K);

qmaxArray2=zeros(numLinks,1,noEn,K);

% Duplicate EQ object
newEQ2=newEQ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If you dont do this overwriting, qmax and rhoj will both be 0, which will
% cause distDraw0 to go in an infinite loop.
mapLinksApp=setInitialParams(mapLinksApp,numLinks);
mapLinksApp2=setInitialParams(mapLinksApp2,numLinks);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RUN FILTER

redrawCount=0;

% Initial ensembles draws

for i=1:noEn
    
    % distDraw0 gives only the density draws
    [xi_p(1:totCells,1,i,1) redrawCount]=distDraw0(xk_p(:,1,1),Pk_p(:,:,1),mapLinksApp,...
        numLinks,totCells,redrawCount,i);
    
    % Open loop solution has same inital ensembles
    xi_p2(:,1,i,1)=xi_p(:,1,i,1);
    
end

redrawCount

% Initialize X5 and filtered ensembles arrays. This allows every ensemble
% to be reproduced at every step, even after the simulation has been
% performed.
X5_array=cell(K,1);
Af_array=cell(K,1);
D_array=cell(K,1);
Gam_array=cell(K,1);

% Innovation and Kalman gain arrays
innov_array=cell(K,1);
kalGain_array=cell(K,1);

% Initial filtered ensembles in array form
Af_i=composeEnsToMat(xi_p,totCells,noEn,1);
Af_array{1}=Af_i;

for k=1:K-1

%     % If you want to make sure code is processing
%     if mod(k,10)==0
%         k
%         toc
%     end
    
    % MODEL PREDICTION (PRIOR)
    
    for i=1:noEn
        
        % FRAGILITY FUNCTIONS
        
        % If the EQ can occur
        if k>=keq
            
            % Get the EQ characteristics. Note that distribution parameters
            % are defined in the function itself
            [m xloc yloc]=getEQParams(isEQInp,newEQ);
            [m2 xloc2 yloc2]=getEQParams(isEQInp,newEQ);
            
            % Create the earthquakes
            nEQ=editEQ(m,xloc,yloc);
            nEQ2=editEQ(m2,xloc2,yloc2);
            
            % Compute bridge location w.r.t. epicenter
            mapLinksApp=getDistanceEQ(mapLinksApp,numLinks,nEQ);
            mapLinksApp2=getDistanceEQ(mapLinksApp2,numLinks,nEQ2);
            
            % Construct vector of distances, R
            distVec=getDistVec(mapLinksApp,keyLinks,bridgeLinks,numLinks);
            distVec2=getDistVec(mapLinksApp2,keyLinks,bridgeLinks,numLinks);
 
            % Get spectral accelerations
            PGA=getPGA(nEQ.mag,distVec,F,S_SR,S_HR);
            PGA2=getPGA(nEQ2.mag,distVec2,F,S_SR,S_HR);
            
            % Takes the Sa from the attenuation law and gives the corresponding
            % probability of exceedance
            Pex=getPex(PGA,p,SF,PGAVec);
            Pex2=getPex(PGA2,p,SF,PGAVec);
            
            % Take the probability of exceedance and compute the probability of damage
            % states as a vector
            ProbDS=getProbVec(Pex);
            ProbDS2=getProbVec(Pex2);
            
            % Get the range of the four limit states
            DSlims=getDSRange(ProbDS);
            DSlims2=getDSRange(ProbDS2);
            
            % Assign the parameter based on the probability table (for each
            % ensemble)
            mapLinksApp=getMapParams(mapLinksApp,numLinks,keyLinks,DSlims);
            mapLinksApp2=getMapParams(mapLinksApp2,numLinks,keyLinks,DSlims2);
            
        end
        
        % Code to extract qmax values from each link (even those that arent
        % bridges) for each ensemble
        qmaxArray(:,1,i,k+1)=getqmaxMat(mapLinksApp,keyLinks);
        qmaxArray2(:,1,i,k+1)=getqmaxMat(mapLinksApp2,keyLinks);
        
        % DENSITY AND PARAM MODELS
        
        % Forward progression of densities via CTM
        xi_m(1:totCells,1,i,k+1)=forwardPredictEns(DT,xi_p(:,1,i,k),mapLinksApp,...
            mapNodes,totCells,startApp,endApp,1,err_S,err_Q,err_R,Q_S,Q_Q,Q_R,noiseBC);
        xi_m2(1:totCells,1,i,k+1)=forwardPredictEns(DT,xi_p2(:,1,i,k),mapLinksApp2,...
            mapNodes,totCells,startApp,endApp,1,err_S,err_Q,err_R,Q_S,Q_Q,Q_R,noiseBC);
        
    end
    
    % Compute A, and unbiased covariance matrix (Evensen)
    
    A=composeEnsToMat(xi_m,totCells,noEn,k+1);
    
    % Compute ensemble means
    Abar=A*1/noEn*ones(noEn); 
    
    % Put means into desired notation
    xk_m(:,1,k+1)=Abar(:,1);
   
    Apr=A-Abar;
    
    Pk_m(:,:,k+1)=(1/(noEn-1))*(Apr*Apr');
       
    A2=composeEnsToMat(xi_m2,totCells,noEn,k+1);
    
    Abar2=A2*1/noEn*ones(noEn);
    xk_m2(:,1,k+1)=Abar2(:,1);
    Apr2=A2-Abar2;
    
    Pk_m2(:,:,k+1)=(1/(noEn-1))*(Apr2*Apr2');

    % ANALYSIS SCHEME (Evensen)
    
    D=zeros(numSensors,noEn);
    Gam=zeros(numSensors,noEn);
    
    for i=1:noEn
        
        meas=y(k+1,:)';
       
        % Get ensemble noise and noisy measurement
        [ensMeas ensMeasNoise]=getEnsMeasNoise(meas,mapSensors,keySensors);
        
        % Add to D and Gam matrices
        D(:,i)=ensMeas;
        Gam(:,i)=ensMeasNoise;
        
    end
    
    % Compute Re
    Re=(1/(noEn-1))*(Gam*Gam');
    
    % Compute D' 
    Dpr=D-H*A;
    
    % Compute innovation and Kalman gain
    [innov_array{k+1},kalGain_array{k+1},D_array{k+1},Gam_array{k+1}]=...
        getUpdateInfo(H,Pk_m(:,:,k+1),Re,Dpr,D,Gam);
   
    % Singular value decomposition (SVD)
    
    [U,Sig,V]=svd(H*Apr+Gam);
    Lam=Sig*Sig';
    
    % Compute update matrix, Aa, from X scheme 
    
    X1=Lam\U';
    X2=X1*Dpr;
    X3=U*X2;
    X4=(H*Apr)'*X3;
    X5=eye(noEn)+X4;
    
    Aa=A*X5;
    
    % Add X5 and Aa to their respective arrays (for smoother)
    X5_array{k+1}=X5;
    Af_array{k+1}=Aa;
    
    % Do not do analysis scheme!
    Aa2=A2;
    
    % Decompose Aa into ensemble vectors
    
    for i=1:noEn
        
        xi_p(:,1,i,k+1)=decompMatToEns(Aa,i);
        xi_p2(:,1,i,k+1)=decompMatToEns(Aa2,i);
   
    end

    % Compute the unbiased covariance matrix for updated ensembles
    
    Aabar=Aa*1/noEn*ones(noEn);
    
    % Put means into desired notation
    xk_p(:,1,k+1)=Aabar(:,1);
    Aapr=Aa-Aabar;
    Pk_p(:,:,k+1)=(1/(noEn-1))*(Aapr*Aapr');
   
    Aabar2=Aa2*1/noEn*ones(noEn);
    
    % Put means into desired notation
    xk_p2(:,1,k+1)=Abar2(:,1);
    Aapr2=Aa2-Aabar2;
    Pk_p2(:,:,k+1)=(1/(noEn-1))*(Aapr2*Aapr2');
       
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert to matrix form of density means

xMatm=conv2Mat(xk_m);
xMatp=conv2Mat(xk_p);
xMatp2=conv2Mat(xk_p2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc