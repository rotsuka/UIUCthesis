% FundDiaPlot2: plots fundamental diagram and sending and receiving diagrams
% with noise
%
% NOTE: Not a function. Uses an rng variable to get similar noise
% everytime.

close all; clear all; clc; 
load myrng.mat
rng(notrand)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs

numLanes=2; % 1.5, 1, 0
qmax=2400;
rhoj=144;
vmax=100;

% Vector
a=0;
b=rhoj*numLanes;
rho=a:.1:b;

% Sending and receiving function uncertainty
mu_s=150;
mu_q=200;
mu_r=300;
var_s=100^2;
var_q=150^2;
var_r=250^2;

% Lognormal params
err_S=log(mu_s^2/sqrt(mu_s^2+var_s)) % lambda
err_Q=log(mu_q^2/sqrt(mu_q^2+var_q)) % lambda
err_R=log(mu_r^2/sqrt(mu_r^2+var_r)) % lambda
Q_S=log(1+var_s/mu_s^2) % xi^2
Q_Q=log(1+var_q/mu_q^2) % xi^2
Q_R=log(1+var_r/mu_r^2) % xi^2

% Amount of noise entries  (total)
num=5000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Upper bounds

qm=qmax*numLanes;
rj=rhoj*numLanes;
rc=qm/vmax;
    
for i=1:length(rho)
    
    % Separate
    R=receiving_q(rho(i),rhoj,qmax,vmax,numLanes,err_S,err_Q,Q_S,Q_Q,0);
    S=sending_q(rho(i),rhoj,qmax,vmax,numLanes,err_R,err_Q,Q_S,Q_Q,0);
    rec_array(i)=R;
    send_array(i)=S;
    
    % Together
    if rho(i)<=rc
        q=vmax*rho(i);
    else
        if rj==0
            q=0;
        else
            q=(qm/(rj-rc))*(rj-rho(i));
        end
    end
    fund_array(i)=q;
    
end

% Noisy

randrhos=a+(b-a).*rand(num,1);
randrhor=a+(b-a).*rand(num,1);
s_noise=zeros(size(randrhos));
r_noise=zeros(size(randrhor));

count=1;

for i=1:num
    
    s_noise(i)=sending_q_neg(randrhos(i),rhoj,qmax,vmax,numLanes,err_S,err_Q,Q_S,Q_Q,1);
    r_noise(i)=receiving_q_neg(randrhor(i),rhoj,qmax,vmax,numLanes,err_R,err_Q,Q_R,Q_Q,1);
        
    if randrhos(i)<rc
        fund_noise(count)=s_noise(i);
        randrhof(count)=randrhos(i);
        count=count+1;
    end
    
    if randrhor(i)>=rc
        fund_noise(count)=r_noise(i);
        randrhof(count)=randrhor(i);
        count=count+1;
    end
     
end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot
   
figure 

% subplot(1,3,1)
plot(randrhos,s_noise,'r.','MarkerSize',4)
hold on
plot(rho,send_array)
hold off
set(gcf,'defaulttextinterpreter','latex');
xlabel('$\rho$')
ylabel('$S(\rho)$')
legend('Noisy','Deterministic')
ylim([0 5000])

figure 

% subplot(1,3,2)
plot(randrhor,r_noise,'r.','MarkerSize',4)
hold on
plot(rho,rec_array)
hold off
set(gcf,'defaulttextinterpreter','latex');
xlabel('$\rho$')
ylabel('$R(\rho)$')
legend('Noisy','Deterministic')
ylim([0 5000])

figure 

% subplot(1,3,3)
plot(randrhof,fund_noise,'r.','MarkerSize',4)
hold on
plot(rho,fund_array)
hold off
set(gcf,'defaulttextinterpreter','latex');
xlabel('$\rho$')
ylabel('$Q(\rho)$')
legend('Noisy','Deterministic')
ylim([0 5000])