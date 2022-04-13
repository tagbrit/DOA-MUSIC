% Main code for generating signals and plotting DOA estimation
% Includes the MUSIC function at the end

clear;
clc;
close all; 

%% Signal Generation and Setup

N = 8; % number of sensors
K = 1; % number of sources
T = 1000; % number of snapshots
R = 20; % SNR
theta = [-50,0,60]*pi/180; % angle of arrival
d = 0.5; % distance between sources

A = zeros(N,K); % generation of the steering matrix 
for n=1:N
    for k=1:K
     A(n,k)=exp(-1i*2*pi*(n-1)*(d)*sin(theta(k)));
    end
end

S = randn(K,T); % randomly generated source signals
AS = A*S; % multiplying source signals and steering matrix
X = awgn(AS,R); % adds white gaussian noise to AS with certain SNR

[music] = MUSIC(X,K);  % call for MUSIC function


%% Plotting

theta_scan_deg = -90:0.05:90; 
theta_scan = theta_scan_deg*pi/180;
figure(1); 
plot(theta_scan*180/pi,10*log10(music)); hold on
grid on;
title('MUSIC Spectrum');
xlabel('degree');
ylabel('Magnitude (dB)');
xlim([-70,70])

%% MUSIC Function

function [music] = MUSIC(X,K)
    [~,T]=size(X); % want the number columns
    Rx = (1/T)*X*X'; % covariance matrix
    [V,~] = eig(Rx); % eigen decomp of Rx
    V2 = V(:,K+1:end); % estimate noise subspace
    [N,~] = size(V2); % want the number of rows
    theta_scan_deg = -90:0.05:90;
    theta_scan = theta_scan_deg*pi/180;
    for i=1:length(theta_scan)
        a= exp(-1i*pi*sin(theta_scan(i))*(0:1:(N-1))'); % create steering column vector
        music(i)= 1/(norm(V2'*a,2).^2); % MUSIC spatial spectrum
    end
end
