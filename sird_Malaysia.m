clear model data params options
datafile = readmatrix('sird_data_Malaysia.csv');
Ptotal = 3.153*10^7;%malaysia
rng = 3683;
Idata = datafile(13:end,3);
Rdata = datafile(13:end,2);
Ddata = datafile(13:end,4);
Sdata = datafile(13:end,5);
tdata = datafile(13:end,1);
t = 1:1:length(Idata);
tdays = datetime(2020,3,26) + caldays(1:length(t));
Pastparamestimate = readmatrix('Malaysia-SIRD-param-estimate.csv');
Szero = Sdata(1);
Izero = Idata(1);
Rzero = Rdata(1);
Dzero = Ddata(1);
T = length(t);
N = T;
M = 500;

% sigma = 0;
dt = T/N; R = 1; Dt = R*dt; L = N/R;
% mu = 0.000021;
mu = 0;
%======================================================

%======================================================
% 4-stage SRK Tableus form %Refer Equation 4.55
A = [0 0 0 0;1/2 0 0 0;0 1/2 0 0;0 0 1 0]; %A
B_1 = [0 0 0 0;-0.7242916356 0 0 0;0.4237353406 -0.1994437050 0 0;-1.578475506 0.840100343 1.738375163 0]; %B^(1)
B_2 = [0 0 0 0;2.702000410 0 0 0;1.757261649 0 0 0;-2.918524118 0 0 0]; %B^(2)
G_1 = [-.7800788474;0.07363768240;1.486520013;0.2199211524]; %GAMMA 1
G_2 = [1.693950844;1.636107882;-3.024009558;-0.3060491602]; %GAMMA 2
% alpha =[1/6; 1/3; 1/3; 1/6]; %ALPHA
%======================================================
% plot(Datay(1,:),Datay(2,:),'ro');hold on %plot data
% plot(Datax(1,:),Datax(2,:),'bo');hold on %plot data
%======================================================

% alpha =[1/6; 1/3; 1/3; 1/6]; %ALPHA
%======================================================
% plot(Datay(1,:),Datay(2,:),'ro');hold on %plot data
% plot(Datax(1,:),Datax(2,:),'bo');hold on %plot data
fS = @(S,I,R,D,Beta,r,d,gamma,sigma,mu) mu*(S+I+R)-S.*I./Ptotal.*Beta+gamma.*R - mu.*S;
gS = @(S,I,R,D,Beta,r,d,gamma,sigma,mu) -sigma.*I.*S./Ptotal;
fI = @(S,I,R,D,Beta,r,d,gamma,sigma,mu) S./Ptotal.*Beta.*I-(r+d).*I - mu.*I;
gI = @(S,I,R,D,Beta,r,d,gamma,sigma,mu) sigma.*I.*S./Ptotal;
fR = @(S,I,R,D,Beta,r,d,gamma,sigma,mu) r.*I - gamma.*R - mu.*R;
gR = @(S,I,R,D,Beta,r,d,gamma,sigma,mu) 0;
fD = @(S,I,R,D,Beta,r,d,gamma,sigma,mu) d.*I;
gD = @(S,I,R,D,Beta,r,d,gamma,sigma,mu) 0;
%% SRK4==================================================
sigma = 0.05;
%======================================================%SRK-4 Method
Ssrk4_2 = zeros(M,N);
Isrk4_2 = zeros(M,N);
Rsrk4_2 = zeros(M,N);
Dsrk4_2 = zeros(M,N);
Stemp = Szero*ones(M,1);
Itemp = Izero*ones(M,1);
Rtemp = Rzero*ones(M,1);
Dtemp = Dzero*ones(M,1);
Smodel = ones(3,N);
Imodel = ones(3,N);
Rmodel = ones(3,N);
Dmodel = ones(3,N);
SdW1 = sqrt(dt)*randn(M,N);       %Numerical Brownian increments
SdW2 = sqrt(dt)*randn(M,N);       %Numerical Brownian increments
IdW1 = sqrt(dt)*randn(M,N);       %Numerical Brownian increments
IdW2 = sqrt(dt)*randn(M,N);       %Numerical Brownian increments
RdW1 = sqrt(dt)*randn(M,N);       %Numerical Brownian increments
RdW2 = sqrt(dt)*randn(M,N);       %Numerical Brownian increments
DdW1 = sqrt(dt)*randn(M,N);       %Numerical Brownian increments
DdW2 = sqrt(dt)*randn(M,N);       %Numerical Brownian increments

for j = 1:L
    beta = Pastparamestimate(j,2);
    r = Pastparamestimate(j,3);
    d = Pastparamestimate(j,4);
    gamma = Pastparamestimate(j,5);

    SWinc1 = sum(SdW1(:,R*(j-1)+1:R*j),2);
    SWinc2 = sum(SdW2(:,R*(j-1)+1:R*j),2);
    IWinc1 = sum(IdW1(:,R*(j-1)+1:R*j),2);
    IWinc2 = sum(IdW2(:,R*(j-1)+1:R*j),2);
    RWinc1 = sum(RdW1(:,R*(j-1)+1:R*j),2);
    RWinc2 = sum(RdW2(:,R*(j-1)+1:R*j),2);
    DWinc1 = sum(DdW1(:,R*(j-1)+1:R*j),2);
    DWinc2 = sum(DdW2(:,R*(j-1)+1:R*j),2);

    SJ10 = (0.5*(Dt^(1/2))*(SWinc1+(1/sqrt(3)).*SWinc2));
    IJ10 = (0.5*(Dt^(1/2))*(IWinc1+(1/sqrt(3)).*IWinc2));
    RJ10 = (0.5*(Dt^(1/2))*(RWinc1+(1/sqrt(3)).*RWinc2));
    DJ10 = (0.5*(Dt^(1/2))*(DWinc1+(1/sqrt(3)).*DWinc2));

    S1 = Stemp;
    I1 = Itemp;
    R1 = Rtemp;
    D1 = Dtemp;
    %Refer Equation 4.56
    f1_S = fS(S1, I1, R1, D1, beta, r, d, gamma, sigma,mu);
    g1_S = gS(S1, I1, R1, D1, beta, r, d, gamma, sigma,mu);
    f1_I = fI(S1, I1, R1, D1, beta, r, d, gamma, sigma,mu);
    g1_I = gI(S1, I1, R1, D1, beta, r, d, gamma, sigma,mu);
    f1_R = fR(S1, I1, R1, D1, beta, r, d, gamma, sigma,mu);
    g1_R = gR(S1, I1, R1, D1, beta, r, d, gamma, sigma,mu);
    f1_D = fD(S1, I1, R1, D1, beta, r, d, gamma, sigma,mu);
    g1_D = gD(S1, I1, R1, D1, beta, r, d, gamma, sigma,mu);

    S2 = Stemp + Dt*(A(2)*f1_S) + (B_1(2).*SWinc1 + B_2(2).*SJ10).*g1_S;
    I2 = Itemp + Dt*(A(2)*f1_I) + (B_1(2).*IWinc1 + B_2(2).*IJ10).*g1_I;
    R2 = Rtemp + Dt*(A(2)*f1_R) + (B_1(2).*RWinc1 + B_2(2).*RJ10).*g1_R;
    D2 = Dtemp + Dt*(A(2)*f1_D) + (B_1(2).*DWinc1 + B_2(2).*DJ10).*g1_D;

    %Refer Equation 4.56
    f2_S = fS(S2, I2, R2, D2, beta, r, d, gamma, sigma,mu);
    g2_S = gS(S2, I2, R2, D2, beta, r, d, gamma, sigma,mu);
    f2_I = fI(S2, I2, R2, D2, beta, r, d, gamma, sigma,mu);
    g2_I = gI(S2, I2, R2, D2, beta, r, d, gamma, sigma,mu);
    f2_R = fR(S2, I2, R2, D2, beta, r, d, gamma, sigma,mu);
    g2_R = gR(S2, I2, R2, D2, beta, r, d, gamma, sigma,mu);
    f2_D = fD(S2, I2, R2, D2, beta, r, d, gamma, sigma,mu);
    g2_D = gD(S2, I2, R2, D2, beta, r, d, gamma, sigma,mu);

    S3 = Stemp + Dt*(A(7)*f2_S) + (B_1(3).*SWinc1 + B_2(3).*SJ10).*g2_S + B_1(7)*SJ10.*g2_S;
    I3 = Itemp + Dt*(A(7)*f2_I) + (B_1(3).*IWinc1 + B_2(3).*IJ10).*g2_I + B_1(7)*IJ10.*g2_I;
    R3 = Rtemp + Dt*(A(7)*f2_R) + (B_1(3).*RWinc1 + B_2(3).*RJ10).*g2_R + B_1(7)*RJ10.*g2_R;
    D3 = Dtemp + Dt*(A(7)*f2_D) + (B_1(3).*DWinc1 + B_2(3).*DJ10).*g2_D + B_1(7)*DJ10.*g2_D;

    %Refer Equation 4.56
    f3_S = fS(S3, I3, R3, D3, beta, r, d, gamma, sigma,mu);
    g3_S = gS(S3, I3, R3, D3, beta, r, d, gamma, sigma,mu);
    f3_I = fI(S3, I3, R3, D3, beta, r, d, gamma, sigma,mu);
    g3_I = gI(S3, I3, R3, D3, beta, r, d, gamma, sigma,mu);
    f3_R = fR(S3, I3, R3, D3, beta, r, d, gamma, sigma,mu);
    g3_R = gR(S3, I3, R3, D3, beta, r, d, gamma, sigma,mu);
    f3_D = fD(S3, I3, R3, D3, beta, r, d, gamma, sigma,mu);
    g3_D = gD(S3, I3, R3, D3, beta, r, d, gamma, sigma,mu);

    S4 = Stemp + Dt*(f3_S) + (B_1(4).*SWinc1 + B_2(4).*SJ10).*g3_S + (B_1(8).*SWinc1.*g3_S + B_1(12).*SWinc1.*g3_S);
    I4 = Itemp + Dt*(f3_I) + (B_1(4).*IWinc1 + B_2(4).*IJ10).*g3_I + (B_1(8).*IWinc1.*g3_I + B_1(12).*IWinc1.*g3_I);
    R4 = Rtemp + Dt*(f3_R) + (B_1(4).*RWinc1 + B_2(4).*RJ10).*g3_R + (B_1(8).*RWinc1.*g3_R + B_1(12).*RWinc1.*g3_R);
    D4 = Dtemp + Dt*(f3_D) + (B_1(4).*DWinc1 + B_2(4).*DJ10).*g3_D + (B_1(8).*DWinc1.*g3_D + B_1(12).*DWinc1.*g3_D);
    %Refer Equation 4.56

    f4_S = fS(S4, I4, R4, D4, beta, r, d, gamma, sigma,mu);
    g4_S = gS(S4, I4, R4, D4, beta, r, d, gamma, sigma,mu);
    f4_I = fI(S4, I4, R4, D4, beta, r, d, gamma, sigma,mu);
    g4_I = gI(S4, I4, R4, D4, beta, r, d, gamma, sigma,mu);
    f4_R = fR(S4, I4, R4, D4, beta, r, d, gamma, sigma,mu);
    g4_R = gR(S4, I4, R4, D4, beta, r, d, gamma, sigma,mu);
    f4_D = fD(S4, I4, R4, D4, beta, r, d, gamma, sigma,mu);
    g4_D = gD(S4, I4, R4, D4, beta, r, d, gamma, sigma,mu);

    Stemp = Stemp + Dt.*((1/6)*f1_S + (1/3).*f2_S + (1/3).*f3_S + (1/6).*f4_S) + (G_1(1).*SWinc1 + G_2(1).*SJ10).*g1_S + (G_1(2).*SWinc1 + G_2(2).*SJ10).*g2_S + (G_1(3).*SWinc1 + G_2(3).*SJ10).*g3_S + (G_1(4).*SWinc1 + G_2(4).*SJ10).*g4_S;
    Itemp = Itemp + Dt.*((1/6)*f1_I + (1/3).*f2_I + (1/3).*f3_I + (1/6).*f4_I) + (G_1(1).*IWinc1 + G_2(1).*IJ10).*g1_I + (G_1(2).*IWinc1 + G_2(2).*IJ10).*g2_I + (G_1(3).*IWinc1 + G_2(3).*IJ10).*g3_I + (G_1(4).*IWinc1 + G_2(4).*IJ10).*g4_I;
    Rtemp = Rtemp + Dt.*((1/6)*f1_R + (1/3).*f2_R + (1/3).*f3_R + (1/6).*f4_R) + (G_1(1).*RWinc1 + G_2(1).*RJ10).*g1_R + (G_1(2).*RWinc1 + G_2(2).*RJ10).*g2_R + (G_1(3).*RWinc1 + G_2(3).*RJ10).*g3_R + (G_1(4).*RWinc1 + G_2(4).*RJ10).*g4_R;
    Dtemp = Dtemp + Dt.*((1/6)*f1_D + (1/3).*f2_D + (1/3).*f3_D + (1/6).*f4_D) + (G_1(1).*DWinc1 + G_2(1).*DJ10).*g1_D + (G_1(2).*DWinc1 + G_2(2).*DJ10).*g2_D + (G_1(3).*DWinc1 + G_2(3).*DJ10).*g3_D + (G_1(4).*DWinc1 + G_2(4).*DJ10).*g4_D;

    %Refer Equation 4.56
    Ssrk4_2(:,j) = Stemp.';
    Isrk4_2(:,j) = Itemp.';
    Rsrk4_2(:,j) = Rtemp.';
    Dsrk4_2(:,j) = Dtemp.';
end
for i = 1:L
   Ssrk_3 = Ssrk4_2(:,i);
   Isrk_3 = Isrk4_2(:,i);
   Rsrk_3 = Rsrk4_2(:,i);
   Dsrk_3 = Dsrk4_2(:,i);
   Smodel(1,i) = interp1(linspace(1/numel(Ssrk_3),1,numel(Ssrk_3)), sort(Ssrk_3), 0.25);
   Imodel(1,i) = interp1(linspace(1/numel(Isrk_3),1,numel(Isrk_3)), sort(Isrk_3), 0.25);
   Rmodel(1,i) = interp1(linspace(1/numel(Rsrk_3),1,numel(Rsrk_3)), sort(Rsrk_3), 0.25);
   Dmodel(1,i) = interp1(linspace(1/numel(Dsrk_3),1,numel(Dsrk_3)), sort(Dsrk_3), 0.25);
   Smodel(2,i) = interp1(linspace(1/numel(Ssrk_3),1,numel(Ssrk_3)), sort(Ssrk_3), 0.50);
   Imodel(2,i) = interp1(linspace(1/numel(Isrk_3),1,numel(Isrk_3)), sort(Isrk_3), 0.50);
   Rmodel(2,i) = interp1(linspace(1/numel(Rsrk_3),1,numel(Rsrk_3)), sort(Rsrk_3), 0.50);
   Dmodel(2,i) = interp1(linspace(1/numel(Dsrk_3),1,numel(Dsrk_3)), sort(Dsrk_3), 0.50);
   Smodel(3,i) = interp1(linspace(1/numel(Ssrk_3),1,numel(Ssrk_3)), sort(Ssrk_3), 0.75);
   Imodel(3,i) = interp1(linspace(1/numel(Isrk_3),1,numel(Isrk_3)), sort(Isrk_3), 0.75);
   Rmodel(3,i) = interp1(linspace(1/numel(Rsrk_3),1,numel(Rsrk_3)), sort(Rsrk_3), 0.75);
   Dmodel(3,i) = interp1(linspace(1/numel(Dsrk_3),1,numel(Dsrk_3)), sort(Dsrk_3), 0.75);
end
Save_srk4 = Smodel(2,:);
Iave_srk4 = Imodel(2,:);
Rave_srk4 = Rmodel(2,:);
Dave_srk4 = Dmodel(2,:);
%% RK4==================================================
sigma = 0;
%======================================================%SRK-4 Method
Srk4_2 = zeros(1,N);
Irk4_2 = zeros(1,N);
Rrk4_2 = zeros(1,N);
Drk4_2 = zeros(1,N);
Stemp = Szero*ones(1,1);
Itemp = Izero*ones(1,1);
Rtemp = Rzero*ones(1,1);
Dtemp = Dzero*ones(1,1);

for j = 1:L
    beta = Pastparamestimate(j,2);
    r = Pastparamestimate(j,3);
    d = Pastparamestimate(j,4);
    gamma = Pastparamestimate(j,5);

    S1 = Stemp;
    I1 = Itemp;
    R1 = Rtemp;
    D1 = Dtemp;
    %Refer Equation 4.56
    f1_S = fS(S1, I1, R1, D1, beta, r, d, gamma, sigma,mu);
    f1_I = fI(S1, I1, R1, D1, beta, r, d, gamma, sigma,mu);
    f1_R = fR(S1, I1, R1, D1, beta, r, d, gamma, sigma,mu);
    f1_D = fD(S1, I1, R1, D1, beta, r, d, gamma, sigma,mu);
    
    S2 = Stemp + Dt*(A(2)*f1_S); 
    I2 = Itemp + Dt*(A(2)*f1_I);
    R2 = Rtemp + Dt*(A(2)*f1_R); 
    D2 = Dtemp + Dt*(A(2)*f1_D); 

    %Refer Equation 4.56
    f2_S = fS(S2, I2, R2, D2, beta, r, d, gamma, sigma,mu);
    f2_I = fI(S2, I2, R2, D2, beta, r, d, gamma, sigma,mu);
    f2_R = fR(S2, I2, R2, D2, beta, r, d, gamma, sigma,mu);
    f2_D = fD(S2, I2, R2, D2, beta, r, d, gamma, sigma,mu);
    
    S3 = Stemp + Dt*(A(7)*f2_S);
    I3 = Itemp + Dt*(A(7)*f2_I);
    R3 = Rtemp + Dt*(A(7)*f2_R);
    D3 = Dtemp + Dt*(A(7)*f2_D);

    %Refer Equation 4.56
    f3_S = fS(S3, I3, R3, D3, beta, r, d, gamma, sigma,mu);
    f3_I = fI(S3, I3, R3, D3, beta, r, d, gamma, sigma,mu);
    f3_R = fR(S3, I3, R3, D3, beta, r, d, gamma, sigma,mu);
    f3_D = fD(S3, I3, R3, D3, beta, r, d, gamma, sigma,mu);
    
    S4 = Stemp + Dt*(f3_S);
    I4 = Itemp + Dt*(f3_I);
    R4 = Rtemp + Dt*(f3_R);
    D4 = Dtemp + Dt*(f3_D);
    %Refer Equation 4.56

    f4_S = fS(S4, I4, R4, D4, beta, r, d, gamma, sigma,mu);
    g4_S = gS(S4, I4, R4, D4, beta, r, d, gamma, sigma,mu);
    f4_I = fI(S4, I4, R4, D4, beta, r, d, gamma, sigma,mu);
    g4_I = gI(S4, I4, R4, D4, beta, r, d, gamma, sigma,mu);
    f4_R = fR(S4, I4, R4, D4, beta, r, d, gamma, sigma,mu);
    g4_R = gR(S4, I4, R4, D4, beta, r, d, gamma, sigma,mu);
    f4_D = fD(S4, I4, R4, D4, beta, r, d, gamma, sigma,mu);
    g4_D = gD(S4, I4, R4, D4, beta, r, d, gamma, sigma,mu);

    Stemp = Stemp + Dt.*((1/6)*f1_S + (1/3).*f2_S + (1/3).*f3_S + (1/6).*f4_S);
    Itemp = Itemp + Dt.*((1/6)*f1_I + (1/3).*f2_I + (1/3).*f3_I + (1/6).*f4_I);
    Rtemp = Rtemp + Dt.*((1/6)*f1_R + (1/3).*f2_R + (1/3).*f3_R + (1/6).*f4_R);
    Dtemp = Dtemp + Dt.*((1/6)*f1_D + (1/3).*f2_D + (1/3).*f3_D + (1/6).*f4_D);
    
    %Refer Equation 4.56
    Srk4_2(:,j) = Stemp.';
    Irk4_2(:,j) = Itemp.';
    Rrk4_2(:,j) = Rtemp.';
    Drk4_2(:,j) = Dtemp.';
end
% RMSE_srk4x = rmse(Xsrk4_2,Datax);
% RMSE_srk4y = rmse(Ysrk4_2,Datay);
Save_rk4 = mean(Srk4_2,1);
Iave_rk4 = mean(Irk4_2,1);
Rave_rk4 = mean(Rrk4_2,1);
Dave_rk4 = mean(Drk4_2,1);
% Smodel = Save_srk4;
Smodel_2 = Save_rk4;
Imodel_2 = Iave_rk4;
Rmodel_2 = Rave_rk4;
Dmodel_2 = Dave_rk4;
%% Adam_Bashford_Predictor_Corrector=====================
Sabpc_2 = Szero*ones(1,N+1);
Iabpc_2 = Izero*ones(1,N+1);
Rabpc_2 = Rzero*ones(1,N+1);
Dabpc_2 = Dzero*ones(1,N+1);

pcS1 = 0; pcS2 = 0; pcS3 = 0;
pcI1 = 0; pcI2 = 0; pcI3 = 0;
pcR1 = 0; pcR2 = 0; pcR3 = 0;
pcD1 = 0; pcD2 = 0; pcD3 = 0;

for j = 1:L
    beta = Pastparamestimate(j,2);
    r = Pastparamestimate(j,3);
    d = Pastparamestimate(j,4);
    gamma = Pastparamestimate(j,5);

    S1 = Sabpc_2(j);
    I1 = Iabpc_2(j);
    R1 = Rabpc_2(j);
    D1 = Dabpc_2(j);

    f1_S = fS(S1,I1,R1,D1,beta,r,d,gamma,sigma,mu);
    f1_I = fI(S1,I1,R1,D1,beta,r,d,gamma,sigma,mu);
    f1_R = fR(S1,I1,R1,D1,beta,r,d,gamma,sigma,mu);
    f1_D = fD(S1,I1,R1,D1,beta,r,d,gamma,sigma,mu);

    pcS0 = f1_S;
    pcI0 = f1_I;
    pcR0 = f1_R;
    pcD0 = f1_D;

    S2 = Sabpc_2(j) + dt*(0.5*f1_S);
    I2 = Iabpc_2(j) + dt*(0.5*f1_I);
    R2 = Rabpc_2(j) + dt*(0.5*f1_R);
    D2 = Dabpc_2(j) + dt*(0.5*f1_D);

    f2_S = fS(S2,I2,R2,D2,beta,r,d,gamma,sigma,mu);
    f2_I = fI(S2,I2,R2,D2,beta,r,d,gamma,sigma,mu);
    f2_R = fR(S2,I2,R2,D2,beta,r,d,gamma,sigma,mu);
    f2_D = fD(S2,I2,R2,D2,beta,r,d,gamma,sigma,mu);

    S3 = Sabpc_2(j) + dt*(0.5*f2_S);
    I3 = Iabpc_2(j) + dt*(0.5*f2_I);
    R3 = Rabpc_2(j) + dt*(0.5*f2_R);
    D3 = Dabpc_2(j) + dt*(0.5*f2_D);

    f3_S = fS(S3,I3,R3,D3,beta,r,d,gamma,sigma,mu);
    f3_I = fI(S3,I3,R3,D3,beta,r,d,gamma,sigma,mu);
    f3_R = fR(S3,I3,R3,D3,beta,r,d,gamma,sigma,mu);
    f3_D = fD(S3,I3,R3,D3,beta,r,d,gamma,sigma,mu);

    S4 = Sabpc_2(j) + dt*f3_S;
    I4 = Iabpc_2(j) + dt*f3_I;
    R4 = Rabpc_2(j) + dt*f3_R;
    D4 = Dabpc_2(j) + dt*f3_D;

    f4_S = fS(S4,I4,R4,D4,beta,r,d,gamma,sigma,mu);
    f4_I = fI(S4,I4,R4,D4,beta,r,d,gamma,sigma,mu);
    f4_R = fR(S4,I4,R4,D4,beta,r,d,gamma,sigma,mu);
    f4_D = fD(S4,I4,R4,D4,beta,r,d,gamma,sigma,mu);

    Sabpc_2(j+1) = Sabpc_2(j) + dt*((1/6)*f1_S + (1/3)*f2_S + (1/3)*f3_S + (1/6)*f4_S);
    Iabpc_2(j+1) = Iabpc_2(j) + dt*((1/6)*f1_I + (1/3)*f2_I + (1/3)*f3_I + (1/6)*f4_I);
    Rabpc_2(j+1) = Rabpc_2(j) + dt*((1/6)*f1_R + (1/3)*f2_R + (1/3)*f3_R + (1/6)*f4_R);
    Dabpc_2(j+1) = Dabpc_2(j) + dt*((1/6)*f1_D + (1/3)*f2_D + (1/3)*f3_D + (1/6)*f4_D);

    pcS1 = pcS0; pcS2 = pcS1; pcS3 = pcS2;
    pcI1 = pcI0; pcI2 = pcI1; pcI3 = pcI2;
    pcR1 = pcR0; pcR2 = pcR1; pcR3 = pcR2;
    pcD1 = pcD0; pcD2 = pcD1; pcD3 = pcD2;
end
        
for j = 1:L
    beta = Pastparamestimate(j,2);
    r = Pastparamestimate(j,3);
    d = Pastparamestimate(j,4);
    gamma = Pastparamestimate(j,5);
        
    S1 = Sabpc_2(j);
    I1 = Iabpc_2(j);
    R1 = Rabpc_2(j);
    D1 = Dabpc_2(j);

    f1_S = fS(S1,I1,R1,D1,beta,r,d,gamma,sigma,mu);
    f1_I = fI(S1,I1,R1,D1,beta,r,d,gamma,sigma,mu);
    f1_R = fR(S1,I1,R1,D1,beta,r,d,gamma,sigma,mu);
    f1_D = fD(S1,I1,R1,D1,beta,r,d,gamma,sigma,mu);

    pcS0 = f1_S;
    pcI0 = f1_I;
    pcR0 = f1_R;
    pcD0 = f1_D;

    wS = Sabpc_2(j) + dt*( 55.0 * pcS0 - 59.0 * pcS1 + 37.0 * pcS2 - 9.0 * pcS3 ) / 24.0;
    wI = Iabpc_2(j) + dt*( 55.0 * pcI0 - 59.0 * pcI1 + 37.0 * pcI2 - 9.0 * pcI3 ) / 24.0;
    wR = Rabpc_2(j) + dt*( 55.0 * pcR0 - 59.0 * pcR1 + 37.0 * pcR2 - 9.0 * pcR3 ) / 24.0;
    wD = Dabpc_2(j) + dt*( 55.0 * pcD0 - 59.0 * pcD1 + 37.0 * pcD2 - 9.0 * pcD3 ) / 24.0;

    fw_S = fS(wS,wI,wR,wD,beta,r,d,gamma,sigma,mu);
    fw_I = fI(wS,wI,wR,wD,beta,r,d,gamma,sigma,mu);
    fw_R = fR(wS,wI,wR,wD,beta,r,d,gamma,sigma,mu);
    fw_D = fD(wS,wI,wR,wD,beta,r,d,gamma,sigma,mu);
    
    Sabpc_2(j+1) = Sabpc_2(j) + dt*( 9.0 * fw_S + 19.0 * pcS0 - 5.0 * pcS1 + pcS2 ) / 24.0;
    Iabpc_2(j+1) = Iabpc_2(j) + dt*( 9.0 * fw_I + 19.0 * pcI0 - 5.0 * pcI1 + pcI2 ) / 24.0;
    Rabpc_2(j+1) = Rabpc_2(j) + dt*( 9.0 * fw_R + 19.0 * pcR0 - 5.0 * pcR1 + pcR2 ) / 24.0;
    Dabpc_2(j+1) = Dabpc_2(j) + dt*( 9.0 * fw_D + 19.0 * pcD0 - 5.0 * pcD1 + pcD2 ) / 24.0;

    pcS1 = pcS0; pcS2 = pcS1; pcS3 = pcS2;
    pcI1 = pcI0; pcI2 = pcI1; pcI3 = pcI2;
    pcR1 = pcR0; pcR2 = pcR1; pcR3 = pcR2;
    pcD1 = pcD0; pcD2 = pcD1; pcD3 = pcD2;
end
% RMSE_srk4x = rmse(Xsrk4_2,Datax);
% RMSE_srk4y = rmse(Ysrk4_2,Datay);
Save_abpc = mean(Sabpc_2,1);
Iave_abpc = mean(Iabpc_2,1);
Rave_abpc = mean(Rabpc_2,1);
Dave_abpc = mean(Dabpc_2,1);
% Smodel = Save_srk4;
SDet = Save_abpc(2:end);
IDet = Iave_abpc(2:end);
RDet = Rave_abpc(2:end);
DDet = Dave_abpc(2:end);
%% Plot
figure('Position', [1, 1, 600, 1000])
tlayout = tiledlayout(4,1);
tlayout.TileSpacing = 'compact';
tlayout.Padding = 'compact';
    ax1 = nexttile;
    plot(tdays,Sdata,"-s",'Color',"#000000",'MarkerIndices',1:20:length(tdays)); hold on;
%     plot(tdays,Smodel_2,':', 'Color',"#2F4F4F"); hold on;
    plot(tdays,SDet,"--", 'Color',"#778899");
%     plot(tdays,Smodel(2,:),'-', 'Color',"#696969"); hold on;
    cax1 = [105/255, 105/255, 105/255];
%     fill([tdays, fliplr(tdays)], [Smodel(1,:), fliplr(Smodel(3,:))], cax1, 'LineStyle', "None",'FaceAlpha',0.2); hold on;
    title(ax1,'Susceptible Population, S')
    ax1.FontSize = 8;
    hold off;
    legend('Susceptible Data','Susceptible Simulation','Location','best')
    xlabel('(e)', 'HorizontalAlignment', 'center'); 
    text(max(xlim), min(ylim) - 0.1 * range(ylim), 'date (days)', 'HorizontalAlignment', 'right');
%     ylabel('Population')    
    
    ax2 = nexttile; 
    plot(tdays,Idata,"-s",'Color',"#2C00D2",'MarkerIndices',1:20:length(tdays)); hold on;
%     plot(tdays,Imodel_2,':','Color',"#273DE3"); hold on;
    plot(tdays,IDet,'--', 'Color',"#4BA3E5")
%     plot(tdays,Imodel(2,:),'-','Color',"#77DFE9"); hold on;
    cax2 = [39/255,61/255,227/255];
%     fill([tdays, fliplr(tdays)], [Imodel(1,:), fliplr(Imodel(3,:))], cax2, 'LineStyle', "None",'FaceAlpha',0.2); hold on;
    title(ax2,'Infected Population')
    ax2.FontSize = 8;
    legend('Infected Data','Infected Simulation','Location','best')
    hold off;
    xlabel('(f)', 'HorizontalAlignment', 'center'); 
    text(max(xlim), min(ylim) - 0.1 * range(ylim), 'date (days)', 'HorizontalAlignment', 'right');
    %ylabel('Infected Population')    
    
    ax3 = nexttile;
    plot(tdays,Rdata,'-s','Color',"#008080",'MarkerIndices',1:20:length(tdays)); hold on;
%     plot(tdays,Rmodel_2,':','Color',"#556B2F");hold on;
    plot(tdays,RDet,'--', 'Color',"#3CB371")
%     plot(tdays,Rmodel(2,:),'-','Color',"#00FF00");hold on;
    title(ax3,'Recovered Population, R')
    cax3 = [70/255,192/255,30/255];
%     fill([tdays, fliplr(tdays)], [Rmodel(1,:), fliplr(Rmodel(3,:))], cax3, 'LineStyle', "None",'FaceAlpha',0.2); hold on;
    ax3.FontSize = 8;
    legend('Recovered Data','Recovered Simulation','Location','best')
    hold off;
    xlabel('(g)', 'HorizontalAlignment', 'center'); 
    text(max(xlim), min(ylim) - 0.1 * range(ylim), 'date (days)', 'HorizontalAlignment', 'right');
    %ylabel('Recovered Population')    
    
    ax4 = nexttile;
    plot(tdays,Ddata,'-s','Color',"#8B0000",'MarkerIndices',1:20:length(tdays)); hold on;
%     plot(tdays,Dmodel_2,':','Color',"#FF0000");hold on;
    plot(tdays,DDet,'--', 'Color',"#8B4513")
%     plot(tdays,Dmodel(2,:),'-','Color',"#F08080");hold on;
    title(ax4,'Death Population, D')
    cax4 = [217/255,46/255,157/255];
%     fill([tdays, fliplr(tdays)], [Dmodel(1,:), fliplr(Dmodel(3,:))], cax4, 'LineStyle', "None",'FaceAlpha',0.2); hold on;
    ax4.FontSize = 8;
    legend('Death Data','Death Simulation','Location','best')
    hold off;
    xlabel('(h)', 'HorizontalAlignment', 'center'); 
    text(max(xlim), min(ylim) - 0.1 * range(ylim), 'date (days)', 'HorizontalAlignment', 'right');
    %ylabel('Death Population') 
    
    print -depsc sirdMalaysia.eps
%% RMSE
RMSE = zeros(3,4);
RMSE(3,1) = sqrt(sum((Save_srk4(:) - Sdata(:)).^2)/length(Save_srk4));
RMSE(1,1) = sqrt(sum((Smodel_2(:) - Sdata(:)).^2)/length(Save_srk4));
RMSE(2,1) = sqrt(sum((SDet(:) - Sdata(:)).^2)/length(Save_srk4));
RMSE(3,2) = sqrt(sum((Iave_srk4(:) - Idata(:)).^2)/length(Save_srk4));
RMSE(1,2) = sqrt(sum((Imodel_2(:) - Idata(:)).^2)/length(Save_srk4));
RMSE(2,2) = sqrt(sum((IDet(:) - Idata(:)).^2)/length(Save_srk4));
RMSE(3,3) = sqrt(sum((Rave_srk4(:) - Rdata(:)).^2)/length(Save_srk4));
RMSE(1,3) = sqrt(sum((Rmodel_2(:) - Rdata(:)).^2)/length(Save_srk4));
RMSE(2,3) = sqrt(sum((RDet(:) - Rdata(:)).^2)/length(Save_srk4));
RMSE(3,4) = sqrt(sum((Dave_srk4(:) - Ddata(:)).^2)/length(Save_srk4));
RMSE(1,4) = sqrt(sum((Dmodel_2(:) - Ddata(:)).^2)/length(Save_srk4));
RMSE(2,4) = sqrt(sum((DDet(:) - Ddata(:)).^2)/length(Save_srk4));

