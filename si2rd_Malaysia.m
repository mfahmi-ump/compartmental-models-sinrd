clear model data params options
datafile = readmatrix('si3rd_data_Malaysia.csv');
Ptotal = 3.153*10^7;%malaysia
rng = 3683;
I1data = datafile(11:end,3);
I2data = datafile(11:end,5)+ datafile(11:end,6);
Rdata = datafile(11:end,2);
Ddata = datafile(11:end,4);
Sdata = datafile(11:end,7);
tdata = datafile(11:end,1);
t = 1:1:length(I1data);
tdays = datetime(2020,3,24) + caldays(1:length(t));
Pastparamestimate = readmatrix('Malaysia-SI2RD-param-estimate.csv');
Szero = Sdata(1);
I1zero = I1data(1);
I2zero = I2data(1);
Rzero = Rdata(1);
Dzero = Ddata(1);
T = length(t);
N = T;
M = 500;

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
fS  = @(S,I1,I2,R,D,Beta1,Beta2,r1,r2,d1,d2,p1,p2,gamma,sigma,mu) mu*(S+I1+I2+R)-S./Ptotal.*(Beta1.*I1+Beta2.*I2)+gamma.*R - mu.*S;
gS  = @(S,I1,I2,R,D,Beta1,Beta2,r1,r2,d1,d2,p1,p2,gamma,sigma,mu) -sigma.*I1.*S./Ptotal;
fI1 = @(S,I1,I2,R,D,Beta1,Beta2,r1,r2,d1,d2,p1,p2,gamma,sigma,mu) p1.*S./Ptotal.*(Beta1.*I1+Beta2.*I2)-(r1+d1).*I1 - mu.*I1;
gI1 = @(S,I1,I2,R,D,Beta1,Beta2,r1,r2,d1,d2,p1,p2,gamma,sigma,mu) p1.*sigma.*I1.*S./Ptotal;
fI2 = @(S,I1,I2,R,D,Beta1,Beta2,r1,r2,d1,d2,p1,p2,gamma,sigma,mu) p2.*S./Ptotal.*(Beta1.*I1+Beta2.*I2)-(r2+d2).*I2 - mu.*I2;
gI2 = @(S,I1,I2,R,D,Beta1,Beta2,r1,r2,d1,d2,p1,p2,gamma,sigma,mu) p2.*sigma.*I1.*S./Ptotal;
fR  = @(S,I1,I2,R,D,Beta1,Beta2,r1,r2,d1,d2,p1,p2,gamma,sigma,mu) r1.*I1 + r2.*I2 - gamma.*R - mu.*R;
gR  = @(S,I1,I2,R,D,Beta1,Beta2,r1,r2,d1,d2,p1,p2,gamma,sigma,mu) 0;
fD  = @(S,I1,I2,R,D,Beta1,Beta2,r1,r2,d1,d2,p1,p2,gamma,sigma,mu) d1.*I1 + d2.*I2;
gD  = @(S,I1,I2,R,D,Beta1,Beta2,r1,r2,d1,d2,p1,p2,gamma,sigma,mu) 0;
%% SRK4==================================================
sigma = 0.05;
%======================================================%SRK-4 Method
Ssrk4_2 = zeros(M,N);
I1srk4_2 = zeros(M,N);
I2srk4_2 = zeros(M,N);
Rsrk4_2 = zeros(M,N);
Dsrk4_2 = zeros(M,N);
Stemp = Szero*ones(M,1);
I1temp = I1zero*ones(M,1);
I2temp = I2zero*ones(M,1);
Rtemp = Rzero*ones(M,1);
Dtemp = Dzero*ones(M,1);
Smodel = ones(3,N);
I1model = ones(3,N);
I2model = ones(3,N);
Rmodel = ones(3,N);
Dmodel = ones(3,N);
SdW1 = sqrt(dt)*randn(M,N);       %Numerical Brownian increments
SdW2 = sqrt(dt)*randn(M,N);       %Numerical Brownian increments
I1dW1 = sqrt(dt)*randn(M,N);       %Numerical Brownian increments
I1dW2 = sqrt(dt)*randn(M,N);       %Numerical Brownian increments
I2dW1 = sqrt(dt)*randn(M,N);       %Numerical Brownian increments
I2dW2 = sqrt(dt)*randn(M,N);       %Numerical Brownian increments
RdW1 = sqrt(dt)*randn(M,N);       %Numerical Brownian increments
RdW2 = sqrt(dt)*randn(M,N);       %Numerical Brownian increments
DdW1 = sqrt(dt)*randn(M,N);       %Numerical Brownian increments
DdW2 = sqrt(dt)*randn(M,N);       %Numerical Brownian increments

for j = 1:L
    Beta1 = Pastparamestimate(j,2);
    Beta2 = Pastparamestimate(j,3);
    r1 = Pastparamestimate(j,4);
    r2 = Pastparamestimate(j,5);
    d1 = Pastparamestimate(j,6);
    d2 = Pastparamestimate(j,7);
    p1 = Pastparamestimate(j,8);
    p2 = Pastparamestimate(j,9);
    gamma = Pastparamestimate(j,10);

    SWinc1 = sum(SdW1(:,R*(j-1)+1:R*j),2);
    SWinc2 = sum(SdW2(:,R*(j-1)+1:R*j),2);
    I1Winc1 = sum(I1dW1(:,R*(j-1)+1:R*j),2);
    I1Winc2 = sum(I1dW2(:,R*(j-1)+1:R*j),2);
    I2Winc1 = sum(I2dW1(:,R*(j-1)+1:R*j),2);
    I2Winc2 = sum(I2dW2(:,R*(j-1)+1:R*j),2);
    RWinc1 = sum(RdW1(:,R*(j-1)+1:R*j),2);
    RWinc2 = sum(RdW2(:,R*(j-1)+1:R*j),2);
    DWinc1 = sum(DdW1(:,R*(j-1)+1:R*j),2);
    DWinc2 = sum(DdW2(:,R*(j-1)+1:R*j),2);

    SJ10 = (0.5*(Dt^(1/2))*(SWinc1+(1/sqrt(3)).*SWinc2));
    I1J10 = (0.5*(Dt^(1/2))*(I1Winc1+(1/sqrt(3)).*I1Winc2));
    I2J10 = (0.5*(Dt^(1/2))*(I2Winc1+(1/sqrt(3)).*I2Winc2));
    RJ10 = (0.5*(Dt^(1/2))*(RWinc1+(1/sqrt(3)).*RWinc2));
    DJ10 = (0.5*(Dt^(1/2))*(DWinc1+(1/sqrt(3)).*DWinc2));

    S1 = Stemp;
    I11 = I1temp;
    I21 = I2temp;
    R1 = Rtemp;
    D1 = Dtemp;
    %Refer Equation 4.56
    f1_S = fS(S1, I11, I21, R1, D1, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    g1_S = gS(S1, I11, I21, R1, D1, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    f1_I1 = fI1(S1, I11, I21, R1, D1, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    g1_I1 = gI1(S1, I11, I21, R1, D1, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    f1_I2 = fI2(S1, I11, I21, R1, D1, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    g1_I2 = gI2(S1, I11, I21, R1, D1, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    f1_R = fR(S1, I11, I21, R1, D1, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    g1_R = gR(S1, I11, I21, R1, D1, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    f1_D = fD(S1, I11, I21, R1, D1, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    g1_D = gD(S1, I11, I21, R1, D1, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);


    S2 = Stemp + Dt*(A(2)*f1_S) + (B_1(2).*SWinc1 + B_2(2).*SJ10).*g1_S;
    I12 = I1temp + Dt*(A(2)*f1_I1) + (B_1(2).*I1Winc1 + B_2(2).*I1J10).*g1_I1;
    I22 = I2temp + Dt*(A(2)*f1_I2) + (B_1(2).*I2Winc1 + B_2(2).*I2J10).*g1_I2;
    R2 = Rtemp + Dt*(A(2)*f1_R) + (B_1(2).*RWinc1 + B_2(2).*RJ10).*g1_R;
    D2 = Dtemp + Dt*(A(2)*f1_D) + (B_1(2).*DWinc1 + B_2(2).*DJ10).*g1_D;

    %Refer Equation 4.56
    f2_S = fS(S2, I12, I22, R2, D2, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    g2_S = gS(S2, I12, I22, R2, D2, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    f2_I1 = fI1(S2, I12, I22, R2, D2, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    g2_I1 = gI1(S2, I12, I22, R2, D2, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    f2_I2 = fI2(S2, I12, I22, R2, D2, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    g2_I2 = gI2(S2, I12, I22, R2, D2, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    f2_R = fR(S2, I12, I22, R2, D2, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    g2_R = gR(S2, I12, I22, R2, D2, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    f2_D = fD(S2, I12, I22, R2, D2, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    g2_D = gD(S2, I12, I22, R2, D2, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);

    S3 = Stemp + Dt*(A(7)*f2_S) + (B_1(3).*SWinc1 + B_2(3).*SJ10).*g2_S + B_1(7)*SJ10.*g2_S;
    I13 = I1temp + Dt*(A(7)*f2_I1) + (B_1(3).*I1Winc1 + B_2(3).*I1J10).*g2_I1 + B_1(7)*I1J10.*g2_I1;
    I23 = I2temp + Dt*(A(7)*f2_I2) + (B_1(3).*I2Winc1 + B_2(3).*I2J10).*g2_I2 + B_1(7)*I2J10.*g2_I2;
    R3 = Rtemp + Dt*(A(7)*f2_R) + (B_1(3).*RWinc1 + B_2(3).*RJ10).*g2_R + B_1(7)*RJ10.*g2_R;
    D3 = Dtemp + Dt*(A(7)*f2_D) + (B_1(3).*DWinc1 + B_2(3).*DJ10).*g2_D + B_1(7)*DJ10.*g2_D;

    %Refer Equation 4.56
    f3_S = fS(S3, I13, I23, R3, D3, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    g3_S = gS(S3, I13, I23, R3, D3, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    f3_I1 = fI1(S3, I13, I23, R3, D3, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    g3_I1 = gI1(S3, I13, I23, R3, D3, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    f3_I2 = fI2(S3, I13, I23, R3, D3, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    g3_I2 = gI2(S3, I13, I23, R3, D3, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    f3_R = fR(S3, I13, I23, R3, D3, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    g3_R = gR(S3, I13, I23, R3, D3, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    f3_D = fD(S3, I13, I23, R3, D3, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    g3_D = gD(S3, I13, I23, R3, D3, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);

    
    S4 = Stemp + Dt*(f3_S) + (B_1(4).*SWinc1 + B_2(4).*SJ10).*g3_S + (B_1(8).*SWinc1.*g3_S + B_1(12).*SWinc1.*g3_S);
    I14 = I1temp + Dt*(f3_I1) + (B_1(4).*I1Winc1 + B_2(4).*I1J10).*g3_I1 + (B_1(8).*I1Winc1.*g3_I1 + B_1(12).*I1Winc1.*g3_I1);
    I24 = I2temp + Dt*(f3_I2) + (B_1(4).*I2Winc1 + B_2(4).*I2J10).*g3_I2 + (B_1(8).*I2Winc1.*g3_I2 + B_1(12).*I2Winc1.*g3_I2);
    R4 = Rtemp + Dt*(f3_R) + (B_1(4).*RWinc1 + B_2(4).*RJ10).*g3_R + (B_1(8).*RWinc1.*g3_R + B_1(12).*RWinc1.*g3_R);
    D4 = Dtemp + Dt*(f3_D) + (B_1(4).*DWinc1 + B_2(4).*DJ10).*g3_D + (B_1(8).*DWinc1.*g3_D + B_1(12).*DWinc1.*g3_D);
    %Refer Equation 4.56

    f4_S = fS(S4, I14, I24, R3, D3, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    g4_S = gS(S4, I14, I24, R3, D3, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    f4_I1 = fI1(S4, I14, I24, R3, D3, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    g4_I1 = gI1(S4, I14, I24, R3, D3, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    f4_I2 = fI2(S4, I14, I24, R3, D3, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    g4_I2 = gI2(S4, I14, I24, R3, D3, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    f4_R = fR(S4, I14, I24, R3, D3, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    g4_R = gR(S4, I14, I24, R3, D3, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    f4_D = fD(S4, I14, I24, R3, D3, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    g4_D = gD(S4, I14, I24, R3, D3, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);


    Stemp = Stemp + Dt.*((1/6)*f1_S + (1/3).*f2_S + (1/3).*f3_S + (1/6).*f4_S) + (G_1(1).*SWinc1 + G_2(1).*SJ10).*g1_S + (G_1(2).*SWinc1 + G_2(2).*SJ10).*g2_S + (G_1(3).*SWinc1 + G_2(3).*SJ10).*g3_S + (G_1(4).*SWinc1 + G_2(4).*SJ10).*g4_S;
    I1temp = I1temp + Dt.*((1/6)*f1_I1 + (1/3).*f2_I1 + (1/3).*f3_I1 + (1/6).*f4_I1) + (G_1(1).*I1Winc1 + G_2(1).*I1J10).*g1_I1 + (G_1(2).*I1Winc1 + G_2(2).*I1J10).*g2_I1 + (G_1(3).*I1Winc1 + G_2(3).*I1J10).*g3_I1 + (G_1(4).*I1Winc1 + G_2(4).*I1J10).*g4_I1;
    I2temp = I2temp + Dt.*((1/6)*f1_I2 + (1/3).*f2_I2 + (1/3).*f3_I2 + (1/6).*f4_I2) + (G_1(1).*I2Winc1 + G_2(1).*I2J10).*g1_I2 + (G_1(2).*I2Winc1 + G_2(2).*I2J10).*g2_I2 + (G_1(3).*I2Winc1 + G_2(3).*I2J10).*g3_I2 + (G_1(4).*I2Winc1 + G_2(4).*I2J10).*g4_I2;
    Rtemp = Rtemp + Dt.*((1/6)*f1_R + (1/3).*f2_R + (1/3).*f3_R + (1/6).*f4_R) + (G_1(1).*RWinc1 + G_2(1).*RJ10).*g1_R + (G_1(2).*RWinc1 + G_2(2).*RJ10).*g2_R + (G_1(3).*RWinc1 + G_2(3).*RJ10).*g3_R + (G_1(4).*RWinc1 + G_2(4).*RJ10).*g4_R;
    Dtemp = Dtemp + Dt.*((1/6)*f1_D + (1/3).*f2_D + (1/3).*f3_D + (1/6).*f4_D) + (G_1(1).*DWinc1 + G_2(1).*DJ10).*g1_D + (G_1(2).*DWinc1 + G_2(2).*DJ10).*g2_D + (G_1(3).*DWinc1 + G_2(3).*DJ10).*g3_D + (G_1(4).*DWinc1 + G_2(4).*DJ10).*g4_D;

    %Refer Equation 4.56
    Ssrk4_2(:,j) = Stemp.';
    I1srk4_2(:,j) = I1temp.';
    I2srk4_2(:,j) = I2temp.';
    Rsrk4_2(:,j) = Rtemp.';
    Dsrk4_2(:,j) = Dtemp.';
end
for i = 1:L
   Ssrk_3 = Ssrk4_2(:,i);
   I1srk_3 = I1srk4_2(:,i);
   I2srk_3 = I2srk4_2(:,i);
   Rsrk_3 = Rsrk4_2(:,i);
   Dsrk_3 = Dsrk4_2(:,i);
   Smodel(1,i) = interp1(linspace(1/numel(Ssrk_3),1,numel(Ssrk_3)), sort(Ssrk_3), 0.25);
   I1model(1,i) = interp1(linspace(1/numel(I1srk_3),1,numel(I1srk_3)), sort(I1srk_3), 0.25);
   I2model(1,i) = interp1(linspace(1/numel(I2srk_3),1,numel(I2srk_3)), sort(I2srk_3), 0.25);
   Rmodel(1,i) = interp1(linspace(1/numel(Rsrk_3),1,numel(Rsrk_3)), sort(Rsrk_3), 0.25);
   Dmodel(1,i) = interp1(linspace(1/numel(Dsrk_3),1,numel(Dsrk_3)), sort(Dsrk_3), 0.25);
   Smodel(2,i) = interp1(linspace(1/numel(Ssrk_3),1,numel(Ssrk_3)), sort(Ssrk_3), 0.50);
   I1model(2,i) = interp1(linspace(1/numel(I1srk_3),1,numel(I1srk_3)), sort(I1srk_3), 0.50);
   I2model(2,i) = interp1(linspace(1/numel(I2srk_3),1,numel(I2srk_3)), sort(I2srk_3), 0.50);
   Rmodel(2,i) = interp1(linspace(1/numel(Rsrk_3),1,numel(Rsrk_3)), sort(Rsrk_3), 0.50);
   Dmodel(2,i) = interp1(linspace(1/numel(Dsrk_3),1,numel(Dsrk_3)), sort(Dsrk_3), 0.50);
   Smodel(3,i) = interp1(linspace(1/numel(Ssrk_3),1,numel(Ssrk_3)), sort(Ssrk_3), 0.75);
   I1model(3,i) = interp1(linspace(1/numel(I1srk_3),1,numel(I1srk_3)), sort(I1srk_3), 0.75);
   I2model(3,i) = interp1(linspace(1/numel(I2srk_3),1,numel(I2srk_3)), sort(I2srk_3), 0.75);
   Rmodel(3,i) = interp1(linspace(1/numel(Rsrk_3),1,numel(Rsrk_3)), sort(Rsrk_3), 0.75);
   Dmodel(3,i) = interp1(linspace(1/numel(Dsrk_3),1,numel(Dsrk_3)), sort(Dsrk_3), 0.75);
end
Save_srk4 = Smodel(2,:);
I1ave_srk4 = I1model(2,:);
I2ave_srk4 = I2model(2,:);
Rave_srk4 = Rmodel(2,:);
Dave_srk4 = Dmodel(2,:);
%% RK4==================================================
sigma = 0;
Srk4_2 = zeros(1,N);
I1rk4_2 = zeros(1,N);
I2rk4_2 = zeros(1,N);
Rrk4_2 = zeros(1,N);
Drk4_2 = zeros(1,N);
Stemp = Szero*ones(1,1);
I1temp = I1zero*ones(1,1);
I2temp = I2zero*ones(1,1);
Rtemp = Rzero*ones(1,1);
Dtemp = Dzero*ones(1,1);

for j = 1:L
    Beta1 = Pastparamestimate(j,2);
    Beta2 = Pastparamestimate(j,3);
    r1 = Pastparamestimate(j,4);
    r2 = Pastparamestimate(j,5);
    d1 = Pastparamestimate(j,6);
    d2 = Pastparamestimate(j,7);
    p1_1 = Pastparamestimate(j,8);
    p2_1 = Pastparamestimate(j,9);
    gamma = Pastparamestimate(j,10);
    
    p1 = p1_1/(p1_1+p2_1);
    p2 = p2_1/(p1_1+p2_1);

    S1 = Stemp;
    I11 = I1temp;
    I21 = I2temp;
    R1 = Rtemp;
    D1 = Dtemp;
    %Refer Equation 4.56
    f1_S = fS(S1, I11, I21, R1, D1, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    g1_S = gS(S1, I11, I21, R1, D1, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    f1_I1 = fI1(S1, I11, I21, R1, D1, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    g1_I1 = gI1(S1, I11, I21, R1, D1, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    f1_I2 = fI2(S1, I11, I21, R1, D1, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    g1_I2 = gI2(S1, I11, I21, R1, D1, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    f1_R = fR(S1, I11, I21, R1, D1, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    g1_R = gR(S1, I11, I21, R1, D1, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    f1_D = fD(S1, I11, I21, R1, D1, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    g1_D = gD(S1, I11, I21, R1, D1, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);

    S2 = Stemp + Dt*(A(2)*f1_S);
    I12 = I1temp + Dt*(A(2)*f1_I1);
    I22 = I2temp + Dt*(A(2)*f1_I2);
    R2 = Rtemp + Dt*(A(2)*f1_R);
    D2 = Dtemp + Dt*(A(2)*f1_D);

    %Refer Equation 4.56
    f2_S = fS(S2, I12, I22, R2, D2, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    g2_S = gS(S2, I12, I22, R2, D2, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    f2_I1 = fI1(S2, I12, I22, R2, D2, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    g2_I1 = gI1(S2, I12, I22, R2, D2, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    f2_I2 = fI2(S2, I12, I22, R2, D2, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    g2_I2 = gI2(S2, I12, I22, R2, D2, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    f2_R = fR(S2, I12, I22, R2, D2, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    g2_R = gR(S2, I12, I22, R2, D2, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    f2_D = fD(S2, I12, I22, R2, D2, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    g2_D = gD(S2, I12, I22, R2, D2, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);

    S3 = Stemp + Dt*(A(7)*f2_S);
    I13 = I1temp + Dt*(A(7)*f2_I1);
    I23 = I2temp + Dt*(A(7)*f2_I2);
    R3 = Rtemp + Dt*(A(7)*f2_R);
    D3 = Dtemp + Dt*(A(7)*f2_D);

    %Refer Equation 4.56
    f3_S = fS(S3, I13, I23, R3, D3, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    g3_S = gS(S3, I13, I23, R3, D3, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    f3_I1 = fI1(S3, I13, I23, R3, D3, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    g3_I1 = gI1(S3, I13, I23, R3, D3, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    f3_I2 = fI2(S3, I13, I23, R3, D3, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    g3_I2 = gI2(S3, I13, I23, R3, D3, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    f3_R = fR(S3, I13, I23, R3, D3, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    g3_R = gR(S3, I13, I23, R3, D3, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    f3_D = fD(S3, I13, I23, R3, D3, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    g3_D = gD(S3, I13, I23, R3, D3, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);

    S4 = Stemp + Dt*(f3_S);
    I14 = I1temp + Dt*(f3_I1);
    I24 = I2temp + Dt*(f3_I2);
    R4 = Rtemp + Dt*(f3_R);
    D4 = Dtemp + Dt*(f3_D);
    %Refer Equation 4.56

    f4_S = fS(S4, I14, I24, R3, D3, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    g4_S = gS(S4, I14, I24, R3, D3, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    f4_I1 = fI1(S4, I14, I24, R3, D3, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    g4_I1 = gI1(S4, I14, I24, R3, D3, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    f4_I2 = fI2(S4, I14, I24, R3, D3, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    g4_I2 = gI2(S4, I14, I24, R3, D3, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    f4_R = fR(S4, I14, I24, R3, D3, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    g4_R = gR(S4, I14, I24, R3, D3, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    f4_D = fD(S4, I14, I24, R3, D3, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);
    g4_D = gD(S4, I14, I24, R3, D3, Beta1, Beta2, r1, r2, d1, d2, p1, p2, gamma, sigma, mu);

    Stemp = Stemp + Dt.*((1/6)*f1_S + (1/3).*f2_S + (1/3).*f3_S + (1/6).*f4_S);
    I1temp = I1temp + Dt.*((1/6)*f1_I1 + (1/3).*f2_I1 + (1/3).*f3_I1 + (1/6).*f4_I1);
    I2temp = I2temp + Dt.*((1/6)*f1_I2 + (1/3).*f2_I2 + (1/3).*f3_I2 + (1/6).*f4_I2);
    Rtemp = Rtemp + Dt.*((1/6)*f1_R + (1/3).*f2_R + (1/3).*f3_R + (1/6).*f4_R);
    Dtemp = Dtemp + Dt.*((1/6)*f1_D + (1/3).*f2_D + (1/3).*f3_D + (1/6).*f4_D);
    
    %Refer Equation 4.56
    Srk4_2(:,j) = Stemp.';
    I1rk4_2(:,j) = I1temp.';
    I2rk4_2(:,j) = I2temp.';
    Rrk4_2(:,j) = Rtemp.';
    Drk4_2(:,j) = Dtemp.';
end
% RMSE_srk4x = rmse(Xsrk4_2,Datax);
% RMSE_srk4y = rmse(Ysrk4_2,Datay);
Smodel_2 = mean(Srk4_2,1);
I1model_2 = mean(I1rk4_2,1);
I2model_2 = mean(I2rk4_2,1);
Rmodel_2 = mean(Rrk4_2,1);
Dmodel_2 = mean(Drk4_2,1);
%% Adam_Bashford_Predictor_Corrector=====================
% Sabpc_2 = Szero*ones(1,N+1);
% I1abpc_2 = I1zero*ones(1,N+1);
% I2abpc_2 = I2zero*ones(1,N+1);
% I3abpc_2 = I3zero*ones(1,N+1);
% Rabpc_2 = Rzero*ones(1,N+1);
% Dabpc_2 = Dzero*ones(1,N+1);
% 
% pcS1 = 0; pcS2 = 0; pcS3 = 0;
% pcI11 = 0; pcI12 = 0; pcI13 = 0;
% pcI21 = 0; pcI22 = 0; pcI23 = 0;
% pcI31 = 0; pcI32 = 0; pcI33 = 0;
% pcR1 = 0; pcR2 = 0; pcR3 = 0;
% pcD1 = 0; pcD2 = 0; pcD3 = 0;
% 
% for j = 1:L
%     Beta1 = Pastparamestimate(j,2);
%     Beta2 = Pastparamestimate(j,3);
%     Beta3 = Pastparamestimate(j,4);
%     r1 = Pastparamestimate(j,5);
%     r2 = Pastparamestimate(j,6);
%     r3 = Pastparamestimate(j,7);
%     d1 = Pastparamestimate(j,8);
%     d2 = Pastparamestimate(j,9);
%     d3 = Pastparamestimate(j,10);
%     p1 = Pastparamestimate(j,11);
%     p2 = Pastparamestimate(j,12);
%     p3 = Pastparamestimate(j,13);
%     gamma = Pastparamestimate(j,14);
% 
%     S1 = Sabpc_2(j);
%     I11 = I1abpc_2(j);
%     I21 = I2abpc_2(j);
%     I31 = I3abpc_2(j);
%     R1 = Rabpc_2(j);
%     D1 = Dabpc_2(j);
% 
%     f1_S = fS(S1,I11,I21,I31,R1,D1,Beta1,Beta2,Beta3,r1,r2,r3,d1,d2,d3,p1,p2,p3,gamma,sigma,mu);
%     f1_I1 = fI1(S1,I11,I21,I31,R1,D1,Beta1,Beta2,Beta3,r1,r2,r3,d1,d2,d3,p1,p2,p3,gamma,sigma,mu);
%     f1_I2 = fI2(S1,I11,I21,I31,R1,D1,Beta1,Beta2,Beta3,r1,r2,r3,d1,d2,d3,p1,p2,p3,gamma,sigma,mu);
%     f1_I3 = fI3(S1,I11,I21,I31,R1,D1,Beta1,Beta2,Beta3,r1,r2,r3,d1,d2,d3,p1,p2,p3,gamma,sigma,mu);
%     f1_R = fR(S1,I11,I21,I31,R1,D1,Beta1,Beta2,Beta3,r1,r2,r3,d1,d2,d3,p1,p2,p3,gamma,sigma,mu);
%     f1_D = fD(S1,I11,I21,I31,R1,D1,Beta1,Beta2,Beta3,r1,r2,r3,d1,d2,d3,p1,p2,p3,gamma,sigma,mu);
% 
%     pcS0 = f1_S;
%     pcI10 = f1_I1;
%     pcI20 = f1_I2;
%     pcI30 = f1_I3;
%     pcR0 = f1_R;
%     pcD0 = f1_D;
% 
%     S2 = Sabpc_2(j) + dt*(0.5*f1_S);
%     I12 = I1abpc_2(j) + dt*(0.5*f1_I1);
%     I22 = I2abpc_2(j) + dt*(0.5*f1_I2);
%     I32 = I3abpc_2(j) + dt*(0.5*f1_I3);
%     R2 = Rabpc_2(j) + dt*(0.5*f1_R);
%     D2 = Dabpc_2(j) + dt*(0.5*f1_D);
% 
%     f2_S = fS(S2,I12,I22,I32,R2,D2,Beta1,Beta2,Beta3,r1,r2,r3,d1,d2,d3,p1,p2,p3,gamma,sigma,mu);
%     f2_I1 = fI1(S2,I12,I22,I32,R2,D2,Beta1,Beta2,Beta3,r1,r2,r3,d1,d2,d3,p1,p2,p3,gamma,sigma,mu);
%     f2_I2 = fI2(S2,I12,I22,I32,R2,D2,Beta1,Beta2,Beta3,r1,r2,r3,d1,d2,d3,p1,p2,p3,gamma,sigma,mu);
%     f2_I3 = fI3(S2,I12,I22,I32,R2,D2,Beta1,Beta2,Beta3,r1,r2,r3,d1,d2,d3,p1,p2,p3,gamma,sigma,mu);
%     f2_R = fR(S2,I12,I22,I32,R2,D2,Beta1,Beta2,Beta3,r1,r2,r3,d1,d2,d3,p1,p2,p3,gamma,sigma,mu);
%     f2_D = fD(S2,I12,I22,I32,R2,D2,Beta1,Beta2,Beta3,r1,r2,r3,d1,d2,d3,p1,p2,p3,gamma,sigma,mu);
% 
%     S3 = Sabpc_2(j) + dt*(0.5*f2_S);
%     I13 = I1abpc_2(j) + dt*(0.5*f2_I1);
%     I23 = I2abpc_2(j) + dt*(0.5*f2_I2);
%     I33 = I3abpc_2(j) + dt*(0.5*f2_I3);
%     R3 = Rabpc_2(j) + dt*(0.5*f2_R);
%     D3 = Dabpc_2(j) + dt*(0.5*f2_D);
% 
%     f3_S = fS(S3,I13,I23,I33,R3,D3,Beta1,Beta2,Beta3,r1,r2,r3,d1,d2,d3,p1,p2,p3,gamma,sigma,mu);
%     f3_I1 = fI1(S3,I13,I23,I33,R3,D3,Beta1,Beta2,Beta3,r1,r2,r3,d1,d2,d3,p1,p2,p3,gamma,sigma,mu);
%     f3_I2 = fI2(S3,I13,I23,I33,R3,D3,Beta1,Beta2,Beta3,r1,r2,r3,d1,d2,d3,p1,p2,p3,gamma,sigma,mu);
%     f3_I3 = fI3(S3,I13,I23,I33,R3,D3,Beta1,Beta2,Beta3,r1,r2,r3,d1,d2,d3,p1,p2,p3,gamma,sigma,mu);
%     f3_R = fR(S3,I13,I23,I33,R3,D3,Beta1,Beta2,Beta3,r1,r2,r3,d1,d2,d3,p1,p2,p3,gamma,sigma,mu);
%     f3_D = fD(S3,I13,I23,I33,R3,D3,Beta1,Beta2,Beta3,r1,r2,r3,d1,d2,d3,p1,p2,p3,gamma,sigma,mu);
% 
%     S4 = Sabpc_2(j) + dt*f3_S;
%     I14 = I1abpc_2(j) + dt*f3_I1;
%     I24 = I2abpc_2(j) + dt*f3_I2;
%     I34 = I3abpc_2(j) + dt*f3_I3;
%     R4 = Rabpc_2(j) + dt*f3_R;
%     D4 = Dabpc_2(j) + dt*f3_D;
% 
%     f4_S = fS(S4,I14,I24,I34,R4,D4,Beta1,Beta2,Beta3,r1,r2,r3,d1,d2,d3,p1,p2,p3,gamma,sigma,mu);
%     f4_I1 = fI1(S4,I14,I24,I34,R4,D4,Beta1,Beta2,Beta3,r1,r2,r3,d1,d2,d3,p1,p2,p3,gamma,sigma,mu);
%     f4_I2 = fI2(S4,I14,I24,I34,R4,D4,Beta1,Beta2,Beta3,r1,r2,r3,d1,d2,d3,p1,p2,p3,gamma,sigma,mu);
%     f4_I3 = fI3(S4,I14,I24,I34,R4,D4,Beta1,Beta2,Beta3,r1,r2,r3,d1,d2,d3,p1,p2,p3,gamma,sigma,mu);
%     f4_R = fR(S4,I14,I24,I34,R4,D4,Beta1,Beta2,Beta3,r1,r2,r3,d1,d2,d3,p1,p2,p3,gamma,sigma,mu);
%     f4_D = fD(S4,I14,I24,I34,R4,D4,Beta1,Beta2,Beta3,r1,r2,r3,d1,d2,d3,p1,p2,p3,gamma,sigma,mu);
% 
%     Sabpc_2(j+1) = Sabpc_2(j) + dt*((1/6)*f1_S + (1/3)*f2_S + (1/3)*f3_S + (1/6)*f4_S);
%     I1abpc_2(j+1) = I1abpc_2(j) + dt*((1/6)*f1_I1 + (1/3)*f2_I1 + (1/3)*f3_I1 + (1/6)*f4_I1);
%     I2abpc_2(j+1) = I2abpc_2(j) + dt*((1/6)*f1_I2 + (1/3)*f2_I2 + (1/3)*f3_I2 + (1/6)*f4_I2);
%     I3abpc_2(j+1) = I3abpc_2(j) + dt*((1/6)*f1_I3 + (1/3)*f2_I3 + (1/3)*f3_I3 + (1/6)*f4_I3);
%     Rabpc_2(j+1) = Rabpc_2(j) + dt*((1/6)*f1_R + (1/3)*f2_R + (1/3)*f3_R + (1/6)*f4_R);
%     Dabpc_2(j+1) = Dabpc_2(j) + dt*((1/6)*f1_D + (1/3)*f2_D + (1/3)*f3_D + (1/6)*f4_D);
% 
%     pcS1 = pcS0; pcS2 = pcS1; pcS3 = pcS2;
%     pcI11 = pcI10; pcI12 = pcI11; pcI13 = pcI12;
%     pcI21 = pcI20; pcI22 = pcI21; pcI23 = pcI22;
%     pcI31 = pcI30; pcI32 = pcI31; pcI33 = pcI32;
%     pcR1 = pcR0; pcR2 = pcR1; pcR3 = pcR2;
%     pcD1 = pcD0; pcD2 = pcD1; pcD3 = pcD2;
% end
%         
% for j = 1:L
%     Beta1 = Pastparamestimate(j,2);
%     Beta2 = Pastparamestimate(j,3);
%     Beta3 = Pastparamestimate(j,4);
%     r1 = Pastparamestimate(j,5);
%     r2 = Pastparamestimate(j,6);
%     r3 = Pastparamestimate(j,7);
%     d1 = Pastparamestimate(j,8);
%     d2 = Pastparamestimate(j,9);
%     d3 = Pastparamestimate(j,10);
%     p1 = Pastparamestimate(j,11);
%     p2 = Pastparamestimate(j,12);
%     p3 = Pastparamestimate(j,13);
%     gamma = Pastparamestimate(j,14);
%     
%     S1 = Sabpc_2(j);
%     I11 = I1abpc_2(j);
%     I21 = I2abpc_2(j);
%     I31 = I3abpc_2(j);
%     R1 = Rabpc_2(j);
%     D1 = Dabpc_2(j);
% 
%     f1_S = fS(S1,I11,I21,I31,R1,D1,Beta1,Beta2,Beta3,r1,r2,r3,d1,d2,d3,p1,p2,p3,gamma,sigma,mu);
%     f1_I1 = fI1(S1,I11,I21,I31,R1,D1,Beta1,Beta2,Beta3,r1,r2,r3,d1,d2,d3,p1,p2,p3,gamma,sigma,mu);
%     f1_I2 = fI2(S1,I11,I21,I31,R1,D1,Beta1,Beta2,Beta3,r1,r2,r3,d1,d2,d3,p1,p2,p3,gamma,sigma,mu);
%     f1_I3 = fI3(S1,I11,I21,I31,R1,D1,Beta1,Beta2,Beta3,r1,r2,r3,d1,d2,d3,p1,p2,p3,gamma,sigma,mu);
%     f1_R = fR(S1,I11,I21,I31,R1,D1,Beta1,Beta2,Beta3,r1,r2,r3,d1,d2,d3,p1,p2,p3,gamma,sigma,mu);
%     f1_D = fD(S1,I11,I21,I31,R1,D1,Beta1,Beta2,Beta3,r1,r2,r3,d1,d2,d3,p1,p2,p3,gamma,sigma,mu);
% 
%     pcS0 = f1_S;
%     pcI10 = f1_I1;
%     pcI20 = f1_I2;
%     pcI30 = f1_I3;
%     pcR0 = f1_R;
%     pcD0 = f1_D;
% 
%     wS = Sabpc_2(j) + dt*( 55.0 * pcS0 - 59.0 * pcS1 + 37.0 * pcS2 - 9.0 * pcS3 ) / 24.0;
%     wI1 = I1abpc_2(j) + dt*( 55.0 * pcI10 - 59.0 * pcI11 + 37.0 * pcI12 - 9.0 * pcI13 ) / 24.0;
%     wI2 = I2abpc_2(j) + dt*( 55.0 * pcI20 - 59.0 * pcI21 + 37.0 * pcI22 - 9.0 * pcI23 ) / 24.0;
%     wI3 = I3abpc_2(j) + dt*( 55.0 * pcI30 - 59.0 * pcI31 + 37.0 * pcI32 - 9.0 * pcI33 ) / 24.0;
%     wR = Rabpc_2(j) + dt*( 55.0 * pcR0 - 59.0 * pcR1 + 37.0 * pcR2 - 9.0 * pcR3 ) / 24.0;
%     wD = Dabpc_2(j) + dt*( 55.0 * pcD0 - 59.0 * pcD1 + 37.0 * pcD2 - 9.0 * pcD3 ) / 24.0;
% 
%     fw_S = fS(wS,wI1,wI2,wI3,wR,wD,Beta1,Beta2,Beta3,r1,r2,r3,d1,d2,d3,p1,p2,p3,gamma,sigma,mu);
%     fw_I1 = fI1(wS,wI1,wI2,wI3,wR,wD,Beta1,Beta2,Beta3,r1,r2,r3,d1,d2,d3,p1,p2,p3,gamma,sigma,mu);
%     fw_I2 = fI2(wS,wI1,wI2,wI3,wR,wD,Beta1,Beta2,Beta3,r1,r2,r3,d1,d2,d3,p1,p2,p3,gamma,sigma,mu);
%     fw_I3 = fI3(wS,wI1,wI2,wI3,wR,wD,Beta1,Beta2,Beta3,r1,r2,r3,d1,d2,d3,p1,p2,p3,gamma,sigma,mu);
%     fw_R = fR(wS,wI1,wI2,wI3,wR,wD,Beta1,Beta2,Beta3,r1,r2,r3,d1,d2,d3,p1,p2,p3,gamma,sigma,mu);
%     fw_D = fD(wS,wI1,wI2,wI3,wR,wD,Beta1,Beta2,Beta3,r1,r2,r3,d1,d2,d3,p1,p2,p3,gamma,sigma,mu);
%     
%     Sabpc_2(j+1) = Sabpc_2(j) + dt*( 9.0 * fw_S + 19.0 * pcS0 - 5.0 * pcS1 + pcS2 ) / 24.0;
%     I1abpc_2(j+1) = I1abpc_2(j) + dt*( 9.0 * fw_I1 + 19.0 * pcI10 - 5.0 * pcI11 + pcI12 ) / 24.0;
%     I2abpc_2(j+1) = I2abpc_2(j) + dt*( 9.0 * fw_I2 + 19.0 * pcI20 - 5.0 * pcI21 + pcI22 ) / 24.0;
%     I3abpc_2(j+1) = I3abpc_2(j) + dt*( 9.0 * fw_I3 + 19.0 * pcI30 - 5.0 * pcI31 + pcI32 ) / 24.0;
%     Rabpc_2(j+1) = Rabpc_2(j) + dt*( 9.0 * fw_R + 19.0 * pcR0 - 5.0 * pcR1 + pcR2 ) / 24.0;
%     Dabpc_2(j+1) = Dabpc_2(j) + dt*( 9.0 * fw_D + 19.0 * pcD0 - 5.0 * pcD1 + pcD2 ) / 24.0;
% 
%     pcS1 = pcS0; pcS2 = pcS1; pcS3 = pcS2;
%     pcI11 = pcI10; pcI12 = pcI11; pcI13 = pcI12;
%     pcI21 = pcI20; pcI22 = pcI21; pcI23 = pcI22;
%     pcI31 = pcI30; pcI32 = pcI31; pcI33 = pcI32;
%     pcR1 = pcR0; pcR2 = pcR1; pcR3 = pcR2;
%     pcD1 = pcD0; pcD2 = pcD1; pcD3 = pcD2;
% end
% % RMSE_srk4x = rmse(Xsrk4_2,Datax);
% % RMSE_srk4y = rmse(Ysrk4_2,Datay);
% Save_abpc = mean(Sabpc_2,1);
% I1ave_abpc = mean(I1abpc_2,1);
% I2ave_abpc = mean(I2abpc_2,1);
% I3ave_abpc = mean(I3abpc_2,1);
% Rave_abpc = mean(Rabpc_2,1);
% Dave_abpc = mean(Dabpc_2,1);
% % Smodel = Save_srk4;
% SDet = Save_abpc(2:end);
% I1Det = I1ave_abpc(2:end);
% I2Det = I2ave_abpc(2:end);
% I3Det = I3ave_abpc(2:end);
% RDet = Rave_abpc(2:end);
% DDet = Dave_abpc(2:end);
%% Plot
%     Top plot
%     ax1 = nexttile;
%     plot(tdays,Sdata,"-s",'Color',"#000000",'MarkerIndices',1:20:length(tdays)); hold on;
%     plot(tdays,Smodel_2,':', 'Color',"#2F4F4F"); hold on;
%     plot(tdays,SDet,"--", 'Color',"#778899");
%     plot(tdays,Smodel(2,:),'-', 'Color',"#696969"); hold on;
%     cax1 = [105/255, 105/255, 105/255];
%     fill([tdays, fliplr(tdays)], [Smodel(1,:), fliplr(Smodel(3,:))], cax1, 'LineStyle', "None",'FaceAlpha',0.2); hold on;
%     title(ax1,'Susceptible Population, S')
%     ax1.FontSize = 8;
%     hold off;
%     legend('Susceptible Data','RK4','Adam-Moulton PC','SRK4','Location','best')
% %     legend('Susceptible Data','Simulation RK4')
%     xlabel('time (date)'); ylabel('Susceptible Population')    
    

%     ax4 = nexttile; 
%     plot(tdays,I3data,'-s','Color',"#FF4500",'MarkerIndices',1:20:length(tdays)); hold on;
% %     plot(tdays,I3model_2,':','Color',"#FF7F50"); hold on;
%     plot(tdays,I3Det,'--', 'Color',"#FFA500")
% %     plot(tdays,I3model(2,:),'-','Color',"#FFD700"); hold on;
%     title(ax4,'Ventilated Population, I3')
%     cax4 = [204/255,141/255,61/255];
% %     fill([tdays, fliplr(tdays)], [I3model(1,:), fliplr(I3model(3,:))], cax4, 'LineStyle', "None",'FaceAlpha',0.2); hold on;
%     ax4.FontSize = 8;
% %     legend('Infected Data','RK4','Adam-Moulton PC','SRK4','Location','best')
%     hold off;
%     xlabel('time (date)'); ylabel('Ventilated Population')    
    
%     ax5 = nexttile;
%     plot(tdays,Rdata,'-s','Color',"#008080",'MarkerIndices',1:20:length(tdays)); hold on;
%     plot(tdays,Rmodel_2,':','Color',"#556B2F");hold on;
%     plot(tdays,RDet,'--', 'Color',"#3CB371")
%     plot(tdays,Rmodel(2,:),'-','Color',"#00FF00");hold on;
%     title(ax5,'Recovered Population, R')
%     cax5 = [70/255,192/255,30/255];
%     fill([tdays, fliplr(tdays)], [Rmodel(1,:), fliplr(Rmodel(3,:))], cax5, 'LineStyle', "None",'FaceAlpha',0.2); hold on;
%     ax5.FontSize = 8;
%     legend('Recovered Data','RK4','Adam-Moulton PC','SRK4','Location','best')
%     hold off;
%     xlabel('time (date)'); ylabel('Recovered Population')    
%     
%     ax6 = nexttile;
%     plot(tdays,Ddata,'-s','Color',"#8B0000",'MarkerIndices',1:20:length(tdays)); hold on;
%     plot(tdays,Dmodel_2,':','Color',"#FF0000");hold on;
%     plot(tdays,DDet,'--', 'Color',"#8B4513")
%     plot(tdays,Dmodel(2,:),'-','Color',"#F08080");hold on;
%     title(ax6,'Death Population, D')
%     cax6 = [217/255,46/255,157/255];
%     fill([tdays, fliplr(tdays)], [Dmodel(1,:), fliplr(Dmodel(3,:))], cax6, 'LineStyle', "None",'FaceAlpha',0.2); hold on;
%     ax6.FontSize = 8;
%     legend('Death Data','RK4','Adam-Moulton PC','SRK4','Location','best')
%     hold off;
%     xlabel('time (date)'); ylabel('Death Population') 

% figure('Position', [1, 1, 1400, 900])
datafile1 = readtable('Malaysia-SI2RD-param-estimate.csv');
% datafile2 = readtable('Canada-SI3RD-param-estimate.csv');
% datafile3 = readtable('Italy-SI3RD-param-estimate.csv');
% datafile4 = readtable('Slovakia-SI3RD-param-estimate.csv');
% datafile5 = readtable('South-Africa-SI3RD-param-estimate.csv');
% datafile9 = readtable('Curvefitting_SIRD.csv');

SIRDtMalaysia = datafile1{1:end,1};
SIRDbeta1Malaysia = datafile1{1:end,2};
SIRDbeta2Malaysia = datafile1{1:end,3};
SIRDr1Malaysia = datafile1{1:end,4};
SIRDr2Malaysia = datafile1{1:end,5};
SIRDd1Malaysia = datafile1{1:end,6};
SIRDd2Malaysia = datafile1{1:end,7};
SIRDp1Malaysia = datafile1{1:end,8};
SIRDp2Malaysia = datafile1{1:end,9};

SIRDp1Malaysia(SIRDp1Malaysia < 0) = 0.00001;
SIRDp2Malaysia(SIRDp2Malaysia < 0.00001) = 0.00001;
SIRDp1Malaysia = SIRDp1Malaysia./(SIRDp1Malaysia + SIRDp2Malaysia);
SIRDp2Malaysia = SIRDp2Malaysia./(SIRDp1Malaysia + SIRDp2Malaysia);

SIRDbeta1Malaysia(SIRDbeta1Malaysia < 0) = 0.00001;
SIRDbeta2Malaysia(SIRDbeta2Malaysia < 0) = 0.00001;
SIRDr1Malaysia(SIRDr1Malaysia < 0) = 0.00001;
SIRDr2Malaysia(SIRDr2Malaysia < 0) = 0.00001;
SIRDd1Malaysia(SIRDd1Malaysia < 0) = 0.00001;
SIRDd2Malaysia(SIRDd2Malaysia < 0) = 0.00001;
 
% SIRDtCanada = datafile2{1:end,1};
% SIRDbeta1Canada = datafile2{1:end,2};
% SIRDbeta2Canada = datafile2{1:end,3};
% SIRDbeta3Canada = datafile2{1:end,4};
% SIRDr1Canada = datafile2{1:end,5};
% SIRDr2Canada = datafile2{1:end,6};
% SIRDr3Canada = datafile2{1:end,7};
% SIRDd1Canada = datafile2{1:end,8};
% SIRDd2Canada = datafile2{1:end,9};
% SIRDd3Canada = datafile2{1:end,10};
% SIRDp1Canada = datafile2{1:end,11};
% SIRDp2Canada = datafile2{1:end,12};
% SIRDp3Canada = datafile2{1:end,13};
%  
% SIRDtItaly = datafile3{1:end,1};
% SIRDbeta1Italy = datafile3{1:end,2};
% SIRDbeta2Italy = datafile3{1:end,3};
% SIRDbeta3Italy = datafile3{1:end,4};
% SIRDr1Italy = datafile3{1:end,5};
% SIRDr2Italy = datafile3{1:end,6};
% SIRDr3Italy = datafile3{1:end,7};
% SIRDd1Italy = datafile3{1:end,8};
% SIRDd2Italy = datafile3{1:end,9};
% SIRDd3Italy = datafile3{1:end,10};
% SIRDp1Italy = datafile3{1:end,11};
% SIRDp2Italy = datafile3{1:end,12};
% SIRDp3Italy = datafile3{1:end,13};
 
% SIRDtSlovakia = datafile4{1:end,1};
% SIRDbeta1Slovakia = datafile4{1:end,2};
% SIRDbeta2Slovakia = datafile4{1:end,3};
% SIRDbeta3Slovakia = datafile4{1:end,4};
% SIRDr1Slovakia = datafile4{1:end,5};
% SIRDr2Slovakia = datafile4{1:end,6};
% SIRDr3Slovakia = datafile4{1:end,7};
% SIRDd1Slovakia = datafile4{1:end,8};
% SIRDd2Slovakia = datafile4{1:end,9};
% SIRDd3Slovakia = datafile4{1:end,10};
% SIRDp1Slovakia = datafile4{1:end,11};
% SIRDp2Slovakia = datafile4{1:end,12};
% SIRDp3Slovakia = datafile4{1:end,13};
 
% SIRDtSouthAfrica = datafile5{1:end,1};
% SIRDbeta1SouthAfrica = datafile5{1:end,2};
% SIRDbeta2SouthAfrica = datafile5{1:end,3};
% SIRDbeta3SouthAfrica = datafile5{1:end,4};
% SIRDr1SouthAfrica = datafile5{1:end,5};
% SIRDr2SouthAfrica = datafile5{1:end,6};
% SIRDr3SouthAfrica = datafile5{1:end,7};
% SIRDd1SouthAfrica = datafile5{1:end,8};
% SIRDd2SouthAfrica = datafile5{1:end,9};
% SIRDd3SouthAfrica = datafile5{1:end,10};
% SIRDp1SouthAfrica = datafile5{1:end,11};
% SIRDp2SouthAfrica = datafile5{1:end,12};
% SIRDp3SouthAfrica = datafile5{1:end,13};

% Curvefit_Malaysia_r1 = datafile9{:,25}(~all(isnan(datafile9{:,25}), 2));
% Curvefit_Malaysia_r2 = datafile9{:,26}(~all(isnan(datafile9{:,26}), 2));
% Curvefit_Malaysia_r3 = datafile9{:,27}(~all(isnan(datafile9{:,27}), 2));
% Curvefit_Malaysia_d1 = datafile9{:,28}(~all(isnan(datafile9{:,28}), 2));
% Curvefit_Malaysia_d2 = datafile9{:,29}(~all(isnan(datafile9{:,29}), 2));
% Curvefit_Malaysia_d3 = datafile9{:,30}(~all(isnan(datafile9{:,30}), 2));

% Curvefit_Canada_r1 = datafile9{:,31}(~all(isnan(datafile9{:,31}), 2));
% Curvefit_Canada_r2 = datafile9{:,32}(~all(isnan(datafile9{:,32}), 2));
% Curvefit_Canada_r3 = datafile9{:,33}(~all(isnan(datafile9{:,33}), 2));
% Curvefit_Canada_d1 = datafile9{:,34}(~all(isnan(datafile9{:,34}), 2));
% Curvefit_Canada_d2 = datafile9{:,35}(~all(isnan(datafile9{:,35}), 2));
% Curvefit_Canada_d3 = datafile9{:,36}(~all(isnan(datafile9{:,36}), 2));
% 
% Curvefit_Italy_r1 = datafile9{:,37}(~all(isnan(datafile9{:,37}), 2));
% Curvefit_Italy_r2 = datafile9{:,38}(~all(isnan(datafile9{:,38}), 2));
% Curvefit_Italy_r3 = datafile9{:,39}(~all(isnan(datafile9{:,39}), 2));
% Curvefit_Italy_d1 = datafile9{:,40}(~all(isnan(datafile9{:,40}), 2));
% Curvefit_Italy_d2 = datafile9{:,41}(~all(isnan(datafile9{:,41}), 2));
% Curvefit_Italy_d3 = datafile9{:,42}(~all(isnan(datafile9{:,42}), 2));

% Curvefit_Slovakia_r1 = datafile9{:,43}(~all(isnan(datafile9{:,43}), 2));
% Curvefit_Slovakia_r2 = datafile9{:,44}(~all(isnan(datafile9{:,44}), 2));
% Curvefit_Slovakia_r3 = datafile9{:,45}(~all(isnan(datafile9{:,45}), 2));
% Curvefit_Slovakia_d1 = datafile9{:,46}(~all(isnan(datafile9{:,46}), 2));
% Curvefit_Slovakia_d2 = datafile9{:,47}(~all(isnan(datafile9{:,47}), 2));
% Curvefit_Slovakia_d3 = datafile9{:,48}(~all(isnan(datafile9{:,48}), 2));

% Curvefit_SouthAfrica_r1 = datafile9{:,49}(~all(isnan(datafile9{:,49}), 2));
% Curvefit_SouthAfrica_r2 = datafile9{:,50}(~all(isnan(datafile9{:,50}), 2));
% Curvefit_SouthAfrica_r3 = datafile9{:,51}(~all(isnan(datafile9{:,51}), 2));
% Curvefit_SouthAfrica_d1 = datafile9{:,52}(~all(isnan(datafile9{:,52}), 2));
% Curvefit_SouthAfrica_d2 = datafile9{:,53}(~all(isnan(datafile9{:,53}), 2));
% Curvefit_SouthAfrica_d3 = datafile9{:,54}(~all(isnan(datafile9{:,54}), 2));

% Curvefit_Malaysia_p1 = datafile9{:,55}(~all(isnan(datafile9{:,55}), 2));
% Curvefit_Malaysia_p2 = datafile9{:,56}(~all(isnan(datafile9{:,56}), 2));
% Curvefit_Malaysia_p3 = datafile9{:,57}(~all(isnan(datafile9{:,57}), 2));
% Curvefit_Canada_p1 = datafile9{:,58}(~all(isnan(datafile9{:,58}), 2));
% Curvefit_Canada_p2 = datafile9{:,59}(~all(isnan(datafile9{:,59}), 2));
% Curvefit_Canada_p3 = datafile9{:,60}(~all(isnan(datafile9{:,60}), 2));
% Curvefit_Italy_p1 = datafile9{:,61}(~all(isnan(datafile9{:,61}), 2));
% Curvefit_Italy_p2 = datafile9{:,62}(~all(isnan(datafile9{:,62}), 2));
% Curvefit_Italy_p3 = datafile9{:,63}(~all(isnan(datafile9{:,63}), 2));
% Curvefit_Slovakia_p1 = datafile9{:,64}(~all(isnan(datafile9{:,64}), 2));
% Curvefit_Slovakia_p2 = datafile9{:,65}(~all(isnan(datafile9{:,65}), 2));
% Curvefit_Slovakia_p3 = datafile9{:,66}(~all(isnan(datafile9{:,66}), 2));
% Curvefit_SouthAfrica_p1 = datafile9{:,67}(~all(isnan(datafile9{:,67}), 2));
% Curvefit_SouthAfrica_p2 = datafile9{:,68}(~all(isnan(datafile9{:,68}), 2));
% Curvefit_SouthAfrica_p3 = datafile9{:,69}(~all(isnan(datafile9{:,69}), 2));

figure('Units','centimeters','Position',[2 2 20 18]); % 20cm width x 18cm height (~2/3 A4)

% --- Layout parameters ---
rowsSmall = 4; colsSmall = 2;         
gapX = 0.07; gapY = 0.06;             
marginLeft = 0.08; marginRight = 0.05;
marginTop = 0.05; marginBottom = 0.05;

widthSmall = (1 - marginLeft - marginRight - (colsSmall-1)*gapX)/colsSmall;
heightSmall = (0.45 - marginTop - (rowsSmall-1)*gapY)/rowsSmall + 0.03;  % shrink top half

% --- Small plots ---
smallVars = { ...
    SIRDbeta1Malaysia, SIRDbeta2Malaysia; ...
    SIRDr1Malaysia, SIRDr2Malaysia; ...
    SIRDd1Malaysia, SIRDd2Malaysia; ...
    SIRDp1Malaysia, SIRDp2Malaysia ...
};

smallTitles = { ...
    'Infectious rate, b_1','Infectious rate, b_2'; ...
    'Recovery rate, r_1','Recovery rate, r_2'; ...
    'Fatality rate, d_1','Fatality rate, d_2'; ...
    'Proportion of Population, p_1','Proportion of Population, p_2' ...
};

xlabelLetters = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)'};
count = 1;

for row = 1:rowsSmall
    for col = 1:colsSmall
        left = marginLeft + (col-1)*(widthSmall + gapX);
        bottom = 1 - marginTop - row*heightSmall - (row-1)*gapY;
        ax = axes('Position',[left bottom widthSmall heightSmall]);
        plot(tdays, smallVars{row,col}, '-', 'Color', "#191970"); hold on;
        title(smallTitles{row,col},'FontSize',8);
        ylabel([smallTitles{row,col}(1:end-5) ' Trend'],'FontSize',7); 
        xlabel(['                                        ' xlabelLetters{count} '                            date (days)'],'FontSize',7);
        ax.FontSize = 10;
        hold off;
        count = count + 1;
    end
end

% --- R_t plot ---
left = marginLeft; width = 1 - marginLeft - marginRight;
bottom = 0.37 - heightSmall;  % slightly lower
height = 0.08;
axR = axes('Position',[left bottom width height]);
r0 = SIRDp1Malaysia.*(SIRDbeta1Malaysia./(SIRDr1Malaysia+SIRDd1Malaysia)) + ...
     SIRDp2Malaysia.*(SIRDbeta2Malaysia./(SIRDr2Malaysia+SIRDd2Malaysia));
plot(tdays,r0,'-','Color',"#191970"); hold on;
plot(tdays, ones(size(r0)),'--','Color',"#191970");
title('Real-time Reproductive Number, R_t','FontSize',8);
xlabel('                                                                                                           (i)                                                                                           date (days)','FontSize',7);
ylabel('R_t','FontSize',7);
axR.FontSize = 10;
hold off;

% --- Large plots (bottom) ---
largeHeight = 0.19; gapLarge = 0.05;
widthLarge = (1 - marginLeft - marginRight - gapLarge)/2;
bottomLarge = marginBottom;

% Large plot 1: I1
axI1 = axes('Position',[marginLeft bottomLarge widthLarge largeHeight]);
plot(tdays, I1data, '-s','Color',"#2C00D2",'MarkerIndices',1:20:length(tdays)); hold on;
plot(tdays, I1model_2, ':','Color',"#273DE3");
title('Infected Population, I1','FontSize',8);
xlabel('                                                    (j)                                  date (days)','FontSize',7); 
ylabel('Population','FontSize',7);
legend('Infected Data','Infected Simulation','Location','best'); axI1.FontSize = 10;
hold off;

% Large plot 2: I2
axI2 = axes('Position',[marginLeft + widthLarge + gapLarge bottomLarge widthLarge largeHeight]);
plot(tdays, I2data, '-s','Color',"#4B0082",'MarkerIndices',1:20:length(tdays)); hold on;
plot(tdays, I2model_2, ':','Color',"#7B68EE");
title('Hospitalised Population, I2','FontSize',10);
xlabel('                                                    (k)                                  date (days)','FontSize',7); 
ylabel('Population','FontSize',7);
legend('Hospitalized Data','Hospitalized Simulation','Location','best'); 
axI2.FontSize=10; hold off;

%% Save
print('-depsc','SI3RDparamanalysis.eps');

%% RMSE

%% RMSE
% RMSE = zeros(3,6);
% RMSE(3,1) = sqrt(sum((Save_srk4(:) - Sdata(:)).^2)/length(Save_srk4));
% RMSE(1,1) = sqrt(sum((Smodel_2(:) - Sdata(:)).^2)/length(Save_srk4));
% RMSE(2,1) = sqrt(sum((SDet(:) - Sdata(:)).^2)/length(Save_srk4));
% RMSE(3,2) = sqrt(sum((I1ave_srk4(:) - I1data(:)).^2)/length(Save_srk4));
% RMSE(1,2) = sqrt(sum((I1model_2(:) - I1data(:)).^2)/length(Save_srk4));
% RMSE(2,2) = sqrt(sum((I1Det(:) - I1data(:)).^2)/length(Save_srk4));
% RMSE(3,3) = sqrt(sum((I2ave_srk4(:) - I2data(:)).^2)/length(Save_srk4));
% RMSE(1,3) = sqrt(sum((I2model_2(:) - I2data(:)).^2)/length(Save_srk4));
% RMSE(2,3) = sqrt(sum((I2Det(:) - I2data(:)).^2)/length(Save_srk4));
% RMSE(3,4) = sqrt(sum((I3ave_srk4(:) - I3data(:)).^2)/length(Save_srk4));
% RMSE(1,4) = sqrt(sum((I3model_2(:) - I3data(:)).^2)/length(Save_srk4));
% RMSE(2,4) = sqrt(sum((I3Det(:) - I3data(:)).^2)/length(Save_srk4));
% RMSE(3,5) = sqrt(sum((Rave_srk4(:) - Rdata(:)).^2)/length(Save_srk4));
% RMSE(1,5) = sqrt(sum((Rmodel_2(:) - Rdata(:)).^2)/length(Save_srk4));
% RMSE(2,5) = sqrt(sum((RDet(:) - Rdata(:)).^2)/length(Save_srk4));
% RMSE(3,6) = sqrt(sum((Dave_srk4(:) - Ddata(:)).^2)/length(Save_srk4));
% RMSE(1,6) = sqrt(sum((Dmodel_2(:) - Ddata(:)).^2)/length(Save_srk4));
% RMSE(2,6) = sqrt(sum((DDet(:) - Ddata(:)).^2)/length(Save_srk4));
% 
% Outputdata = cat(1, RDet, I1Det, DDet, I2Det, I3Det, SDet, Rmodel_2, I1model_2, Dmodel_2, I2model_2, I3model_2, Smodel_2,Rave_srk4, I1ave_srk4, Dave_srk4, I2ave_srk4, I3ave_srk4, Save_srk4);


