clear model data params options
datafile = readmatrix('si4rd_data_Malaysia.csv');
Ptotal = 3.153*10^7;%malaysia
rng = 3683;

I1data = datafile(13:end,2);
I2data = datafile(13:end,5);
I3data = datafile(13:end,6);
I4data = datafile(13:end,7);
Rdata = datafile(13:end,3);
Ddata = datafile(13:end,4);
Sdata = datafile(13:end,8);
tdata = datafile(13:end,1);
Hdata = datafile(13:end,9);
p_vax = 1;

t = 1:1:length(I1data);
tdays = datetime(2020,3,24) + caldays(1:length(t));
Pastparamestimate = readmatrix('Malaysia-SI4RD-param-estimate.csv');
Szero = Sdata(1);
I1zero = I1data(1);
I2zero = I2data(1);
I3zero = I3data(1);
I4zero = I4data(1);
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
fS  = @(S,I1,I2,I3,I4,R,D,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu) mu*(S+I1+I2+I3+14+R)-S./Ptotal.*(Beta1.*I1+Beta2.*I2+Beta3.*I3+Beta4.*I4)+gamma.*R - mu.*S;
gS  = @(S,I1,I2,I3,I4,R,D,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu) -sigma.*I1.*S./Ptotal;
fI1 = @(S,I1,I2,I3,I4,R,D,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu) p1.*S./Ptotal.*(Beta1.*I1+Beta2.*I2+Beta3.*I3+Beta4.*I4)-(r1+d1).*I1 - mu.*I1;
gI1 = @(S,I1,I2,I3,I4,R,D,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu) p1.*sigma.*I1.*S./Ptotal;
fI2 = @(S,I1,I2,I3,I4,R,D,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu) p2.*S./Ptotal.*(Beta1.*I1+Beta2.*I2+Beta3.*I3+Beta4.*I4)-(r2+d2).*I2 - mu.*I2;
gI2 = @(S,I1,I2,I3,I4,R,D,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu) p2.*sigma.*I1.*S./Ptotal;
fI3 = @(S,I1,I2,I3,I4,R,D,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu) p3.*S./Ptotal.*(Beta1.*I1+Beta2.*I2+Beta3.*I3+Beta4.*I4)-(r3+d3).*I3 - mu.*I3;
gI3 = @(S,I1,I2,I3,I4,R,D,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu) p3.*sigma.*I1.*S./Ptotal;
fI4 = @(S,I1,I2,I3,I4,R,D,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu) p4.*S./Ptotal.*(Beta1.*I1+Beta2.*I2+Beta3.*I3+Beta4.*I4)-(r4+d4).*I4 - mu.*I4;
gI4 = @(S,I1,I2,I3,I4,R,D,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu) p4.*sigma.*I1.*S./Ptotal;
fR  = @(S,I1,I2,I3,I4,R,D,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu) r1.*I1 + r2.*I2 + r3.*I3 + r4.*I4 - gamma.*R - mu.*R;
gR  = @(S,I1,I2,I3,I4,R,D,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu) 0;
fD  = @(S,I1,I2,I3,I4,R,D,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu) d1.*I1 + d2.*I2 + d3.*I3 + d4.*I4;
gD  = @(S,I1,I2,I3,I4,R,D,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu) 0;
%% SRK4==================================================
sigma = 0.05;
%======================================================%SRK-4 Method
Ssrk4_2 = zeros(M,N);
I1srk4_2 = zeros(M,N);
I2srk4_2 = zeros(M,N);
I3srk4_2 = zeros(M,N);
I4srk4_2 = zeros(M,N);
Rsrk4_2 = zeros(M,N);
Dsrk4_2 = zeros(M,N);
Stemp = Szero*ones(M,1);
I1temp = I1zero*ones(M,1);
I2temp = I2zero*ones(M,1);
I3temp = I3zero*ones(M,1);
I4temp = I4zero*ones(M,1);
Rtemp = Rzero*ones(M,1);
Dtemp = Dzero*ones(M,1);
Smodel = ones(3,N);
I1model = ones(3,N);
I2model = ones(3,N);
I3model = ones(3,N);
I4model = ones(3,N);
Rmodel = ones(3,N);
Dmodel = ones(3,N);
SdW1 = sqrt(dt)*randn(M,N);       %Numerical Brownian increments
SdW2 = sqrt(dt)*randn(M,N);       %Numerical Brownian increments
I1dW1 = sqrt(dt)*randn(M,N);       %Numerical Brownian increments
I1dW2 = sqrt(dt)*randn(M,N);       %Numerical Brownian increments
I2dW1 = sqrt(dt)*randn(M,N);       %Numerical Brownian increments
I2dW2 = sqrt(dt)*randn(M,N);       %Numerical Brownian increments
I3dW1 = sqrt(dt)*randn(M,N);       %Numerical Brownian increments
I3dW2 = sqrt(dt)*randn(M,N);       %Numerical Brownian increments
I4dW1 = sqrt(dt)*randn(M,N);       %Numerical Brownian increments
I4dW2 = sqrt(dt)*randn(M,N);       %Numerical Brownian increments
RdW1 = sqrt(dt)*randn(M,N);       %Numerical Brownian increments
RdW2 = sqrt(dt)*randn(M,N);       %Numerical Brownian increments
DdW1 = sqrt(dt)*randn(M,N);       %Numerical Brownian increments
DdW2 = sqrt(dt)*randn(M,N);       %Numerical Brownian increments

for j = 1:L
    p1 = Pastparamestimate(j,2);
    p2 = Pastparamestimate(j,3);
    p3 = Pastparamestimate(j,4);
    p4 = Pastparamestimate(j,5);
    Beta1 = Pastparamestimate(j,6);
    Beta2 = Pastparamestimate(j,7);
    Beta3 = Pastparamestimate(j,8);
    Beta4 = Pastparamestimate(j,9);
    r1 = Pastparamestimate(j,10);
    r2 = Pastparamestimate(j,11);
    r3 = Pastparamestimate(j,12);
    r4 = Pastparamestimate(j,13);
    d1 = Pastparamestimate(j,14);
    d2 = Pastparamestimate(j,15);
    d3 = Pastparamestimate(j,16);
    d4 = Pastparamestimate(j,17);
    gamma = Pastparamestimate(j,18);
    
    ptotal = p1 + p2 + p3 + p4;
    p1 = p1/ptotal*p_vax;
    p2 = p2/ptotal*p_vax;
    p3 = p3/ptotal*p_vax;
    p4 = p4/ptotal*p_vax;

    SWinc1 = sum(SdW1(:,R*(j-1)+1:R*j),2);
    SWinc2 = sum(SdW2(:,R*(j-1)+1:R*j),2);
    I1Winc1 = sum(I1dW1(:,R*(j-1)+1:R*j),2);
    I1Winc2 = sum(I1dW2(:,R*(j-1)+1:R*j),2);
    I2Winc1 = sum(I2dW1(:,R*(j-1)+1:R*j),2);
    I2Winc2 = sum(I2dW2(:,R*(j-1)+1:R*j),2);
    I3Winc1 = sum(I3dW1(:,R*(j-1)+1:R*j),2);
    I3Winc2 = sum(I3dW2(:,R*(j-1)+1:R*j),2);
    I4Winc1 = sum(I4dW1(:,R*(j-1)+1:R*j),2);
    I4Winc2 = sum(I4dW2(:,R*(j-1)+1:R*j),2);
    RWinc1 = sum(RdW1(:,R*(j-1)+1:R*j),2);
    RWinc2 = sum(RdW2(:,R*(j-1)+1:R*j),2);
    DWinc1 = sum(DdW1(:,R*(j-1)+1:R*j),2);
    DWinc2 = sum(DdW2(:,R*(j-1)+1:R*j),2);

    SJ10 = (0.5*(Dt^(1/2))*(SWinc1+(1/sqrt(3)).*SWinc2));
    I1J10 = (0.5*(Dt^(1/2))*(I1Winc1+(1/sqrt(3)).*I1Winc2));
    I2J10 = (0.5*(Dt^(1/2))*(I2Winc1+(1/sqrt(3)).*I2Winc2));
    I3J10 = (0.5*(Dt^(1/2))*(I3Winc1+(1/sqrt(3)).*I3Winc2));
    I4J10 = (0.5*(Dt^(1/2))*(I4Winc1+(1/sqrt(3)).*I4Winc2));
    RJ10 = (0.5*(Dt^(1/2))*(RWinc1+(1/sqrt(3)).*RWinc2));
    DJ10 = (0.5*(Dt^(1/2))*(DWinc1+(1/sqrt(3)).*DWinc2));

    S1 = Stemp;
    I11 = I1temp;
    I21 = I2temp;
    I31 = I3temp;
    I41 = I4temp;
    R1 = Rtemp;
    D1 = Dtemp;
    %Refer Equation 4.56
    f1_S = fS(S1,I11,I21,I31,I41,R1,D1,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g1_S = gS(S1,I11,I21,I31,I41,R1,D1,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f1_I1 = fI1(S1,I11,I21,I31,I41,R1,D1,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g1_I1 = gI1(S1,I11,I21,I31,I41,R1,D1,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f1_I2 = fI2(S1,I11,I21,I31,I41,R1,D1,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g1_I2 = gI2(S1,I11,I21,I31,I41,R1,D1,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f1_I3 = fI3(S1,I11,I21,I31,I41,R1,D1,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g1_I3 = gI3(S1,I11,I21,I31,I41,R1,D1,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f1_I4 = fI4(S1,I11,I21,I31,I41,R1,D1,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g1_I4 = gI4(S1,I11,I21,I31,I41,R1,D1,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f1_R = fR(S1,I11,I21,I31,I41,R1,D1,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g1_R = gR(S1,I11,I21,I31,I41,R1,D1,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f1_D = fD(S1,I11,I21,I31,I41,R1,D1,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g1_D = gD(S1,I11,I21,I31,I41,R1,D1,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);

    S2 = Stemp + Dt*(A(2)*f1_S) + (B_1(2).*SWinc1 + B_2(2).*SJ10).*g1_S;
    I12 = I1temp + Dt*(A(2)*f1_I1) + (B_1(2).*I1Winc1 + B_2(2).*I1J10).*g1_I1;
    I22 = I2temp + Dt*(A(2)*f1_I2) + (B_1(2).*I2Winc1 + B_2(2).*I2J10).*g1_I2;
    I32 = I3temp + Dt*(A(2)*f1_I3) + (B_1(2).*I3Winc1 + B_2(2).*I3J10).*g1_I3;
    I42 = I4temp + Dt*(A(2)*f1_I4) + (B_1(2).*I4Winc1 + B_2(2).*I4J10).*g1_I4;
    R2 = Rtemp + Dt*(A(2)*f1_R) + (B_1(2).*RWinc1 + B_2(2).*RJ10).*g1_R;
    D2 = Dtemp + Dt*(A(2)*f1_D) + (B_1(2).*DWinc1 + B_2(2).*DJ10).*g1_D;

    %Refer Equation 4.56
    f2_S = fS(S2,I12,I22,I32,I42,R2,D2,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g2_S = gS(S2,I12,I22,I32,I42,R2,D2,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f2_I1 = fI1(S2,I12,I22,I32,I42,R2,D2,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g2_I1 = gI1(S2,I12,I22,I32,I42,R2,D2,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f2_I2 = fI2(S2,I12,I22,I32,I42,R2,D2,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g2_I2 = gI2(S2,I12,I22,I32,I42,R2,D2,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f2_I3 = fI3(S2,I12,I22,I32,I42,R2,D2,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g2_I3 = gI3(S2,I12,I22,I32,I42,R2,D2,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f2_I4 = fI4(S2,I12,I22,I32,I42,R2,D2,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g2_I4 = gI4(S2,I12,I22,I32,I42,R2,D2,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f2_R = fR(S2,I12,I22,I32,I42,R2,D2,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g2_R = gR(S2,I12,I22,I32,I42,R2,D2,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f2_D = fD(S2,I12,I22,I32,I42,R2,D2,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g2_D = gD(S2,I12,I22,I32,I42,R2,D2,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);

    S3 = Stemp + Dt*(A(7)*f2_S) + (B_1(3).*SWinc1 + B_2(3).*SJ10).*g2_S + B_1(7)*SJ10.*g2_S;
    I13 = I1temp + Dt*(A(7)*f2_I1) + (B_1(3).*I1Winc1 + B_2(3).*I1J10).*g2_I1 + B_1(7)*I1J10.*g2_I1;
    I23 = I2temp + Dt*(A(7)*f2_I2) + (B_1(3).*I2Winc1 + B_2(3).*I2J10).*g2_I2 + B_1(7)*I2J10.*g2_I2;
    I33 = I3temp + Dt*(A(7)*f2_I3) + (B_1(3).*I3Winc1 + B_2(3).*I3J10).*g2_I3 + B_1(7)*I3J10.*g2_I3;
    I43 = I4temp + Dt*(A(7)*f2_I3) + (B_1(3).*I4Winc1 + B_2(3).*I4J10).*g2_I4 + B_1(7)*I4J10.*g2_I4;
    R3 = Rtemp + Dt*(A(7)*f2_R) + (B_1(3).*RWinc1 + B_2(3).*RJ10).*g2_R + B_1(7)*RJ10.*g2_R;
    D3 = Dtemp + Dt*(A(7)*f2_D) + (B_1(3).*DWinc1 + B_2(3).*DJ10).*g2_D + B_1(7)*DJ10.*g2_D;

    %Refer Equation 4.56
    f3_S = fS(S3,I13,I23,I33,I43,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g3_S = gS(S3,I13,I23,I33,I43,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f3_I1 = fI1(S3,I13,I23,I33,I43,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g3_I1 = gI1(S3,I13,I23,I33,I43,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f3_I2 = fI2(S3,I13,I23,I33,I43,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g3_I2 = gI2(S3,I13,I23,I33,I43,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f3_I3 = fI3(S3,I13,I23,I33,I43,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g3_I3 = gI3(S3,I13,I23,I33,I43,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f3_I4 = fI4(S3,I13,I23,I33,I43,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g3_I4 = gI4(S3,I13,I23,I33,I43,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f3_R = fR(S3,I13,I23,I33,I43,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g3_R = gR(S3,I13,I23,I33,I43,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f3_D = fD(S3,I13,I23,I33,I43,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g3_D = gD(S3,I13,I23,I33,I43,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);

    S4 = Stemp + Dt*(f3_S) + (B_1(4).*SWinc1 + B_2(4).*SJ10).*g3_S + (B_1(8).*SWinc1.*g3_S + B_1(12).*SWinc1.*g3_S);
    I14 = I1temp + Dt*(f3_I1) + (B_1(4).*I1Winc1 + B_2(4).*I1J10).*g3_I1 + (B_1(8).*I1Winc1.*g3_I1 + B_1(12).*I1Winc1.*g3_I1);
    I24 = I2temp + Dt*(f3_I2) + (B_1(4).*I2Winc1 + B_2(4).*I2J10).*g3_I2 + (B_1(8).*I2Winc1.*g3_I2 + B_1(12).*I2Winc1.*g3_I2);
    I34 = I3temp + Dt*(f3_I3) + (B_1(4).*I3Winc1 + B_2(4).*I3J10).*g3_I3 + (B_1(8).*I3Winc1.*g3_I3 + B_1(12).*I3Winc1.*g3_I3);
    I44 = I4temp + Dt*(f3_I4) + (B_1(4).*I4Winc1 + B_2(4).*I4J10).*g3_I4 + (B_1(8).*I4Winc1.*g3_I4 + B_1(12).*I4Winc1.*g3_I4);
    R4 = Rtemp + Dt*(f3_R) + (B_1(4).*RWinc1 + B_2(4).*RJ10).*g3_R + (B_1(8).*RWinc1.*g3_R + B_1(12).*RWinc1.*g3_R);
    D4 = Dtemp + Dt*(f3_D) + (B_1(4).*DWinc1 + B_2(4).*DJ10).*g3_D + (B_1(8).*DWinc1.*g3_D + B_1(12).*DWinc1.*g3_D);
    %Refer Equation 4.56

    f4_S = fS(S4,I14,I24,I34,I44,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g4_S = gS(S4,I14,I24,I34,I44,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f4_I1 = fI1(S4,I14,I24,I34,I44,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g4_I1 = gI1(S4,I14,I24,I34,I44,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f4_I2 = fI2(S4,I14,I24,I34,I44,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g4_I2 = gI2(S4,I14,I24,I34,I44,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f4_I3 = fI3(S4,I14,I24,I34,I44,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g4_I3 = gI3(S4,I14,I24,I34,I44,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f4_I4 = fI3(S4,I14,I24,I34,I44,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g4_I4 = gI3(S4,I14,I24,I34,I44,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f4_R = fR(S4,I14,I24,I34,I44,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g4_R = gR(S4,I14,I24,I34,I44,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f4_D = fD(S4,I14,I24,I34,I44,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g4_D = gD(S4,I14,I24,I34,I44,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);

    Stemp = Stemp + Dt.*((1/6)*f1_S + (1/3).*f2_S + (1/3).*f3_S + (1/6).*f4_S) + (G_1(1).*SWinc1 + G_2(1).*SJ10).*g1_S + (G_1(2).*SWinc1 + G_2(2).*SJ10).*g2_S + (G_1(3).*SWinc1 + G_2(3).*SJ10).*g3_S + (G_1(4).*SWinc1 + G_2(4).*SJ10).*g4_S;
    I1temp = I1temp + Dt.*((1/6)*f1_I1 + (1/3).*f2_I1 + (1/3).*f3_I1 + (1/6).*f4_I1) + (G_1(1).*I1Winc1 + G_2(1).*I1J10).*g1_I1 + (G_1(2).*I1Winc1 + G_2(2).*I1J10).*g2_I1 + (G_1(3).*I1Winc1 + G_2(3).*I1J10).*g3_I1 + (G_1(4).*I1Winc1 + G_2(4).*I1J10).*g4_I1;
    I2temp = I2temp + Dt.*((1/6)*f1_I2 + (1/3).*f2_I2 + (1/3).*f3_I2 + (1/6).*f4_I2) + (G_1(1).*I2Winc1 + G_2(1).*I2J10).*g1_I2 + (G_1(2).*I2Winc1 + G_2(2).*I2J10).*g2_I2 + (G_1(3).*I2Winc1 + G_2(3).*I2J10).*g3_I2 + (G_1(4).*I2Winc1 + G_2(4).*I2J10).*g4_I2;
    I3temp = I3temp + Dt.*((1/6)*f1_I3 + (1/3).*f2_I3 + (1/3).*f3_I3 + (1/6).*f4_I3) + (G_1(1).*I3Winc1 + G_2(1).*I3J10).*g1_I3 + (G_1(2).*I3Winc1 + G_2(2).*I3J10).*g2_I3 + (G_1(3).*I3Winc1 + G_2(3).*I3J10).*g3_I3 + (G_1(4).*I3Winc1 + G_2(4).*I3J10).*g4_I3;
    I4temp = I4temp + Dt.*((1/6)*f1_I4 + (1/3).*f2_I4 + (1/3).*f3_I4 + (1/6).*f4_I4) + (G_1(1).*I4Winc1 + G_2(1).*I4J10).*g1_I4 + (G_1(2).*I4Winc1 + G_2(2).*I4J10).*g2_I4 + (G_1(3).*I4Winc1 + G_2(3).*I4J10).*g3_I4 + (G_1(4).*I4Winc1 + G_2(4).*I4J10).*g4_I4;
    Rtemp = Rtemp + Dt.*((1/6)*f1_R + (1/3).*f2_R + (1/3).*f3_R + (1/6).*f4_R) + (G_1(1).*RWinc1 + G_2(1).*RJ10).*g1_R + (G_1(2).*RWinc1 + G_2(2).*RJ10).*g2_R + (G_1(3).*RWinc1 + G_2(3).*RJ10).*g3_R + (G_1(4).*RWinc1 + G_2(4).*RJ10).*g4_R;
    Dtemp = Dtemp + Dt.*((1/6)*f1_D + (1/3).*f2_D + (1/3).*f3_D + (1/6).*f4_D) + (G_1(1).*DWinc1 + G_2(1).*DJ10).*g1_D + (G_1(2).*DWinc1 + G_2(2).*DJ10).*g2_D + (G_1(3).*DWinc1 + G_2(3).*DJ10).*g3_D + (G_1(4).*DWinc1 + G_2(4).*DJ10).*g4_D;

    %Refer Equation 4.56
    Ssrk4_2(:,j) = Stemp.';
    I1srk4_2(:,j) = I1temp.';
    I2srk4_2(:,j) = I2temp.';
    I3srk4_2(:,j) = I3temp.';
    I4srk4_2(:,j) = I4temp.';
    Rsrk4_2(:,j) = Rtemp.';
    Dsrk4_2(:,j) = Dtemp.';
end
for i = 1:L
   Ssrk_3 = Ssrk4_2(:,i);
   I1srk_3 = I1srk4_2(:,i);
   I2srk_3 = I2srk4_2(:,i);
   I3srk_3 = I3srk4_2(:,i);
   I4srk_3 = I4srk4_2(:,i);
   Rsrk_3 = Rsrk4_2(:,i);
   Dsrk_3 = Dsrk4_2(:,i);
   Smodel(1,i) = interp1(linspace(1/numel(Ssrk_3),1,numel(Ssrk_3)), sort(Ssrk_3), 0.25);
   I1model(1,i) = interp1(linspace(1/numel(I1srk_3),1,numel(I1srk_3)), sort(I1srk_3), 0.25);
   I2model(1,i) = interp1(linspace(1/numel(I2srk_3),1,numel(I2srk_3)), sort(I2srk_3), 0.25);
   I3model(1,i) = interp1(linspace(1/numel(I3srk_3),1,numel(I3srk_3)), sort(I3srk_3), 0.25);
   I4model(1,i) = interp1(linspace(1/numel(I4srk_3),1,numel(I4srk_3)), sort(I4srk_3), 0.25);
   Rmodel(1,i) = interp1(linspace(1/numel(Rsrk_3),1,numel(Rsrk_3)), sort(Rsrk_3), 0.25);
   Dmodel(1,i) = interp1(linspace(1/numel(Dsrk_3),1,numel(Dsrk_3)), sort(Dsrk_3), 0.25);
   Smodel(2,i) = interp1(linspace(1/numel(Ssrk_3),1,numel(Ssrk_3)), sort(Ssrk_3), 0.50);
   I1model(2,i) = interp1(linspace(1/numel(I1srk_3),1,numel(I1srk_3)), sort(I1srk_3), 0.50);
   I2model(2,i) = interp1(linspace(1/numel(I2srk_3),1,numel(I2srk_3)), sort(I2srk_3), 0.50);
   I3model(2,i) = interp1(linspace(1/numel(I3srk_3),1,numel(I3srk_3)), sort(I3srk_3), 0.50);
   I4model(2,i) = interp1(linspace(1/numel(I4srk_3),1,numel(I4srk_3)), sort(I4srk_3), 0.50);
   Rmodel(2,i) = interp1(linspace(1/numel(Rsrk_3),1,numel(Rsrk_3)), sort(Rsrk_3), 0.50);
   Dmodel(2,i) = interp1(linspace(1/numel(Dsrk_3),1,numel(Dsrk_3)), sort(Dsrk_3), 0.50);
   Smodel(3,i) = interp1(linspace(1/numel(Ssrk_3),1,numel(Ssrk_3)), sort(Ssrk_3), 0.75);
   I1model(3,i) = interp1(linspace(1/numel(I1srk_3),1,numel(I1srk_3)), sort(I1srk_3), 0.75);
   I2model(3,i) = interp1(linspace(1/numel(I2srk_3),1,numel(I2srk_3)), sort(I2srk_3), 0.75);
   I3model(3,i) = interp1(linspace(1/numel(I3srk_3),1,numel(I3srk_3)), sort(I3srk_3), 0.75);
   I4model(3,i) = interp1(linspace(1/numel(I4srk_3),1,numel(I4srk_3)), sort(I4srk_3), 0.75);
   Rmodel(3,i) = interp1(linspace(1/numel(Rsrk_3),1,numel(Rsrk_3)), sort(Rsrk_3), 0.75);
   Dmodel(3,i) = interp1(linspace(1/numel(Dsrk_3),1,numel(Dsrk_3)), sort(Dsrk_3), 0.75);
end
Save_srk4 = Smodel(2,:);
I1ave_srk4 = I1model(2,:);
I2ave_srk4 = I2model(2,:);
I3ave_srk4 = I3model(2,:);
I4ave_srk4 = I4model(2,:);
Rave_srk4 = Rmodel(2,:);
Dave_srk4 = Dmodel(2,:);
%% RK4==================================================
sigma = 0;
Srk4_2 = zeros(1,N);
I1rk4_2 = zeros(1,N);
I2rk4_2 = zeros(1,N);
I3rk4_2 = zeros(1,N);
I4rk4_2 = zeros(1,N);
Rrk4_2 = zeros(1,N);
Drk4_2 = zeros(1,N);
Stemp = Szero*ones(1,1);
I1temp = I1zero*ones(1,1);
I2temp = I2zero*ones(1,1);
I3temp = I3zero*ones(1,1);
I4temp = I4zero*ones(1,1);
Rtemp = Rzero*ones(1,1);
Dtemp = Dzero*ones(1,1);

for j = 1:L
    p1 = Pastparamestimate(j,2);
    p2 = Pastparamestimate(j,3);
    p3 = Pastparamestimate(j,4);
    p4 = Pastparamestimate(j,5);
    Beta1 = Pastparamestimate(j,6);
    Beta2 = Pastparamestimate(j,7);
    Beta3 = Pastparamestimate(j,8);
    Beta4 = Pastparamestimate(j,9);
    r1 = Pastparamestimate(j,10);
    r2 = Pastparamestimate(j,11);
    r3 = Pastparamestimate(j,12);
    r4 = Pastparamestimate(j,13);
    d1 = Pastparamestimate(j,14);
    d2 = Pastparamestimate(j,15);
    d3 = Pastparamestimate(j,16);
    d4 = Pastparamestimate(j,17);
    gamma = Pastparamestimate(j,18);
    
    ptotal = p1 + p2 + p3 + p4;
    p1 = p1/ptotal*p_vax;
    p2 = p2/ptotal*p_vax;
    p3 = p3/ptotal*p_vax;
    p4 = p4/ptotal*p_vax;

    S1 = Stemp;
    I11 = I1temp;
    I21 = I2temp;
    I31 = I3temp;
    I41 = I4temp;
    R1 = Rtemp;
    D1 = Dtemp;
    %Refer Equation 4.56
    f1_S = fS(S1,I11,I21,I31,I41,R1,D1,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g1_S = gS(S1,I11,I21,I31,I41,R1,D1,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f1_I1 = fI1(S1,I11,I21,I31,I41,R1,D1,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g1_I1 = gI1(S1,I11,I21,I31,I41,R1,D1,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f1_I2 = fI2(S1,I11,I21,I31,I41,R1,D1,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g1_I2 = gI2(S1,I11,I21,I31,I41,R1,D1,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f1_I3 = fI3(S1,I11,I21,I31,I41,R1,D1,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g1_I3 = gI3(S1,I11,I21,I31,I41,R1,D1,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f1_I4 = fI4(S1,I11,I21,I31,I41,R1,D1,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g1_I4 = gI4(S1,I11,I21,I31,I41,R1,D1,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f1_R = fR(S1,I11,I21,I31,I41,R1,D1,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g1_R = gR(S1,I11,I21,I31,I41,R1,D1,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f1_D = fD(S1,I11,I21,I31,I41,R1,D1,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g1_D = gD(S1,I11,I21,I31,I41,R1,D1,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);

    S2 = Stemp + Dt*(A(2)*f1_S);
    I12 = I1temp + Dt*(A(2)*f1_I1);
    I22 = I2temp + Dt*(A(2)*f1_I2);
    I32 = I3temp + Dt*(A(2)*f1_I3);
    I42 = I4temp + Dt*(A(2)*f1_I4);
    R2 = Rtemp + Dt*(A(2)*f1_R);
    D2 = Dtemp + Dt*(A(2)*f1_D);

    %Refer Equation 4.56
    f2_S = fS(S2,I12,I22,I32,I42,R2,D2,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g2_S = gS(S2,I12,I22,I32,I42,R2,D2,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f2_I1 = fI1(S2,I12,I22,I32,I42,R2,D2,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g2_I1 = gI1(S2,I12,I22,I32,I42,R2,D2,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f2_I2 = fI2(S2,I12,I22,I32,I42,R2,D2,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g2_I2 = gI2(S2,I12,I22,I32,I42,R2,D2,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f2_I3 = fI3(S2,I12,I22,I32,I42,R2,D2,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g2_I3 = gI3(S2,I12,I22,I32,I42,R2,D2,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f2_I4 = fI4(S2,I12,I22,I32,I42,R2,D2,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g2_I4 = gI4(S2,I12,I22,I32,I42,R2,D2,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f2_R = fR(S2,I12,I22,I32,I42,R2,D2,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g2_R = gR(S2,I12,I22,I32,I42,R2,D2,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f2_D = fD(S2,I12,I22,I32,I42,R2,D2,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g2_D = gD(S2,I12,I22,I32,I42,R2,D2,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);

    S3 = Stemp + Dt*(A(7)*f2_S);
    I13 = I1temp + Dt*(A(7)*f2_I1);
    I23 = I2temp + Dt*(A(7)*f2_I2);
    I33 = I3temp + Dt*(A(7)*f2_I3);
    I43 = I4temp + Dt*(A(7)*f2_I4);
    R3 = Rtemp + Dt*(A(7)*f2_R);
    D3 = Dtemp + Dt*(A(7)*f2_D);

    %Refer Equation 4.56
    f3_S = fS(S3,I13,I23,I33,I43,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g3_S = gS(S3,I13,I23,I33,I43,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f3_I1 = fI1(S3,I13,I23,I33,I43,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g3_I1 = gI1(S3,I13,I23,I33,I43,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f3_I2 = fI2(S3,I13,I23,I33,I43,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g3_I2 = gI2(S3,I13,I23,I33,I43,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f3_I3 = fI3(S3,I13,I23,I33,I43,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g3_I3 = gI3(S3,I13,I23,I33,I43,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f3_I4 = fI4(S3,I13,I23,I33,I43,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g3_I4 = gI4(S3,I13,I23,I33,I43,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f3_R = fR(S3,I13,I23,I33,I43,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g3_R = gR(S3,I13,I23,I33,I43,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f3_D = fD(S3,I13,I23,I33,I43,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g3_D = gD(S3,I13,I23,I33,I43,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);

    S4 = Stemp + Dt*(f3_S);
    I14 = I1temp + Dt*(f3_I1);
    I24 = I2temp + Dt*(f3_I2);
    I34 = I3temp + Dt*(f3_I3);
    I44 = I4temp + Dt*(f3_I4);
    R4 = Rtemp + Dt*(f3_R);
    D4 = Dtemp + Dt*(f3_D);
    %Refer Equation 4.56

    f4_S = fS(S4,I14,I24,I34,I44,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g4_S = gS(S4,I14,I24,I34,I44,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f4_I1 = fI1(S4,I14,I24,I34,I44,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g4_I1 = gI1(S4,I14,I24,I34,I44,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f4_I2 = fI2(S4,I14,I24,I34,I44,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g4_I2 = gI2(S4,I14,I24,I34,I44,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f4_I3 = fI3(S4,I14,I24,I34,I44,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g4_I3 = gI3(S4,I14,I24,I34,I44,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f4_I4 = fI3(S4,I14,I24,I34,I44,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g4_I4 = gI3(S4,I14,I24,I34,I44,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f4_R = fR(S4,I14,I24,I34,I44,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g4_R = gR(S4,I14,I24,I34,I44,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f4_D = fD(S4,I14,I24,I34,I44,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    g4_D = gD(S4,I14,I24,I34,I44,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);

    Stemp = Stemp + Dt.*((1/6)*f1_S + (1/3).*f2_S + (1/3).*f3_S + (1/6).*f4_S);
    I1temp = I1temp + Dt.*((1/6)*f1_I1 + (1/3).*f2_I1 + (1/3).*f3_I1 + (1/6).*f4_I1);
    I2temp = I2temp + Dt.*((1/6)*f1_I2 + (1/3).*f2_I2 + (1/3).*f3_I2 + (1/6).*f4_I2);
    I3temp = I3temp + Dt.*((1/6)*f1_I3 + (1/3).*f2_I3 + (1/3).*f3_I3 + (1/6).*f4_I3);
    I4temp = I4temp + Dt.*((1/6)*f1_I4 + (1/3).*f2_I4 + (1/3).*f3_I4 + (1/6).*f4_I4);
    Rtemp = Rtemp + Dt.*((1/6)*f1_R + (1/3).*f2_R + (1/3).*f3_R + (1/6).*f4_R);
    Dtemp = Dtemp + Dt.*((1/6)*f1_D + (1/3).*f2_D + (1/3).*f3_D + (1/6).*f4_D);
    
    %Refer Equation 4.56
    Srk4_2(:,j) = Stemp.';
    I1rk4_2(:,j) = I1temp.';
    I2rk4_2(:,j) = I2temp.';
    I3rk4_2(:,j) = I3temp.';
    I4rk4_2(:,j) = I4temp.';
    Rrk4_2(:,j) = Rtemp.';
    Drk4_2(:,j) = Dtemp.';
end
% RMSE_srk4x = rmse(Xsrk4_2,Datax);
% RMSE_srk4y = rmse(Ysrk4_2,Datay);
Smodel_2 = mean(Srk4_2,1);
I1model_2 = mean(I1rk4_2,1);
I2model_2 = mean(I2rk4_2,1);
I3model_2 = mean(I3rk4_2,1);
I4model_2 = mean(I4rk4_2,1);
Rmodel_2 = mean(Rrk4_2,1);
Dmodel_2 = mean(Drk4_2,1);
%% Adam_Bashford_Predictor_Corrector=====================
Sabpc_2 = Szero*ones(1,N+1);
I1abpc_2 = I1zero*ones(1,N+1);
I2abpc_2 = I2zero*ones(1,N+1);
I3abpc_2 = I3zero*ones(1,N+1);
I4abpc_2 = I4zero*ones(1,N+1);
Rabpc_2 = Rzero*ones(1,N+1);
Dabpc_2 = Dzero*ones(1,N+1);

pcS1 = 0; pcS2 = 0; pcS3 = 0;
pcI11 = 0; pcI12 = 0; pcI13 = 0;
pcI21 = 0; pcI22 = 0; pcI23 = 0;
pcI31 = 0; pcI32 = 0; pcI33 = 0;
pcI41 = 0; pcI42 = 0; pcI43 = 0;
pcR1 = 0; pcR2 = 0; pcR3 = 0;
pcD1 = 0; pcD2 = 0; pcD3 = 0;

for j = 1:L
    p1 = Pastparamestimate(j,2);
    p2 = Pastparamestimate(j,3);
    p3 = Pastparamestimate(j,4);
    p4 = Pastparamestimate(j,5);
    Beta1 = Pastparamestimate(j,6);
    Beta2 = Pastparamestimate(j,7);
    Beta3 = Pastparamestimate(j,8);
    Beta4 = Pastparamestimate(j,9);
    r1 = Pastparamestimate(j,10);
    r2 = Pastparamestimate(j,11);
    r3 = Pastparamestimate(j,12);
    r4 = Pastparamestimate(j,13);
    d1 = Pastparamestimate(j,14);
    d2 = Pastparamestimate(j,15);
    d3 = Pastparamestimate(j,16);
    d4 = Pastparamestimate(j,17);
    gamma = Pastparamestimate(j,18);
    
    ptotal = p1 + p2 + p3 + p4;
    p1 = p1/ptotal*p_vax;
    p2 = p2/ptotal*p_vax;
    p3 = p3/ptotal*p_vax;
    p4 = p4/ptotal*p_vax;
    
    S1 = Sabpc_2(j);
    I11 = I1abpc_2(j);
    I21 = I2abpc_2(j);
    I31 = I3abpc_2(j);
    I41 = I4abpc_2(j);
    R1 = Rabpc_2(j);
    D1 = Dabpc_2(j);

    f1_S = fS(S1,I11,I21,I31,I41,R1,D1,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f1_I1 = fI1(S1,I11,I21,I31,I41,R1,D1,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f1_I2 = fI2(S1,I11,I21,I31,I41,R1,D1,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f1_I3 = fI3(S1,I11,I21,I31,I41,R1,D1,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f1_I4 = fI4(S1,I11,I21,I31,I41,R1,D1,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f1_R = fR(S1,I11,I21,I31,I41,R1,D1,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f1_D = fD(S1,I11,I21,I31,I41,R1,D1,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);

    pcS0 = f1_S;
    pcI10 = f1_I1;
    pcI20 = f1_I2;
    pcI30 = f1_I3;
    pcI40 = f1_I4;
    pcR0 = f1_R;
    pcD0 = f1_D;

    S2 = Sabpc_2(j) + dt*(0.5*f1_S);
    I12 = I1abpc_2(j) + dt*(0.5*f1_I1);
    I22 = I2abpc_2(j) + dt*(0.5*f1_I2);
    I32 = I3abpc_2(j) + dt*(0.5*f1_I3);
    I42 = I4abpc_2(j) + dt*(0.5*f1_I4);
    R2 = Rabpc_2(j) + dt*(0.5*f1_R);
    D2 = Dabpc_2(j) + dt*(0.5*f1_D);

    f2_S = fS(S2,I12,I22,I32,I42,R2,D2,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f2_I1 = fI1(S2,I12,I22,I32,I42,R2,D2,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f2_I2 = fI2(S2,I12,I22,I32,I42,R2,D2,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f2_I3 = fI3(S2,I12,I22,I32,I42,R2,D2,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f2_I4 = fI4(S2,I12,I22,I32,I42,R2,D2,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f2_R = fR(S2,I12,I22,I32,I42,R2,D2,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f2_D = fD(S2,I12,I22,I32,I42,R2,D2,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);

    S3 = Sabpc_2(j) + dt*(0.5*f2_S);
    I13 = I1abpc_2(j) + dt*(0.5*f2_I1);
    I23 = I2abpc_2(j) + dt*(0.5*f2_I2);
    I33 = I3abpc_2(j) + dt*(0.5*f2_I3);
    I43 = I4abpc_2(j) + dt*(0.5*f2_I4);
    R3 = Rabpc_2(j) + dt*(0.5*f2_R);
    D3 = Dabpc_2(j) + dt*(0.5*f2_D);

    f3_S = fS(S3,I13,I23,I33,I43,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f3_I1 = fI1(S3,I13,I23,I33,I43,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f3_I2 = fI2(S3,I13,I23,I33,I43,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f3_I3 = fI3(S3,I13,I23,I33,I43,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f3_I4 = fI4(S3,I13,I23,I33,I43,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f3_R = fR(S3,I13,I23,I33,I43,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f3_D = fD(S3,I13,I23,I33,I43,R3,D3,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);

    S4 = Sabpc_2(j) + dt*f3_S;
    I14 = I1abpc_2(j) + dt*f3_I1;
    I24 = I2abpc_2(j) + dt*f3_I2;
    I34 = I3abpc_2(j) + dt*f3_I3;
    I44 = I4abpc_2(j) + dt*f3_I4;
    R4 = Rabpc_2(j) + dt*f3_R;
    D4 = Dabpc_2(j) + dt*f3_D;

    f4_S = fS(S4,I14,I24,I34,I44,R4,D4,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f4_I1 = fI1(S4,I14,I24,I34,I44,R4,D4,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f4_I2 = fI2(S4,I14,I24,I34,I44,R4,D4,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f4_I3 = fI3(S4,I14,I24,I34,I44,R4,D4,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f4_I4 = fI4(S4,I14,I24,I34,I44,R4,D4,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f4_R = fR(S4,I14,I24,I34,I44,R4,D4,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f4_D = fD(S4,I14,I24,I34,I44,R4,D4,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);

    Sabpc_2(j+1) = Sabpc_2(j) + dt*((1/6)*f1_S + (1/3)*f2_S + (1/3)*f3_S + (1/6)*f4_S);
    I1abpc_2(j+1) = I1abpc_2(j) + dt*((1/6)*f1_I1 + (1/3)*f2_I1 + (1/3)*f3_I1 + (1/6)*f4_I1);
    I2abpc_2(j+1) = I2abpc_2(j) + dt*((1/6)*f1_I2 + (1/3)*f2_I2 + (1/3)*f3_I2 + (1/6)*f4_I2);
    I3abpc_2(j+1) = I3abpc_2(j) + dt*((1/6)*f1_I3 + (1/3)*f2_I3 + (1/3)*f3_I3 + (1/6)*f4_I3);
    I4abpc_2(j+1) = I4abpc_2(j) + dt*((1/6)*f1_I4 + (1/3)*f2_I4 + (1/3)*f3_I4 + (1/6)*f4_I4);
    Rabpc_2(j+1) = Rabpc_2(j) + dt*((1/6)*f1_R + (1/3)*f2_R + (1/3)*f3_R + (1/6)*f4_R);
    Dabpc_2(j+1) = Dabpc_2(j) + dt*((1/6)*f1_D + (1/3)*f2_D + (1/3)*f3_D + (1/6)*f4_D);

    pcS1 = pcS0; pcS2 = pcS1; pcS3 = pcS2;
    pcI11 = pcI10; pcI12 = pcI11; pcI13 = pcI12;
    pcI21 = pcI20; pcI22 = pcI21; pcI23 = pcI22;
    pcI31 = pcI30; pcI32 = pcI31; pcI33 = pcI32;
    pcI41 = pcI40; pcI42 = pcI41; pcI43 = pcI42;
    pcR1 = pcR0; pcR2 = pcR1; pcR3 = pcR2;
    pcD1 = pcD0; pcD2 = pcD1; pcD3 = pcD2;
end
        
for j = 1:L
    p1 = Pastparamestimate(j,2);
    p2 = Pastparamestimate(j,3);
    p3 = Pastparamestimate(j,4);
    p4 = Pastparamestimate(j,5);
    Beta1 = Pastparamestimate(j,6);
    Beta2 = Pastparamestimate(j,7);
    Beta3 = Pastparamestimate(j,8);
    Beta4 = Pastparamestimate(j,9);
    r1 = Pastparamestimate(j,10);
    r2 = Pastparamestimate(j,11);
    r3 = Pastparamestimate(j,12);
    r4 = Pastparamestimate(j,13);
    d1 = Pastparamestimate(j,14);
    d2 = Pastparamestimate(j,15);
    d3 = Pastparamestimate(j,16);
    d4 = Pastparamestimate(j,17);
    gamma = Pastparamestimate(j,18);
    
    p1 = p1/(p1 + p2 + p3 + p4);
    p2 = p2/(p1 + p2 + p3 + p4);
    p3 = p3/(p1 + p2 + p3 + p4);
    p4 = p4/(p1 + p2 + p3 + p4);
    
    S1 = Sabpc_2(j);
    I11 = I1abpc_2(j);
    I21 = I2abpc_2(j);
    I31 = I3abpc_2(j);
    I41 = I4abpc_2(j);
    R1 = Rabpc_2(j);
    D1 = Dabpc_2(j);

    f1_S = fS(S1,I11,I21,I31,I41,R1,D1,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f1_I1 = fI1(S1,I11,I21,I31,I41,R1,D1,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f1_I2 = fI2(S1,I11,I21,I31,I41,R1,D1,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f1_I3 = fI3(S1,I11,I21,I31,I41,R1,D1,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f1_I4 = fI4(S1,I11,I21,I31,I41,R1,D1,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f1_R = fR(S1,I11,I21,I31,I41,R1,D1,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    f1_D = fD(S1,I11,I21,I31,I41,R1,D1,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);

    pcS0 = f1_S;
    pcI10 = f1_I1;
    pcI20 = f1_I2;
    pcI30 = f1_I3;
    pcI40 = f1_I4;
    pcR0 = f1_R;
    pcD0 = f1_D;

    wS = Sabpc_2(j) + dt*( 55.0 * pcS0 - 59.0 * pcS1 + 37.0 * pcS2 - 9.0 * pcS3 ) / 24.0;
    wI1 = I1abpc_2(j) + dt*( 55.0 * pcI10 - 59.0 * pcI11 + 37.0 * pcI12 - 9.0 * pcI13 ) / 24.0;
    wI2 = I2abpc_2(j) + dt*( 55.0 * pcI20 - 59.0 * pcI21 + 37.0 * pcI22 - 9.0 * pcI23 ) / 24.0;
    wI3 = I3abpc_2(j) + dt*( 55.0 * pcI30 - 59.0 * pcI31 + 37.0 * pcI32 - 9.0 * pcI33 ) / 24.0;
    wI4 = I4abpc_2(j) + dt*( 55.0 * pcI40 - 59.0 * pcI41 + 37.0 * pcI42 - 9.0 * pcI43 ) / 24.0;
    wR = Rabpc_2(j) + dt*( 55.0 * pcR0 - 59.0 * pcR1 + 37.0 * pcR2 - 9.0 * pcR3 ) / 24.0;
    wD = Dabpc_2(j) + dt*( 55.0 * pcD0 - 59.0 * pcD1 + 37.0 * pcD2 - 9.0 * pcD3 ) / 24.0;

    fw_S = fS(wS,wI1,wI2,wI3,wI4,wR,wD,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    fw_I1 = fI1(wS,wI1,wI2,wI3,wI4,wR,wD,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    fw_I2 = fI2(wS,wI1,wI2,wI3,wI4,wR,wD,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    fw_I3 = fI3(wS,wI1,wI2,wI3,wI4,wR,wD,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    fw_I4 = fI4(wS,wI1,wI2,wI3,wI4,wR,wD,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    fw_R = fR(wS,wI1,wI2,wI3,wI4,wR,wD,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    fw_D = fD(wS,wI1,wI2,wI3,wI4,wR,wD,Beta1,Beta2,Beta3,Beta4,r1,r2,r3,r4,d1,d2,d3,d4,p1,p2,p3,p4,gamma,sigma,mu);
    
    Sabpc_2(j+1) = Sabpc_2(j) + dt*( 9.0 * fw_S + 19.0 * pcS0 - 5.0 * pcS1 + pcS2 ) / 24.0;
    I1abpc_2(j+1) = I1abpc_2(j) + dt*( 9.0 * fw_I1 + 19.0 * pcI10 - 5.0 * pcI11 + pcI12 ) / 24.0;
    I2abpc_2(j+1) = I2abpc_2(j) + dt*( 9.0 * fw_I2 + 19.0 * pcI20 - 5.0 * pcI21 + pcI22 ) / 24.0;
    I3abpc_2(j+1) = I3abpc_2(j) + dt*( 9.0 * fw_I3 + 19.0 * pcI30 - 5.0 * pcI31 + pcI32 ) / 24.0;
    I4abpc_2(j+1) = I4abpc_2(j) + dt*( 9.0 * fw_I4 + 19.0 * pcI40 - 5.0 * pcI41 + pcI42 ) / 24.0;
    Rabpc_2(j+1) = Rabpc_2(j) + dt*( 9.0 * fw_R + 19.0 * pcR0 - 5.0 * pcR1 + pcR2 ) / 24.0;
    Dabpc_2(j+1) = Dabpc_2(j) + dt*( 9.0 * fw_D + 19.0 * pcD0 - 5.0 * pcD1 + pcD2 ) / 24.0;

    pcS1 = pcS0; pcS2 = pcS1; pcS3 = pcS2;
    pcI11 = pcI10; pcI12 = pcI11; pcI13 = pcI12;
    pcI21 = pcI20; pcI22 = pcI21; pcI23 = pcI22;
    pcI31 = pcI30; pcI32 = pcI31; pcI33 = pcI32;
    pcI41 = pcI40; pcI32 = pcI41; pcI43 = pcI42;
    pcR1 = pcR0; pcR2 = pcR1; pcR3 = pcR2;
    pcD1 = pcD0; pcD2 = pcD1; pcD3 = pcD2;
end
% RMSE_srk4x = rmse(Xsrk4_2,Datax);
% RMSE_srk4y = rmse(Ysrk4_2,Datay);
Save_abpc = mean(Sabpc_2,1);
I1ave_abpc = mean(I1abpc_2,1);
I2ave_abpc = mean(I2abpc_2,1);
I3ave_abpc = mean(I3abpc_2,1);
I4ave_abpc = mean(I4abpc_2,1);
Rave_abpc = mean(Rabpc_2,1);
Dave_abpc = mean(Dabpc_2,1);
% Smodel = Save_srk4;
SDet = Save_abpc(2:end);
I1Det = I1ave_abpc(2:end);
I2Det = I2ave_abpc(2:end);
I3Det = I3ave_abpc(2:end);
I4Det = I4ave_abpc(2:end);
RDet = Rave_abpc(2:end);
DDet = Dave_abpc(2:end);
%% Plot
figure('Position', [1, 1, 1600, 1000])
tlayout = tiledlayout(1,1);
tlayout.TileSpacing = 'compact';
tlayout.Padding = 'compact';
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
%     legend('Susceptible Data','Simulation RK4')
%     xlabel('time (date)'); ylabel('Susceptible Population')    
%     
%     ax2 = nexttile; 
%     plot(tdays,I1data,"-s",'Color',"#2C00D2",'MarkerIndices',1:20:length(tdays)); hold on;
%     plot(tdays,I1model_2,':','Color',"#273DE3"); hold on;
%     plot(tdays,I1Det,'--', 'Color',"#4BA3E5")
%     plot(tdays,I1model(2,:),'-','Color',"#77DFE9"); hold on;
%     cax2 = [39/255,61/255,227/255];
%     fill([tdays, fliplr(tdays)], [I1model(1,:), fliplr(I1model(3,:))], cax2, 'LineStyle', "None",'FaceAlpha',0.2); hold on;
%     title(ax2,'Infected Population, I1')
%     ax2.FontSize = 8;
%     legend('Infected Data','RK4','Adam-Moulton PC','SRK4','Location','best')
%     hold off;
%     xlabel('time (date)'); ylabel('Infected Population')    
    
%     ax3 = nexttile; 
%     plot(tdays,I2data,'-s','Color',"#4B0082",'MarkerIndices',1:20:length(tdays)); hold on;
%     plot(tdays,I2model_2,':','Color',"#7B68EE"); hold on;
%     plot(tdays,I2Det,'--', 'Color',"#8B008B")
%     plot(tdays,I2model(2,:),'-','Color',"#FF00FF"); hold on;
%     cax3 = [255/255, 0/255, 255/255];
%     fill([tdays, fliplr(tdays)], [I2model(1,:), fliplr(I2model(3,:))], cax3, 'LineStyle', "None",'FaceAlpha',0.2); hold on;
%     title(ax3,'Hospitalised Population, I2')
%     ax3.FontSize = 8;
%     legend('Susceptible Data','RK4','Adam-Moulton PC','SRK4','Location','best')
%     hold off;
%     xlabel('time (date)'); ylabel('Hospitalised Population')    
    
    ax4 = nexttile; 
    plot(tdays,I3data,'-s','Color',"#FF4500",'MarkerIndices',1:20:length(tdays)); hold on;
    plot(tdays,I3model_2,':','Color',"#FF7F50"); hold on;
    plot(tdays,I3Det,'--', 'Color',"#FFA500")
    plot(tdays,I3model(2,:),'-','Color',"#FFD700"); hold on;
    title(ax4,'ICU Population, I3')
    cax4 = [204/255,141/255,61/255];
    fill([tdays, fliplr(tdays)], [I3model(1,:), fliplr(I3model(3,:))], cax4, 'LineStyle', "None",'FaceAlpha',0.2); hold on;
    ax4.FontSize = 10;
    legend('Infected Data','RK4','Adam-Moulton PC','SRK4','Location','best')
    hold off;
    xlabel('time (date)'); ylabel('ICU Population')
    
    maxValue = max(I3model(3,end-299:end));

    % Create the inset axes
    insetAx = axes('Position', [0.5, 0.4, 0.43, 0.45]); % Adjust position and size
    box on; % Draw a box around the inset plot
    hold on;

    % Plot the data in the inset
    plot(tdays, I3data, '-s', 'Color', "#FF4500", 'MarkerIndices', 1:20:length(tdays));
    plot(tdays, I3model_2, ':', 'Color', "#FF7F50");
    plot(tdays, I3Det, '--', 'Color', "#FFA500");
    plot(tdays, I3model(2,:), '-', 'Color', "#FFD700");
    fill([tdays, fliplr(tdays)], [I3model(1,:), fliplr(I3model(3,:))], cax4, 'LineStyle', "None", 'FaceAlpha', 0.2);
    xlim([tdays(end-299), tdays(end)]);
    xlabel('time, last 300 days (date)'); ylabel('ICU Population')    
    
    % Add a horizontal dotted line at the maximum value
    yline(maxValue, 'd--', sprintf('Max Q3: %.2f', maxValue), 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom');

    hold off;
    set(insetAx, 'FontSize', 10);
    
%     ax5 = nexttile; 
%     plot(tdays,I4data,'-s','Color',"#731bc4",'MarkerIndices',1:20:length(tdays)); hold on;
%     plot(tdays,I4model_2,':','Color',"#b169c9"); hold on;
%     plot(tdays,I4Det,'--', 'Color',"#50255e")
%     plot(tdays,I4model(2,:),'-','Color',"#470b9c"); hold on;
%     title(ax5,'Ventilated Population, I4')
%     cax4 = [167/255,116/255,237/255];
%     fill([tdays, fliplr(tdays)], [I4model(1,:), fliplr(I4model(3,:))], cax4, 'LineStyle', "None",'FaceAlpha',0.2); hold on;
%     ax5.FontSize = 12;
%     legend('Infected Data','RK4','Adam-Moulton PC','SRK4','Location','best')
%     hold off;
%     xlabel('time (date)'); ylabel('Ventilated Population')    
%     
%     maxValue = max(I4model(3,end-299:end));
% 
%     % Create the inset axes
%     insetAx = axes('Position', [0.5, 0.4, 0.43, 0.45]); % Adjust position and size
%     box on; % Draw a box around the inset plot
%     hold on;
% 
%     % Plot the data in the inset
%     plot(tdays, I4data, '-s', 'Color', "#731bc4", 'MarkerIndices', 1:20:length(tdays));
%     plot(tdays, I4model_2, ':', 'Color', "#b169c9");
%     plot(tdays, I4Det, '--', 'Color', "#50255e");
%     plot(tdays, I4model(2,:), '-', 'Color', "#470b9c");
%     fill([tdays, fliplr(tdays)], [I4model(1,:), fliplr(I4model(3,:))], cax4, 'LineStyle', "None", 'FaceAlpha', 0.2);
%     xlim([tdays(end-299), tdays(end)]);
%     xlabel('time, last 300 days (date)'); ylabel('Ventilated Population')    
%     
%     % Add a horizontal dotted line at the maximum value
%     yline(maxValue, 'd--', sprintf('Max Q3: %.2f', maxValue), 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom');
% 
%     hold off;
%     set(insetAx, 'FontSize', 12);
%     ax6 = nexttile;
%     plot(tdays,Rdata,'-s','Color',"#008080",'MarkerIndices',1:20:length(tdays)); hold on;
%     plot(tdays,Rmodel_2,':','Color',"#556B2F");hold on;
%     plot(tdays,RDet,'--', 'Color',"#3CB371")
%     plot(tdays,Rmodel(2,:),'-','Color',"#00FF00");hold on;
%     title(ax6,'Recovered Population, R')
%     cax5 = [70/255,192/255,30/255];
%     fill([tdays, fliplr(tdays)], [Rmodel(1,:), fliplr(Rmodel(3,:))], cax5, 'LineStyle', "None",'FaceAlpha',0.2); hold on;
%     ax6.FontSize = 8;
%     legend('Recovered Data','RK4','Adam-Moulton PC','SRK4','Location','best')
%     hold off;
%     xlabel('time (date)'); ylabel('Recovered Population')    
    
%     ax7 = nexttile;
%     plot(tdays,Ddata,'-s','Color',"#8B0000",'MarkerIndices',1:20:length(tdays)); hold on;
%     plot(tdays,Dmodel_2,':','Color',"#FF0000");hold on;
%     plot(tdays,DDet,'--', 'Color',"#8B4513")
%     plot(tdays,Dmodel(2,:),'-','Color',"#F08080");hold on;
%     title(ax7,'Death Population, D')
%     cax6 = [217/255,46/255,157/255];
%     fill([tdays, fliplr(tdays)], [Dmodel(1,:), fliplr(Dmodel(3,:))], cax6, 'LineStyle', "None",'FaceAlpha',0.2); hold on;
%     ax7.FontSize = 8;
%     legend('Death Data','RK4','Adam-Moulton PC','SRK4','Location','best')
%     hold off;
%     xlabel('time (date)'); ylabel('Death Population') 
    

    
    print -depsc si3rdMalaysia.eps
%% RMSE
RMSE = zeros(3,6);
RMSE(3,1) = sqrt(sum((Save_srk4(:) - Sdata(:)).^2)/length(Save_srk4));
RMSE(1,1) = sqrt(sum((Smodel_2(:) - Sdata(:)).^2)/length(Save_srk4));
RMSE(2,1) = sqrt(sum((SDet(:) - Sdata(:)).^2)/length(Save_srk4));
RMSE(3,2) = sqrt(sum((I1ave_srk4(:) - I1data(:)).^2)/length(Save_srk4));
RMSE(1,2) = sqrt(sum((I1model_2(:) - I1data(:)).^2)/length(Save_srk4));
RMSE(2,2) = sqrt(sum((I1Det(:) - I1data(:)).^2)/length(Save_srk4));
RMSE(3,3) = sqrt(sum((I2ave_srk4(:) - I2data(:)).^2)/length(Save_srk4));
RMSE(1,3) = sqrt(sum((I2model_2(:) - I2data(:)).^2)/length(Save_srk4));
RMSE(2,3) = sqrt(sum((I2Det(:) - I2data(:)).^2)/length(Save_srk4));
RMSE(3,4) = sqrt(sum((I3ave_srk4(:) - I3data(:)).^2)/length(Save_srk4));
RMSE(1,4) = sqrt(sum((I3model_2(:) - I3data(:)).^2)/length(Save_srk4));
RMSE(2,4) = sqrt(sum((I3Det(:) - I3data(:)).^2)/length(Save_srk4));
RMSE(3,5) = sqrt(sum((Rave_srk4(:) - Rdata(:)).^2)/length(Save_srk4));
RMSE(1,5) = sqrt(sum((Rmodel_2(:) - Rdata(:)).^2)/length(Save_srk4));
RMSE(2,5) = sqrt(sum((RDet(:) - Rdata(:)).^2)/length(Save_srk4));
RMSE(3,6) = sqrt(sum((Dave_srk4(:) - Ddata(:)).^2)/length(Save_srk4));
RMSE(1,6) = sqrt(sum((Dmodel_2(:) - Ddata(:)).^2)/length(Save_srk4));
RMSE(2,6) = sqrt(sum((DDet(:) - Ddata(:)).^2)/length(Save_srk4));

Outputdata = cat(1, RDet, I1Det, DDet, I2Det, I3Det, SDet, Rmodel_2, I1model_2, Dmodel_2, I2model_2, I3model_2, Smodel_2,Rave_srk4, I1ave_srk4, Dave_srk4, I2ave_srk4, I3ave_srk4, Save_srk4);


