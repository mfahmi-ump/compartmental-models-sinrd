clear model data params options
datafile = readmatrix('sir_data_Malaysia.csv');
pa = readmatrix('Malaysia-SIR-param-estimate.csv');
Ptotal = 3.153*10^7;%malaysia
rng = 3683;
RMSE = double.empty;
data_iter = 0;
for dur_length_iter = 3:90
    for iter_start = 1:(length(datafile)-dur_length_iter-9)
        data_iter = data_iter + 1
        Idata = datafile(9+iter_start:9+dur_length_iter+iter_start,3);
        Rdata = datafile(9+iter_start:9+dur_length_iter+iter_start,2) + datafile(9+iter_start:9+dur_length_iter+iter_start,4);
        Sdata = datafile(9+iter_start:9+dur_length_iter+iter_start,5);
        tdata = datafile(9+iter_start:9+dur_length_iter+iter_start,1);
        t = 1:1:length(Idata);
        tdays = datetime(2020,3,26) + iter_start + caldays(1:length(t));
        Pastparamestimate = pa(1+iter_start:1+dur_length_iter+iter_start,:);
        Szero = Sdata(1);
        Izero = Idata(1);
        Rzero = Rdata(1);
        T = length(t);
        N = T;
        M = 500;
        dt = T/N; R = 1; Dt = R*dt; L = N/R;
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
        fS = @(S,I,R,Beta,r,sigma) -S.*I./Ptotal.*Beta;
        gS = @(S,I,R,Beta,r,sigma) -sigma.*I.*S./Ptotal;
        fI = @(S,I,R,Beta,r,sigma) S./Ptotal.*Beta.*I-r.*I;
        gI = @(S,I,R,Beta,r,sigma) sigma.*I.*S./Ptotal;
        fR = @(S,I,R,Beta,r,sigma) r.*I;
        gR = @(S,I,R,Beta,r,sigma) 0;
        %% SRK4==================================================
        sigma = 0.05;
        %======================================================%SRK-4 Method
        Ssrk4_2 = zeros(M,N);
        Isrk4_2 = zeros(M,N);
        Rsrk4_2 = zeros(M,N);
        Stemp = Szero*ones(M,1);
        Itemp = Izero*ones(M,1);
        Rtemp = Rzero*ones(M,1);
        Smodel = ones(3,N);
        Imodel = ones(3,N);
        Rmodel = ones(3,N);
        SdW1 = sqrt(dt)*randn(M,N);       %Numerical Brownian increments
        SdW2 = sqrt(dt)*randn(M,N);       %Numerical Brownian increments
        IdW1 = sqrt(dt)*randn(M,N);       %Numerical Brownian increments
        IdW2 = sqrt(dt)*randn(M,N);       %Numerical Brownian increments
        RdW1 = sqrt(dt)*randn(M,N);       %Numerical Brownian increments
        RdW2 = sqrt(dt)*randn(M,N);       %Numerical Brownian increments

        for j = 1:L
            beta = Pastparamestimate(j,2);
            r = Pastparamestimate(j,3);

            SWinc1 = sum(SdW1(:,R*(j-1)+1:R*j),2);
            SWinc2 = sum(SdW2(:,R*(j-1)+1:R*j),2);
            IWinc1 = sum(IdW1(:,R*(j-1)+1:R*j),2);
            IWinc2 = sum(IdW2(:,R*(j-1)+1:R*j),2);
            RWinc1 = sum(RdW1(:,R*(j-1)+1:R*j),2);
            RWinc2 = sum(RdW2(:,R*(j-1)+1:R*j),2);

            SJ10 = (0.5*(Dt^(1/2))*(SWinc1+(1/sqrt(3)).*SWinc2));
            IJ10 = (0.5*(Dt^(1/2))*(IWinc1+(1/sqrt(3)).*IWinc2));
            RJ10 = (0.5*(Dt^(1/2))*(RWinc1+(1/sqrt(3)).*RWinc2));

            S1 = Stemp;
            I1 = Itemp;
            R1 = Rtemp;
            %Refer Equation 4.56
            f1_S = fS(S1, I1, R1, beta, r ,sigma);
            g1_S = gS(S1, I1, R1, beta, r ,sigma);
            f1_I = fI(S1, I1, R1, beta, r ,sigma);
            g1_I = gI(S1, I1, R1, beta, r ,sigma);
            f1_R = fR(S1, I1, R1, beta, r ,sigma);
            g1_R = gR(S1, I1, R1, beta, r ,sigma);

            S2 = Stemp + Dt*(A(2)*f1_S) + (B_1(2).*SWinc1 + B_2(2).*SJ10).*g1_S;
            I2 = Itemp + Dt*(A(2)*f1_I) + (B_1(2).*IWinc1 + B_2(2).*IJ10).*g1_I;
            R2 = Rtemp + Dt*(A(2)*f1_R) + (B_1(2).*RWinc1 + B_2(2).*RJ10).*g1_R;

            %Refer Equation 4.56
            f2_S = fS(S2, I2, R2, beta, r ,sigma);
            g2_S = gS(S2, I2, R2, beta, r ,sigma);
            f2_I = fI(S2, I2, R2, beta, r ,sigma);
            g2_I = gI(S2, I2, R2, beta, r ,sigma);
            f2_R = fR(S2, I2, R2, beta, r ,sigma);
            g2_R = gR(S2, I2, R2, beta, r ,sigma);

            S3 = Stemp + Dt*(A(7)*f2_S) + (B_1(3).*SWinc1 + B_2(3).*SJ10).*g2_S + B_1(7)*SJ10.*g2_S;
            I3 = Itemp + Dt*(A(7)*f2_I) + (B_1(3).*IWinc1 + B_2(3).*IJ10).*g2_I + B_1(7)*IJ10.*g2_I;
            R3 = Rtemp + Dt*(A(7)*f2_R) + (B_1(3).*RWinc1 + B_2(3).*RJ10).*g2_R + B_1(7)*RJ10.*g2_R;

            %Refer Equation 4.56
            f3_S = fS(S3, I3, R3, beta, r ,sigma);
            g3_S = gS(S3, I3, R3, beta, r ,sigma);
            f3_I = fI(S3, I3, R3, beta, r ,sigma);
            g3_I = gI(S3, I3, R3, beta, r ,sigma);
            f3_R = fR(S3, I3, R3, beta, r ,sigma);
            g3_R = gR(S3, I3, R3, beta, r ,sigma);

            S4 = Stemp + Dt*(f3_S) + (B_1(4).*SWinc1 + B_2(4).*SJ10).*g3_S + (B_1(8).*SWinc1.*g3_S + B_1(12).*SWinc1.*g3_S);
            I4 = Itemp + Dt*(f3_I) + (B_1(4).*IWinc1 + B_2(4).*IJ10).*g3_I + (B_1(8).*IWinc1.*g3_I + B_1(12).*IWinc1.*g3_I);
            R4 = Rtemp + Dt*(f3_R) + (B_1(4).*RWinc1 + B_2(4).*RJ10).*g3_R + (B_1(8).*RWinc1.*g3_R + B_1(12).*RWinc1.*g3_R);
            %Refer Equation 4.56

            f4_S = fS(S4, I4, R4, beta, r ,sigma);
            g4_S = gS(S4, I4, R4, beta, r ,sigma);
            f4_I = fI(S4, I4, R4, beta, r ,sigma);
            g4_I = gI(S4, I4, R4, beta, r ,sigma);
            f4_R = fR(S4, I4, R4, beta, r ,sigma);
            g4_R = gR(S4, I4, R4, beta, r ,sigma);

            Stemp = Stemp + Dt.*((1/6)*f1_S + (1/3).*f2_S + (1/3).*f3_S + (1/6).*f4_S) + (G_1(1).*SWinc1 + G_2(1).*SJ10).*g1_S + (G_1(2).*SWinc1 + G_2(2).*SJ10).*g2_S + (G_1(3).*SWinc1 + G_2(3).*SJ10).*g3_S + (G_1(4).*SWinc1 + G_2(4).*SJ10).*g4_S;
            Itemp = Itemp + Dt.*((1/6)*f1_I + (1/3).*f2_I + (1/3).*f3_I + (1/6).*f4_I) + (G_1(1).*IWinc1 + G_2(1).*IJ10).*g1_I + (G_1(2).*IWinc1 + G_2(2).*IJ10).*g2_I + (G_1(3).*IWinc1 + G_2(3).*IJ10).*g3_I + (G_1(4).*IWinc1 + G_2(4).*IJ10).*g4_I;
            Rtemp = Rtemp + Dt.*((1/6)*f1_R + (1/3).*f2_R + (1/3).*f3_R + (1/6).*f4_R) + (G_1(1).*RWinc1 + G_2(1).*RJ10).*g1_R + (G_1(2).*RWinc1 + G_2(2).*RJ10).*g2_R + (G_1(3).*RWinc1 + G_2(3).*RJ10).*g3_R + (G_1(4).*RWinc1 + G_2(4).*RJ10).*g4_R;

            %Refer Equation 4.56
            Ssrk4_2(:,j) = Stemp.';
            Isrk4_2(:,j) = Itemp.';
            Rsrk4_2(:,j) = Rtemp.';
        end
        for i = 1:L
           Ssrk_3 = Ssrk4_2(:,i);
           Isrk_3 = Isrk4_2(:,i);
           Rsrk_3 = Rsrk4_2(:,i);
           Smodel(1,i) = interp1(linspace(1/numel(Ssrk_3),1,numel(Ssrk_3)), sort(Ssrk_3), 0.25);
           Imodel(1,i) = interp1(linspace(1/numel(Isrk_3),1,numel(Isrk_3)), sort(Isrk_3), 0.25);
           Rmodel(1,i) = interp1(linspace(1/numel(Rsrk_3),1,numel(Rsrk_3)), sort(Rsrk_3), 0.25);
           Smodel(2,i) = interp1(linspace(1/numel(Ssrk_3),1,numel(Ssrk_3)), sort(Ssrk_3), 0.50);
           Imodel(2,i) = interp1(linspace(1/numel(Isrk_3),1,numel(Isrk_3)), sort(Isrk_3), 0.50);
           Rmodel(2,i) = interp1(linspace(1/numel(Rsrk_3),1,numel(Rsrk_3)), sort(Rsrk_3), 0.50);
           Smodel(3,i) = interp1(linspace(1/numel(Ssrk_3),1,numel(Ssrk_3)), sort(Ssrk_3), 0.75);
           Imodel(3,i) = interp1(linspace(1/numel(Isrk_3),1,numel(Isrk_3)), sort(Isrk_3), 0.75);
           Rmodel(3,i) = interp1(linspace(1/numel(Rsrk_3),1,numel(Rsrk_3)), sort(Rsrk_3), 0.75);
        end
        Save_srk4 = mean(Ssrk4_2,1);
        Iave_srk4 = mean(Isrk4_2,1);
        Rave_srk4 = mean(Rsrk4_2,1);
        %% RK4==================================================
        sigma = 0;
        %======================================================%SRK-4 Method
        Srk4_2 = zeros(1,N);
        Irk4_2 = zeros(1,N);
        Rrk4_2 = zeros(1,N);
        Stemp = Szero*ones(1,1);
        Itemp = Izero*ones(1,1);
        Rtemp = Rzero*ones(1,1);

        for j = 1:L
            beta = Pastparamestimate(j,2);
            r = Pastparamestimate(j,3);

            S1 = Stemp;
            I1 = Itemp;
            R1 = Rtemp;
            %Refer Equation 4.56
            f1_S = fS(S1, I1, R1, beta, r ,sigma);
            f1_I = fI(S1, I1, R1, beta, r ,sigma);
            f1_R = fR(S1, I1, R1, beta, r ,sigma);

            S2 = Stemp + Dt*(A(2)*f1_S); 
            I2 = Itemp + Dt*(A(2)*f1_I);
            R2 = Rtemp + Dt*(A(2)*f1_R); 

            %Refer Equation 4.56
            f2_S = fS(S2, I2, R2, beta, r ,sigma);
            f2_I = fI(S2, I2, R2, beta, r ,sigma);
            f2_R = fR(S2, I2, R2, beta, r ,sigma);

            S3 = Stemp + Dt*(A(7)*f2_S);
            I3 = Itemp + Dt*(A(7)*f2_I);
            R3 = Rtemp + Dt*(A(7)*f2_R);

            %Refer Equation 4.56
            f3_S = fS(S3, I3, R3, beta, r ,sigma);
            f3_I = fI(S3, I3, R3, beta, r ,sigma);
            f3_R = fR(S3, I3, R3, beta, r ,sigma);

            S4 = Stemp + Dt*(f3_S);
            I4 = Itemp + Dt*(f3_I);
            R4 = Rtemp + Dt*(f3_R);
            %Refer Equation 4.56

            f4_S = fS(S4, I4, R4, beta, r ,sigma);
            g4_S = gS(S4, I4, R4, beta, r ,sigma);
            f4_I = fI(S4, I4, R4, beta, r ,sigma);
            g4_I = gI(S4, I4, R4, beta, r ,sigma);
            f4_R = fR(S4, I4, R4, beta, r ,sigma);
            g4_R = gR(S4, I4, R4, beta, r ,sigma);

            Stemp = Stemp + Dt.*((1/6)*f1_S + (1/3).*f2_S + (1/3).*f3_S + (1/6).*f4_S);
            Itemp = Itemp + Dt.*((1/6)*f1_I + (1/3).*f2_I + (1/3).*f3_I + (1/6).*f4_I);
            Rtemp = Rtemp + Dt.*((1/6)*f1_R + (1/3).*f2_R + (1/3).*f3_R + (1/6).*f4_R);

            %Refer Equation 4.56
            Srk4_2(:,j) = Stemp.';
            Irk4_2(:,j) = Itemp.';
            Rrk4_2(:,j) = Rtemp.';
        end
        % RMSE_srk4x = rmse(Xsrk4_2,Datax);
        % RMSE_srk4y = rmse(Ysrk4_2,Datay);
        Save_rk4 = mean(Srk4_2,1);
        Iave_rk4 = mean(Irk4_2,1);
        Rave_rk4 = mean(Rrk4_2,1);
        % Smodel = Save_srk4;
        Smodel_2 = Save_rk4;
        Imodel_2 = Iave_rk4;
        Rmodel_2 = Rave_rk4;
        %% Adam_Bashford_Predictor_Corrector=====================
        Sabpc_2 = Szero*ones(1,N+1);
        Iabpc_2 = Izero*ones(1,N+1);
        Rabpc_2 = Rzero*ones(1,N+1);

        pcS1 = 0; pcS2 = 0; pcS3 = 0;
        pcI1 = 0; pcI2 = 0; pcI3 = 0;
        pcR1 = 0; pcR2 = 0; pcR3 = 0;

        for j = 1:L
            beta = Pastparamestimate(j,2);
            r = Pastparamestimate(j,3);

            S1 = Sabpc_2(j);
            I1 = Iabpc_2(j);
            R1 = Rabpc_2(j);

            f1_S = fS(S1,I1,R1,beta,r ,sigma);
            f1_I = fI(S1,I1,R1,beta,r ,sigma);
            f1_R = fR(S1,I1,R1,beta,r ,sigma);

            pcS0 = f1_S;
            pcI0 = f1_I;
            pcR0 = f1_R;

            S2 = Sabpc_2(j) + dt*(0.5*f1_S);
            I2 = Iabpc_2(j) + dt*(0.5*f1_I);
            R2 = Rabpc_2(j) + dt*(0.5*f1_R);

            f2_S = fS(S2,I2,R2,beta,r ,sigma);
            f2_I = fI(S2,I2,R2,beta,r ,sigma);
            f2_R = fR(S2,I2,R2,beta,r ,sigma);

            S3 = Sabpc_2(j) + dt*(0.5*f2_S);
            I3 = Iabpc_2(j) + dt*(0.5*f2_I);
            R3 = Rabpc_2(j) + dt*(0.5*f2_R);

            f3_S = fS(S3,I3,R3,beta,r ,sigma);
            f3_I = fI(S3,I3,R3,beta,r ,sigma);
            f3_R = fR(S3,I3,R3,beta,r ,sigma);

            S4 = Sabpc_2(j) + dt*f3_S;
            I4 = Iabpc_2(j) + dt*f3_I;
            R4 = Rabpc_2(j) + dt*f3_R;

            f4_S = fS(S4,I4,R4,beta,r ,sigma);
            f4_I = fI(S4,I4,R4,beta,r ,sigma);
            f4_R = fR(S4,I4,R4,beta,r ,sigma);

            Sabpc_2(j+1) = Sabpc_2(j) + dt*((1/6)*f1_S + (1/3)*f2_S + (1/3)*f3_S + (1/6)*f4_S);
            Iabpc_2(j+1) = Iabpc_2(j) + dt*((1/6)*f1_I + (1/3)*f2_I + (1/3)*f3_I + (1/6)*f4_I);
            Rabpc_2(j+1) = Rabpc_2(j) + dt*((1/6)*f1_R + (1/3)*f2_R + (1/3)*f3_R + (1/6)*f4_R);

            pcS1 = pcS0; pcS2 = pcS1; pcS3 = pcS2;
            pcI1 = pcI0; pcI2 = pcI1; pcI3 = pcI2;
            pcR1 = pcR0; pcR2 = pcR1; pcR3 = pcR2;
        end

        for j = 1:L
            beta = Pastparamestimate(j,2);
            r = Pastparamestimate(j,3);

            S1 = Sabpc_2(j);
            I1 = Iabpc_2(j);
            R1 = Rabpc_2(j);

            f1_S = fS(S1,I1,R1,beta,r);
            f1_I = fI(S1,I1,R1,beta,r);
            f1_R = fR(S1,I1,R1,beta,r);

            pcS0 = f1_S;
            pcI0 = f1_I;
            pcR0 = f1_R;

            wS = Sabpc_2(j) + dt*( 55.0 * pcS0 - 59.0 * pcS1 + 37.0 * pcS2 - 9.0 * pcS3 ) / 24.0;
            wI = Iabpc_2(j) + dt*( 55.0 * pcI0 - 59.0 * pcI1 + 37.0 * pcI2 - 9.0 * pcI3 ) / 24.0;
            wR = Rabpc_2(j) + dt*( 55.0 * pcR0 - 59.0 * pcR1 + 37.0 * pcR2 - 9.0 * pcR3 ) / 24.0;

            fw_S = fS(wS,wI,wR,beta,r ,sigma);
            fw_I = fI(wS,wI,wR,beta,r ,sigma);
            fw_R = fR(wS,wI,wR,beta,r ,sigma);

            Sabpc_2(j+1) = Sabpc_2(j) + dt*( 9.0 * fw_S + 19.0 * pcS0 - 5.0 * pcS1 + pcS2 ) / 24.0;
            Iabpc_2(j+1) = Iabpc_2(j) + dt*( 9.0 * fw_I + 19.0 * pcI0 - 5.0 * pcI1 + pcI2 ) / 24.0;
            Rabpc_2(j+1) = Rabpc_2(j) + dt*( 9.0 * fw_R + 19.0 * pcR0 - 5.0 * pcR1 + pcR2 ) / 24.0;

            pcS1 = pcS0; pcS2 = pcS1; pcS3 = pcS2;
            pcI1 = pcI0; pcI2 = pcI1; pcI3 = pcI2;
            pcR1 = pcR0; pcR2 = pcR1; pcR3 = pcR2;
        end
        % RMSE_srk4x = rmse(Xsrk4_2,Datax);
        % RMSE_srk4y = rmse(Ysrk4_2,Datay);
        Save_abpc = mean(Sabpc_2,1);
        Iave_abpc = mean(Iabpc_2,1);
        Rave_abpc = mean(Rabpc_2,1);
        % Smodel = Save_srk4;
        SDet = Save_abpc(2:end);
        IDet = Iave_abpc(2:end);
        RDet = Rave_abpc(2:end);
        % %% Plot
        % figure('Position', [1, 3, 900, 300])
        % tlayout = tiledlayout(1,3);
        % tlayout.TileSpacing = 'compact';
        % tlayout.Padding = 'compact';
        % 
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
        %     xlabel('time (date)'); ylabel('Susceptible Population')    
        %     
        %     ax2 = nexttile; 
        %     plot(tdays,Idata,"-s",'Color',"#2C00D2",'MarkerIndices',1:20:length(tdays)); hold on;
        %     plot(tdays,Imodel_2,':','Color',"#273DE3"); hold on;
        %     plot(tdays,IDet,'--', 'Color',"#4BA3E5")
        %     plot(tdays,Imodel(2,:),'-','Color',"#77DFE9"); hold on;
        %     cax2 = [39/255,61/255,227/255];
        %     fill([tdays, fliplr(tdays)], [Imodel(1,:), fliplr(Imodel(3,:))], cax2, 'LineStyle', "None",'FaceAlpha',0.2); hold on;
        %     title(ax2,'Infected Population')
        %     ax2.FontSize = 8;
        %     legend('Infected Data','RK4','Adam-Moulton PC','SRK4','Location','best')
        %     hold off;
        %     xlabel('time (date)'); ylabel('Infected Population')    
        %     
        %     ax3 = nexttile;
        %     plot(tdays,Rdata,'-s','Color',"#008080",'MarkerIndices',1:20:length(tdays)); hold on;
        %     plot(tdays,Rmodel_2,':','Color',"#556B2F");hold on;
        %     plot(tdays,RDet,'--', 'Color',"#3CB371")
        %     plot(tdays,Rmodel(2,:),'-','Color',"#00FF00");hold on;
        %     title(ax3,'Recovered Population, R')
        %     cax3 = [70/255,192/255,30/255];
        %     fill([tdays, fliplr(tdays)], [Rmodel(1,:), fliplr(Rmodel(3,:))], cax3, 'LineStyle', "None",'FaceAlpha',0.2); hold on;
        %     ax3.FontSize = 8;
        %     legend('Recovered Data','RK4','Adam-Moulton PC','SRK4','Location','best')
        %     hold off;
        %     xlabel('time (date)'); ylabel('Recovered Population')    
        %     
        %     print -depsc sirMalaysia.eps
        %% RMSE
        RMSE(data_iter,1) = dur_length_iter;
        RMSE(data_iter,2) = iter_start;
        
        RMSE(data_iter,3) = sqrt(sum((Smodel_2(:) - Sdata(:)).^2)/length(Save_srk4));
        RMSE(data_iter,4) = sqrt(sum((SDet(:) - Sdata(:)).^2)/length(Save_srk4));
        RMSE(data_iter,5) = sqrt(sum((Save_srk4(:) - Sdata(:)).^2)/length(Save_srk4));

        RMSE(data_iter,6) = sqrt(sum((Imodel_2(:) - Idata(:)).^2)/length(Save_srk4));
        RMSE(data_iter,7) = sqrt(sum((IDet(:) - Idata(:)).^2)/length(Save_srk4));
        RMSE(data_iter,8) = sqrt(sum((Iave_srk4(:) - Idata(:)).^2)/length(Save_srk4));

        RMSE(data_iter,9) = sqrt(sum((Rmodel_2(:) - Rdata(:)).^2)/length(Save_srk4));
        RMSE(data_iter,10) = sqrt(sum((RDet(:) - Rdata(:)).^2)/length(Save_srk4));
        RMSE(data_iter,11) = sqrt(sum((Rave_srk4(:) - Rdata(:)).^2)/length(Save_srk4));
    end
end