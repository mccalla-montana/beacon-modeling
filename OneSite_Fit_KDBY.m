clear
clc
close all

%% Load Data

[file,path] = uigetfile('.xlsx');

% Define data ranges
lowtemp = 11;               % Set temperature range based on data input
hightemp = 46; 
f = 8;                      % number of target concentrations
g = hightemp-lowtemp+1;     % number of temperatures 

B_init = 200*10^-9;         % Initial beacon [M]

% Fluorescence, Temperature and Concentration
D = readmatrix(file,'Sheet','NormFluoro');

T = D(:,1); % Temperature
Temperature = T(lowtemp:hightemp,:);

F_NTC = D(:,2);                      % fluorescence of NTC
F_B = D(:,3);                        % fluorescence of 0.2uM beacon samples
F_BT = D(:,4:11);                    % fluorescence of beacon and trigger samples
F_BT_mid = F_BT(lowtemp:hightemp,:); % fluorescence at mid temperature range 30-55C

conc = readmatrix(file,'Sheet','conc'); % trigger concenctration uM

% Initial Dissociation Constants for the chosen temperature range
KD = readmatrix(file,'Sheet','KD');     % dissociation constant inputs [M]
KB = KD(:,1);                           % Beacon dissociation constants 
KB_mid = KB(lowtemp:hightemp,:);
KBY = KD(:,2);                          % Beacon-target dissociation constant
KBY = KBY(lowtemp:hightemp,:);
KYY = KD(:,3);                          % Target dimerization dissociation constant
KYY = KYY(lowtemp:hightemp);

% Weights
w = readmatrix(file,'Sheet','weights');
w = w(lowtemp:hightemp,1:8);

% Target Dimerization Parameters
Y = readmatrix(file,'Sheet','Y');
Y = Y(lowtemp:hightemp,:);
B = readmatrix(file,'Sheet','B');
B = B(lowtemp:hightemp,:);

%% Fluorescent parameters

% [BTT] fluorescent component at 10C,20 uM
beta = F_BT(1,8);

% Beacon fluorescence at 10C (closed beacon)
gamma = F_B(1);

% Random coil fluorescence (open beacon)
delta = F_B(81,1); 

%% Call functions

% Preallocate the optimal KD matrix
KBY_opt = zeros(g,1);        % Optimal KDBY 
minSSR = zeros(g,1);
Yfopt = zeros(g,f);

% Set constraints
A = [];   % linear inequality constraints - none         
b = [];   % NONE
Aeq = []; % linear equality constraints - none
beq=[];   % NONE
nonlcon = []; % NONE
 
optfmin = optimoptions('fmincon','Algorithm','interior-point');

for m = 1:g % loop over temperatures

% Initial upper and lower bounds; Set constraints and constraining regions as needed    

    if m <= 1               
    lb = [50*KBY(m)];   
    ub = [0.01];        
    elseif (m>1) && (m<=15)      
         lb = [KBY_opt(m-1)];  
         ub = [5*KBY_opt(m-1)];  
    else 
        lb = [KBY_opt(m-1)];   
        ub = [1.5*KBY_opt(m-1)];
    end

        % Initial dissociation constants 
        K = [KBY(m)]; % [Molar]  
        k3 = KYY(m);


        fun = @(K) obj_fluoro(K, Y, B,k3, beta, gamma, delta,conc, KB_mid, B_init, F_BT_mid,w,f,m);
        problem = createOptimProblem('fmincon','x0',K,'objective',fun,'lb',lb,'ub',ub,'options',optfmin);
        ms = MultiStart;
        [K,fval]=run(ms,problem,5);


        KBY_opt(m,:) = K(1);
        minSSR(m,1) = fval;

end



%% Calculate optimal fluorescence

%fit_conc = linspace(0,20,100);
for m = 1:g  % Loop over temperatures
    for j = 1:f % Loop over concentrations


        % Free oligo guesses [Molar]
        x = [Y(m,j), B(m,j)];

        optfsolve = optimoptions('fsolve','Algorithm','levenberg-marquardt','Display','iter', ...
            'OptimalityTolerance',1e-10,'FunctionTolerance',1e-60,'StepTolerance',1e-20);

        x = fsolve(@(x)free(x,conc,KYY,KBY_opt,B_init,m,j),x,optfsolve);

            Yfopt(m,j)= x(1);
            Bfopt(m,j) = x(2);

        Fluoro(m,j) = (beta*Yfopt(m,j) + gamma*KBY_opt(m)+delta*KB_mid(m)*KBY_opt(m))/...
        (Yfopt(m,j)+KBY_opt(m)+KBY_opt(m)*KB_mid(m));

    end
end

%% Statistics
for m = 1:g  % Loop over temperatures
    meanY(m) = mean(F_BT_mid(m,:));

    for j = 1:f % Loop over concentrations

         Res(m,j) = F_BT_mid(m,j) - Fluoro(m,j);
         W(m,j) = (F_BT_mid(m,j)-meanY(m));
    end
    
         SSR_opt(m) = sum(Res(m,:).^2); 
         SST(m) = sum(W(m,:).^2); 
         Rsquare(m) = 1 - (SSR_opt(m)/SST(m));
end

%% Calculate fit line 

conc_fit = linspace(0.05*10^-6,20*10^-6,200);

for m = 1:g  % Loop over temperatures
    for j = 1:200 % Loop over fit concentrations
    
        Fluoro_fit(m,j) = (beta*conc_fit(j) + gamma*KBY_opt(m)+delta*KB_mid(m)*KBY_opt(m))/...
        (conc_fit(j)+KBY_opt(m)+KBY_opt(m)*KB_mid(m));

    end
end

%% Export data sheets
filename =  'OneSite_Fit_KDBY.xlsx';
writematrix(F_BT_mid,filename,'Sheet','Data')
writematrix(Fluoro_fit,filename,'Sheet','Fluoro_fit');
writematrix(conc_fit,filename,'Sheet','conc_fit');
writematrix(Yfopt,filename,'Sheet','Yfopt');
writematrix(Bfopt,filename,'Sheet','Bfopt');
writematrix(Fluoro,filename,'Sheet','Fluoro');
writematrix(Temperature,filename,'Sheet','KD_input','Range','A1');
writematrix(KBY,filename,'Sheet','KD_input','Range','B1');
writematrix(Temperature,filename,'Sheet','KD_opt','Range','A1');
writematrix(KBY_opt,filename,'Sheet','KD_opt','Range','B1');
writematrix(Res,filename,'Sheet','Residuals');
writematrix(Rsquare,filename,'Sheet','R2');

%% Plots

figure(1)
tx = tiledlayout(2,2);

nexttile

semilogy(Temperature,KBY_opt)
hold on
plot(Temperature,KBY)
title('Dissociation Constants')
xlabel('Temperature (Celsius)')
ylabel('Molar')
legend('KD_{BY,opt}','KD_{BY}')
hold off

nexttile
hold on
plot(Temperature,Bfopt,'x')
plot(Temperature,B)
xlabel('Temperature (Celsius)')
ylabel('Molar')
legend('Bf_{opt}','Bf_{Nupack}')
hold off

nexttile
hold on
plot(Temperature,Yfopt,'x')
plot(Temperature,Y)
xlabel('Temperature (Celsius)')
ylabel('Molar')
legend('Yf_{opt}','Yf_{Nupack}')
hold off

nexttile
hold on
plot(Temperature,F_BT_mid)
plot(Temperature,Fluoro,'x')
ylim([0,2])
title('Fluorescence fit')
xlabel('Temperature (Celsius)')
ylabel('Fluorescence')
subtitle(['Beta = ', num2str(beta),])
hold off

%%
figure(2)
t = tiledlayout(3,5);

nexttile
q = 1;
hold on 
plot(conc_fit,Fluoro_fit(q,:))
plot(conc,F_BT_mid(q,:),'o')
ylim([0,2])
subtitle(['T = ', num2str(Temperature(q)),'C; R^2 =', num2str(Rsquare(q))]);
hold off

nexttile
q = 3;
hold on 
plot(conc_fit,Fluoro_fit(q,:))
plot(conc,F_BT_mid(q,:),'o')
ylim([0,2])
subtitle(['T = ', num2str(Temperature(q)),'C; R^2 =', num2str(Rsquare(q))]);
hold off

nexttile
q = 5;
hold on 
plot(conc_fit,Fluoro_fit(q,:))
plot(conc,F_BT_mid(q,:),'o')
ylim([0,2])
subtitle(['T = ', num2str(Temperature(q)),'C; R^2 =', num2str(Rsquare(q))]);
hold off

nexttile
q = 7;
hold on 
plot(conc_fit,Fluoro_fit(q,:))
plot(conc,F_BT_mid(q,:),'o')
ylim([0,2])
subtitle(['T = ', num2str(Temperature(q)),'C; R^2 =', num2str(Rsquare(q))]);
hold off

nexttile
q = 9;
hold on 
plot(conc_fit,Fluoro_fit(q,:))
plot(conc,F_BT_mid(q,:),'o')
ylim([0,2])
subtitle(['T = ', num2str(Temperature(q)),'C; R^2 =', num2str(Rsquare(q))]);
hold off

nexttile
q = 11;
hold on 
plot(conc_fit,Fluoro_fit(q,:))
plot(conc,F_BT_mid(q,:),'o')
ylim([0,2])
subtitle(['T = ', num2str(Temperature(q)),'C; R^2 =', num2str(Rsquare(q))]);
hold off

nexttile
q = 13;
hold on 
plot(conc_fit,Fluoro_fit(q,:))
plot(conc,F_BT_mid(q,:),'o')
ylim([0,2])
subtitle(['T = ', num2str(Temperature(q)),'C; R^2 =', num2str(Rsquare(q))]);
hold off

nexttile
q = 15;
hold on 
plot(conc_fit,Fluoro_fit(q,:))
plot(conc,F_BT_mid(q,:),'o')
ylim([0,2])
subtitle(['T = ', num2str(Temperature(q)),'C; R^2 =', num2str(Rsquare(q))]);
hold off

nexttile
q = 17;
hold on 
plot(conc_fit,Fluoro_fit(q,:))
plot(conc,F_BT_mid(q,:),'o')
ylim([0,2])
subtitle(['T = ', num2str(Temperature(q)),'C; R^2 =', num2str(Rsquare(q))]);
hold off

nexttile
q = 20;
hold on 
plot(conc_fit,Fluoro_fit(q,:))
plot(conc,F_BT_mid(q,:),'o')
ylim([0,2])
subtitle(['T = ', num2str(Temperature(q)),'C; R^2 =', num2str(Rsquare(q))]);
hold off

nexttile
q = 23;
hold on 
plot(conc_fit,Fluoro_fit(q,:))
plot(conc,F_BT_mid(q,:),'o')
ylim([0,2])
subtitle(['T = ', num2str(Temperature(q)),'C; R^2 =', num2str(Rsquare(q))]);
hold off

nexttile
q = 26;
hold on 
plot(conc_fit,Fluoro_fit(q,:))
plot(conc,F_BT_mid(q,:),'o')
ylim([0,2])
subtitle(['T = ', num2str(Temperature(q)),'C; R^2 =', num2str(Rsquare(q))]);
hold off

nexttile
q = 30;
hold on 
plot(conc_fit,Fluoro_fit(q,:))
plot(conc,F_BT_mid(q,:),'o')
ylim([0,2])
subtitle(['T = ', num2str(Temperature(q)),'C; R^2 =', num2str(Rsquare(q))]);
hold off

nexttile
q = 33;
hold on 
plot(conc_fit,Fluoro_fit(q,:))
plot(conc,F_BT_mid(q,:),'o')
ylim([0,2])
subtitle(['T = ', num2str(Temperature(q)),'C; R^2 =', num2str(Rsquare(q))]);
hold off

nexttile
q = 36;
hold on 
plot(conc_fit,Fluoro_fit(q,:))
plot(conc,F_BT_mid(q,:),'o')
ylim([0,2])
subtitle(['T = ', num2str(Temperature(q)),'C; R^2 =', num2str(Rsquare(q))]);
hold off


title(t,'Experimental Data and Fit') 
xlabel(t,'[Target] M')
ylabel(t,'Normalized Fluorescence')
t.TileSpacing = 'compact';

%%
figure(3)
hold on
for k = 1:g
   plot(log(Yfopt(k,:)),Res(k,:),'o');
end
hold off


%% Objective function to minimize sum of squares for fluorescent data
function SSR_sum = obj_fluoro(K, Y, B,k3, beta, gamma, delta,conc, KB_mid, B_init, F_BT_mid,w,f,m)

    for j = 1:f
    
        optfsolve = optimoptions('fsolve','Algorithm','levenberg-marquardt','Display','iter', ...
            'OptimalityTolerance',1e-10,'FunctionTolerance',1e-60,'StepTolerance',1e-20);
      
        % x(1) = free target; x(2) = free beacon;  
        x = [Y(m,j), B(m,j)];


        x = fsolve(@(x)model(x,conc,k3,K,B_init,j),x,optfsolve);
        disp(x)

        Yf(m,j)= x(1);                                       
        Bf(m,j) = x(2);

        Fluoro(m,j) = (beta*Yf(m,j) + gamma*K(1)+delta*KB_mid(m)*K(1))/...
        (Yf(m,j)+K(1)+K(1)*KB_mid(m));


     % Sum of Squares  
     SSR = sum(w(m,j).*(F_BT_mid(m,j)-Fluoro(m,j)).^2); 

    end
     
    % Function to minimize (sum of the sums)
    SSR_sum = sum(SSR); 


     function h = model(x,conc,k3,K,B_init,j)


         h(1) = ((2*x(1)^2)/k3) + x(1)*(x(2)/K + 1) - conc(j);
         h(2) = x(2)*(x(1)/K + 1) - B_init;

     end
end


%% Functions

% Function to solve for free target and beacon
function p = free(x,conc,KD_YY,K1_opt,B_init,m,j)

% x(1) = free target; x(2) = free beacon;  

        p(1) = ((2*x(1)^2)/KD_YY(m)) + x(1)*(x(2)/K1_opt(m) + 1) - conc(j);
        p(2) = x(2)*(x(1)/K1_opt(m) + 1) - B_init;

end

