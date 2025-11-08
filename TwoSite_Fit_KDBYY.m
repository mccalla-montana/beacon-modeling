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
F_BT_mid = F_BT(lowtemp:hightemp,:); % fluorescence at temperature range

conc = readmatrix(file,'Sheet','conc'); % trigger concentration uM

% Initial Dissociation Constants
KD = readmatrix(file,'Sheet','KD');     % dissociation constant inputs [M]
KB = KD(:,1);                           % Beacon dissociation constants 
KB_mid = KB(lowtemp:hightemp,:);
KBY = KD(:,2);                          % Beacon-target dissociation constant
KBY = KBY(lowtemp:hightemp,:);
KBYY = KD(:,3);                         % Beacon-target-target dissociation constant
KBYY = KBYY(lowtemp:hightemp,:);
KYY = KD(:,4);                          % Target dimerization dissociation constant
KYY = KYY(lowtemp:hightemp);
                            

% Weights
w = readmatrix(file,'Sheet','weights');
w = w(lowtemp:hightemp,:); 

% Target Dimerization Parameters
Y = readmatrix(file,'Sheet','Y');  % Free target [M]
Y = Y(lowtemp:hightemp,:);     
B = readmatrix(file,'Sheet','B');  % Free beacon [M]
B = B(lowtemp:hightemp,:); % 



%% Fluorescent parameters

% Fluorescence at 10C, 20uM
alpha = F_BT(1,8); % Two targets bound 
beta = 0.955; % Estimate for one target bound

% Beacon fluorescence at 10C (closed beacon)
gamma = F_B(1);

% Random coil fluorescence at 80C (open beacon)
delta = F_B(81,1); 

%% Call functions

% Preallocate the optimal KD matrix
KBYY_opt = zeros(g,1); %KD,BYY
minSSR = zeros(g,1);
Yfopt = zeros(g,f);
Bfopt = zeros(g,f);

% Set constraints
A = [];   % linear inequality constraints - none         
b = [];   % NONE
Aeq = []; % linear equality constraints - none
beq=[];   % NONE

nonlcon = []; % NONE
 
optfmin = optimoptions('fmincon','Algorithm','interior-point');

for m = 1:g % loop over temperatures
 % Initial upper and lower bounds; Set constraints and constraining regions as needed
    if (m == 1)                   
    lb = [0];                
    ub = [1*KBYY(m)];     
    elseif (m>1) && (m<=15)
        lb = [KBYY_opt(m-1)];   
        ub = [1.3*KBYY_opt(m-1)]; 
    else 
        lb = [KBYY_opt(m-1)];   
        ub = [1.7*KBYY_opt(m-1)]; 
    end

        % Initial dissociation constants 
        K = [KBYY(m)]; % [Molar] 
        k2 = KBY(m);
        k3 = KYY(m);
   
        fun = @(K) obj_fluoro(K, Y, B,k2,k3,alpha, beta, gamma, delta,conc, KB_mid, B_init, F_BT_mid,w,f,m);
        problem = createOptimProblem('fmincon','x0',K,'objective',fun,'lb',lb,'ub',ub,'options',optfmin); 
        ms = MultiStart;
        [K,fval]=run(ms,problem,5);


KBYY_opt(m) = K(1);
minSSR(m,1) = fval;

end

%% Calculate optimal fluorescence


for m = 1:g     % Loop over temperatures
    for j = 1:f % Loop over concentrations


        % Free oligo Nupack estimates [Molar]
        x = [Y(m,j), B(m,j)];

        optfsolve = optimoptions('fsolve','Algorithm','levenberg-marquardt','Display','iter', ...
            'OptimalityTolerance',1e-6,'FunctionTolerance',1e-60,'StepTolerance',1e-20);

        x = fsolve(@(x)free(x,conc,k2,k3,KBYY_opt,B_init,m,j),x,optfsolve);

            Yfopt(m,j)= x(1);
            Bfopt(m,j) = x(2);


         Fluoro(m,j) = (alpha*Yfopt(m,j) + beta*KBYY_opt(m))/(Yfopt(m,j) +...
             KBYY_opt(m) + (1/Yfopt(m,j))*(KBYY_opt(m)*k2 +...
             KBYY_opt(m)*k2*KB_mid(m))) +...
             (gamma+delta*KB_mid(m))/(1+KB_mid(m) +... 
             (Yfopt(m,j)/k2)*(1+(Yfopt(m,j)/KBYY_opt(m))));

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

for m = 1:g         % Loop over temperatures
    for j = 1:200   % Loop over fit concentrations
    
        Fluoro_fit(m,j) = (alpha*conc_fit(j) + beta* KBYY_opt(m))/(conc_fit(j) +...
             KBYY_opt(m) + (1/conc_fit(j))*(KBYY_opt(m)*k2 +...
             KBYY_opt(m)*k2*KB_mid(m))) +...
             (gamma+delta*KB_mid(m))/(1+KB_mid(m) +... 
             (conc_fit(j)/k2)*(1+(conc_fit(j)/KBYY_opt(m))));

    end
end

%% Calculate nH

for i = 1:g
    nH(i) = 2/(1+sqrt(KBYY_opt(i)/KBY(i)));
end


%% Export data sheets

filename =  'TwoSite_Fit_KDBYY.xlsx';
writematrix(F_BT_mid,filename,'Sheet','Data');
writematrix(Fluoro_fit,filename,'Sheet','Fluoro_fit');  
writematrix(conc_fit,filename,'Sheet','conc_fit');
writematrix(Yfopt,filename,'Sheet','Yfopt');
writematrix(Bfopt,filename,'Sheet','Bfopt');
writematrix(Fluoro,filename,'Sheet','Fluoro');
writematrix(Temperature,filename,'Sheet','KD_input','Range','A1');
writematrix(KBY,filename,'Sheet','KD_input','Range','B1');
writematrix(KBYY,filename,'Sheet','KD_input','Range','C1');
writematrix(Temperature,filename,'Sheet','KD_opt','Range','A1');
writematrix(KBYY_opt,filename,'Sheet','KD_opt','Range','B1');
writematrix(Res,filename,'Sheet','Residuals');
writematrix(Rsquare,filename,'Sheet','R2');

%% Plots

figure(1)
tx = tiledlayout(2,2);

nexttile
semilogy(Temperature,KBYY_opt)
hold on
plot(Temperature,KBYY)
plot(Temperature,KBY)
title('Dissociation Constants')
xlabel('Temperature (Celsius)')
ylabel('Molar')
legend('KD_{BYY,opt}','KD_{BYY}','KD_{BY}')
hold off

nexttile
hold on
plot(Temperature,Bfopt,'x')
plot(Temperature,B)
xlabel('Temperature (Celsius)')
ylabel('Molar')
legend('B_{opt}','B0')
hold off

nexttile
hold on
plot(Temperature,Yfopt,'x')
plot(Temperature,Y)
xlabel('Temperature (Celsius)')
ylabel('Molar')
legend('Yf_{opt}','Y0')
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
function SSR_sum = obj_fluoro(K,Y,B,k2,k3,alpha, beta, gamma, delta,conc, KB_mid, B_init, F_BT_mid,w,f,m)

for j = 1:f
    
        optfsolve = optimoptions('fsolve','Algorithm','levenberg-marquardt','Display','iter', ...
            'OptimalityTolerance',1e-6,'FunctionTolerance',1e-60,'StepTolerance',1e-20);
      
        % x(1) = free target; x(2) = free beacon;  
        x = [Y(m,j), B(m,j)];

        x = fsolve(@(x)model(x,conc,k2,k3,K,B_init,j),x,optfsolve);
        disp(x)

        Yf(m,j)= x(1);                                       
        Bf(m,j) = x(2);
 

        Fluoro(m,j) = (alpha*Yf(m,j) + beta*K(1))/(Yf(m,j) +...
             K(1) + (1/Yf(m,j))*(K(1)*k2 +...
             K(1)*k2*KB_mid(m))) +...
             (gamma+delta*KB_mid(m))/(1+KB_mid(m) +... 
             (Yf(m,j)/k2)*(1+(Yf(m,j)/K(1))));


% Sum of Squares  
     SSR = sum(w(m,j).*(F_BT_mid(m,j)-Fluoro(m,j)).^2); 

end
     
% Function to minimize (sum of the sums)
    SSR_sum = sum(SSR); 


     function h = model(x,conc,k2,k3,K,B_init,j)

        h(1) = (2*(x(1)^2)/k3) + x(1) + (x(1)*x(2)/k2) + (2*x(2)*x(1)^2)/(k2*K(1))-conc(j);
        h(2) = x(2)*((x(1)^2)/(k2*K(1))+(x(1)/k2)+1)-B_init;

     end
end


%% Functions

% Function to solve for free target and beacon
function p = free(x,conc,k2,k3,K1_opt,B_init,m,j)

% x(1) = free target; x(2) = free beacon;  

        p(1) = (2*(x(1)^2)/k3) + x(1) + ((x(1)*x(2))/k2) +...
            (2*x(2)*(x(1)^2))/(k2*K1_opt(m))-conc(j);
        p(2) = x(2)*((x(1)^2)/(k2*K1_opt(m))+(x(1)/k2)+1)-B_init;

end

