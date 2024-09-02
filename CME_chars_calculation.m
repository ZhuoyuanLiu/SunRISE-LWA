
% 4/23：
% % BDW = 0.37;
% 1st harmonics:
% % BDW = [(63-53)/53 (62-52.5)/52.5 (63-55)/55 (62-53)/53 (58-52)/52 (57-53)/53 (57-52)/52]; % Of the first harmonic band
% BDW = [(10)/30 (10)/29 (10)/28 (10)/27.5 (10)/27 (10)/26.5 (10)/26]; % Of the fundamental band
% fl = [53 52.5 55 53 52 53 52];
% %f0 = 80; % MHz starting frequency of the 1st harmonic
% f0 = 40; % MHz starting frequency of the fundamental
% X = (BDW+1).^2;
% Ma = sqrt( X.*(X+5) ./ ( 2.*(4-X) ) );
% a = 6.13;  % fitted empirical index of density variation over the heliocentric distance,
% r = [2.03 2.03 2.03 2.04 2.04 2.04 2.05]*695700;       % (km) Discard the last timepoint's r data because the drift rate cannot be obtained there 
% drift_rate = [-0.038 -0.015 -0.02 -0.032 -0.027 -0.047 -0.065];
% Vs = -2.*r/a*(1/f0).*drift_rate;   % Shock speed
% Va = Vs./Ma;
% B = 5.1*10^(-5).*Va.*fl;    % G


% Fundamental (all lower band):
BDW = (10)/35; % Of the fundamental band
fl = 30;
%f0 = 80; % MHz starting frequency of the 1st harmonic
f0 = 40; % MHz starting frequency of the fundamental
X = (BDW+1).^2;
Ma = sqrt( X.*(X+5) ./ ( 2.*(4-X) ) );
a = 6.13;  % fitted empirical index of density variation over the heliocentric distance,
r = 1.77*695700;       % (km) Discard the last timepoint's r data because the drift rate cannot be obtained there 
%drift_rate = -0.09;     % Initial
%drift_rate = -0.047;          %(mean)

% Include all the points spanning a band
drift_rate = -[0.09436834 0.06763064 0.04900531 0.0283105  0.03871522 0.03145611 0.03454005 0.03235486 0.01098104 0.01520756 0.0026259  0.01078495 0.00138522];
%drift_rate = -[0.04790007 0.02903641 0.05238406 0.05347539 0.02106828 0.00603957...
%  0.01980188 0.02648936 0.03960376 0.02322913 0.04793313 0.04455423...
%  0.01279571 0.04553647 0.00702276 0.04273283 0.00274526 0.00592115...
%  0.00479331 0.03568839 0.03568839 0.01161457 0.01957269 0.04361914...
%  0.03682667];
drift_rate = mean(drift_rate);
Vs = -2.*r/a*(1/f0).*drift_rate;   % Shock speed
Va = Vs./Ma;
B = 5.1*10^(-5).*Va.*fl;    % G


%% 2024/4/22
% Fundamental:
% BDW = [(18)/52 (18)/50 (18)/45 18/45 (18)/44 (18)/43 (18)/42]; % Of the fundamental band
% fl = [52 50 45 45 44 43 42];
% %f0 = 80; % MHz starting frequency of the 1st harmonic
% f0 = 40; % MHz starting frequency of the fundamental
% X = (BDW+1).^2;
% Ma = sqrt( X.*(X+5) ./ ( 2.*(4-X) ) );
% a = 6.13;  % fitted empirical index of density variation over the heliocentric distance,
% r = [1.565 1.57 1.573 1.577 1.582 1.587 1.592]*695700;       % (km) Discard the last timepoint's r data because the drift rate cannot be obtained there 
% drift_rate = [-0.113 -0.11 -0.1 -0.119 -0.14 -0.127 -0.15];
% Vs = -2.*r/a*(1/f0).*drift_rate;   % Shock speed
% Va = Vs./Ma;
% B = 5.1*10^(-5).*Va.*fl;    % G

fl = 45;
BDW = 15./45;               % Of the fundamental band
%f0 = 80;                   % MHz starting frequency of the 1st harmonic
f0 = 68;                    % MHz starting frequency of the fundamental
X = (BDW+1).^2;
Ma = sqrt( X.*(X+5) ./ ( 2.*(4-X) ) );
a = 6.13;                   % Fitted empirical index of density variation over the heliocentric distance,
r = 1.5*695700;             % (km) Discard the last timepoint's r data because the drift rate cannot be obtained there 
% drift_rate = -0.9;          % Onset drift rate
% drift_rate = -[0.9 0.4 0.28 0.29 0.35 0.24 0.21 0.15 0.19 0.2 0.17];
% drift_rate = mean(drift_rate);% Average drift rate
drift_rate = -0.07657849085086972;
Vs = -2.*r/a*(1/f0).*drift_rate;   % Shock speed
Va = Vs./Ma;
B = 5.1*10^(-5).*Va.*fl;    % G
list = [Vs;Va;Ma;B];

%% 2024/5/3
% Callisto fundamental I:
% clear;
% fl = 35;
% BDW = 10./fl; % Of the fundamental band
% %f0 = 64; % MHz starting frequency of the 1st harmonic
% f0 = 37; % MHz starting frequency of the fundamental
% X = (BDW+1).^2;
% Ma = sqrt( X.*(X+5) ./ ( 2.*(4-X) ) );
% a = 6.13;  % Fitted empirical index of density variation over the heliocentric distance,
% r = 1.83*695700;       % (km) Discard the last timepoint's r data because the drift rate cannot be obtained there 
% % drift_rate = -0.62;    % Initial
% drift_rate = -0.09334678235076059; % Mean
% Vs = -2.*r/a*(1/f0).*drift_rate;   % Shock speed
% Va = Vs./Ma;
% B = 5.1*10^(-5).*Va.*fl;    % G
% list = [Vs;Va;Ma;B];

% Callisto fundamental II:
% clear;
% fl = 33;
% BDW = 10./fl; % Of the fundamental band
% %f0 = 64; % MHz starting frequency of the 1st harmonic
% f0 = 50; % MHz starting frequency of the fundamental
% X = (BDW+1).^2;
% Ma = sqrt( X.*(X+5) ./ ( 2.*(4-X) ) );
% a = 6.13;  % Fitted empirical index of density variation over the heliocentric distance,
% r = 1.69*695700;       % (km) Discard the last timepoint's r data because the drift rate cannot be obtained there 
% % drift_rate = -0.62;    % Initial
% drift_rate = -0.042; % Mean
% Vs = -2.*r/a*(1/f0).*drift_rate;   % Shock speed
% Va = Vs./Ma;
% B = 5.1*10^(-5).*Va.*fl;    % G
% list = [Vs;Va;Ma;B];

% Callisto fundamental III:
clear;
fl = 38;
BDW = 18./fl; % Of the fundamental band
%f0 = 64; % MHz starting frequency of the 1st harmonic
f0 = 60; % MHz starting frequency of the fundamental
X = (BDW+1).^2;
Ma = sqrt( X.*(X+5) ./ ( 2.*(4-X) ) );
a = 6.13;  % Fitted empirical index of density variation over the heliocentric distance,
r = 1.61*695700;       % (km) Discard the last timepoint's r data because the drift rate cannot be obtained there 
% drift_rate = -0.62;    % Initial
drift_rate = -0.24298188488695013; % Mean
Vs = -2.*r/a*(1/f0).*drift_rate;   % Shock speed
Va = Vs./Ma;
B = 5.1*10^(-5).*Va.*fl;    % G
list = [Vs;Va;Ma;B];


%% 2024/5/9
%Callisto fundamental:
clear;
fl = 65;
BDW = 18./fl; % Of the fundamental band
%f0 = 64; % MHz starting frequency of the 1st harmonic
f0 = 87; % MHz starting frequency of the fundamental
X = (BDW+1).^2;
Ma = sqrt( X.*(X+5) ./ ( 2.*(4-X) ) );
a = 6.13;  % Fitted empirical index of density variation over the heliocentric distance,
r = 1.41*695700;       % (km) Discard the last timepoint's r data because the drift rate cannot be obtained there 
% drift_rate = -0.62;    % Initial
drift_rate = -0.1331255671387121; % Mean
Vs = -2.*r/a*(1/f0).*drift_rate;   % Shock speed
Va = Vs./Ma;
B = 5.1*10^(-5).*Va.*fl;    % G
list = [Vs;Va;Ma;B];

%% 2024/5/14
%Callisto fundamental-I:
% clear;
% fl = 40;
% BDW = 8./fl; % Of the fundamental band
% %f0 = 64; % MHz starting frequency of the 1st harmonic
% f0 = 50; % MHz starting frequency of the fundamental
% X = (BDW+1).^2;
% Ma = sqrt( X.*(X+5) ./ ( 2.*(4-X) ) );
% a = 6.13;  % Fitted empirical index of density variation over the heliocentric distance,
% r = 1.65*695700;       % (km) Discard the last timepoint's r data because the drift rate cannot be obtained there 
% % drift_rate = -0.62;    % Initial
% drift_rate = -0.05837987054056449; % Mean
% Vs = -2.*r/a*(1/f0).*drift_rate;   % Shock speed
% Va = Vs./Ma;
% B = 5.1*10^(-5).*Va.*fl;    % G
% list = [Vs;Va;Ma;B];

%Callisto fundamental-II:
% clear;
% fl = 37;
% BDW = 8./fl; % Of the fundamental band
% %f0 = 64; % MHz starting frequency of the 1st harmonic
% f0 = 48; % MHz starting frequency of the fundamental
% X = (BDW+1).^2;
% Ma = sqrt( X.*(X+5) ./ ( 2.*(4-X) ) );
% a = 6.13;  % Fitted empirical index of density variation over the heliocentric distance,
% r = 1.69*695700;       % (km) Discard the last timepoint's r data because the drift rate cannot be obtained there 
% % drift_rate = -0.62;    % Initial
% drift_rate = -0.10411962151939505; % Mean
% Vs = -2.*r/a*(1/f0).*drift_rate;   % Shock speed
% Va = Vs./Ma;
% B = 5.1*10^(-5).*Va.*fl;    % G
% list = [Vs;Va;Ma;B];

%Callisto fundamental-III:
% clear;
% fl = 25;
% BDW = 3./fl; % Of the fundamental band
% %f0 = 64; % MHz starting frequency of the 1st harmonic
% f0 = 29; % MHz starting frequency of the fundamental
% X = (BDW+1).^2;
% Ma = sqrt( X.*(X+5) ./ ( 2.*(4-X) ) );
% a = 6.13;  % Fitted empirical index of density variation over the heliocentric distance,
% r = 2.07*695700;       % (km) Discard the last timepoint's r data because the drift rate cannot be obtained there 
% % drift_rate = -0.62;    % Initial
% drift_rate = -0.15656712969916653; % Mean
% Vs = -2.*r/a*(1/f0).*drift_rate;   % Shock speed
% Va = Vs./Ma;
% B = 5.1*10^(-5).*Va.*fl;    % G
% list = [Vs;Va;Ma;B];
%Callisto fundamental-III update (another part of the signature mix):
% clear;
% fl = 25;
% BDW = 7./fl; % Of the fundamental band
% %f0 = 64; % MHz starting frequency of the 1st harmonic
% f0 = 31; % MHz starting frequency of the fundamental
% X = (BDW+1).^2;
% Ma = sqrt( X.*(X+5) ./ ( 2.*(4-X) ) );
% a = 6.13;  % Fitted empirical index of density variation over the heliocentric distance,
% r = 1.98*695700;       % (km) Discard the last timepoint's r data because the drift rate cannot be obtained there 
% % drift_rate = -0.07677422780635165;    % Initial
% drift_rate = -0.07677422780635165; % Mean
% Vs = -2.*r/a*(1/f0).*drift_rate;   % Shock speed
% Va = Vs./Ma;
% B = 5.1*10^(-5).*Va.*fl;    % G
% list = [Vs;Va;Ma;B];
%Callisto fundamental-III update (upper part on the left side...):
% clear;
% fl = 42;
% BDW = 6./fl; % Of the fundamental band
% f0 = 50; % MHz starting frequency of the fundamental
% X = (BDW+1).^2;
% Ma = sqrt( X.*(X+5) ./ ( 2.*(4-X) ) );
% a = 6.13;  % Fitted empirical index of density variation over the heliocentric distance,
% r = 1.675*695700;       % (km) Discard the last timepoint's r data because the drift rate cannot be obtained there 
% % drift_rate = -0.07677422780635165;    % Initial
% drift_rate = -0.1540041259385037; % Mean
% Vs = -2.*r/a*(1/f0).*drift_rate;   % Shock speed
% Va = Vs./Ma; 
% B = 5.1*10^(-5).*Va.*fl;    % G
% list = [Vs;Va;Ma;B];
% %Callisto fundamental-III update (lower part on the left side...):
% clear;
% fl = 26;
% BDW = 6./fl; % Of the fundamental band
% f0 = 30; % MHz starting frequency of the fundamental
% X = (BDW+1).^2;
% Ma = sqrt( X.*(X+5) ./ ( 2.*(4-X) ) );
% a = 6.13;  % Fitted empirical index of density variation over the heliocentric distance,
% r = 2.01*695700;       % (km) Discard the last timepoint's r data because the drift rate cannot be obtained there 
% % drift_rate = -0.07677422780635165;    % Initial
% drift_rate = -0.08611327382356368; % Mean
% Vs = -2.*r/a*(1/f0).*drift_rate;   % Shock speed
% Va = Vs./Ma;
% B = 5.1*10^(-5).*Va.*fl;    % G
% list = [Vs;Va;Ma;B];

%% 2024/5/29
% LWA Fundamental:
% BDW = [(18)/52 (18)/50 (18)/45 (18)/45 (18)/45]; % Of the fundamental band
% fl = [59 58 57 56 55.5 55 54.5];
% BDW = 18./[59 58 57 56 55.5 55 54.5]; % Of the fundamental band
% %f0 = 80; % MHz starting frequency of the 1st harmonic
% f0 = 87; % MHz starting frequency of the fundamental
% X = (BDW+1).^2;
% Ma = sqrt( X.*(X+5) ./ ( 2.*(4-X) ) );
% a = 6.13;  % Fitted empirical index of density variation over the heliocentric distance,
% r = [1.557 1.56 1.564 1.569 1.575 1.58 1.583]*695700;       % (km) Discard the last timepoint's r data because the drift rate cannot be obtained there 
% drift_rate = [-0.07 -0.07 -0.083 -0.077 -0.085 -0.077 -0.097];
% Vs = -2.*r/a*(1/f0).*drift_rate;   % Shock speed
% Va = Vs./Ma;
% B = 5.1*10^(-5).*Va.*fl;    % G

% Callisto fundamental entire band:
fl = 29;
BDW = 7./35; % Of the fundamental band
%f0 = 80; % MHz starting frequency of the 1st harmonic
f0 = 40; % MHz starting frequency of the fundamental
X = (BDW+1).^2;
Ma = sqrt( X.*(X+5) ./ ( 2.*(4-X) ) );
a = 6.13;  % Fitted empirical index of density variation over the heliocentric distance,
r = 1.8*695700;       % (km) Discard the last timepoint's r data because the drift rate cannot be obtained there 
drift_rate = -0.09575296180415997;    % mean
%drift_rate = mean(-[0.62 0.38 0.28 0.26 0.21 0.2 0.21 0.1 0.03]); % Mean
Vs = -2.*r/a*(1/f0).*drift_rate;   % Shock speed
Va = Vs./Ma;
B = 5.1*10^(-5).*Va.*fl;    % G
list = [Vs;Va;Ma;B];
% 
%% 2024/6/21
% % Callisto fundamental
% fl = [28 27.5 27 26.5 26 25.5 25];
% BDW = 12./[28 27.5 27 26.5 26 25.5 25]; % Of the fundamental band
% %f0 = 80; % MHz starting frequency of the 1st harmonic
% f0 = 35; % MHz starting frequency of the fundamental
% X = (BDW+1).^2;
% Ma = sqrt( X.*(X+5) ./ ( 2.*(4-X) ) );
% a = 6.13;  % Fitted empirical index of density variation over the heliocentric distance,
% r = [1.91 1.92 1.94 1.95 1.96 1.98 2.0]*695700;       % (km) Discard the last timepoint's r data because the drift rate cannot be obtained there 
% drift_rate = [-0.052 -0.04 -0.056 -0.039 -0.036 -0.029 -0.025];
% Vs = -2.*r/a*(1/f0).*drift_rate;   % Shock speed
% Va = Vs./Ma;
% B = 5.1*10^(-5).*Va.*fl;    % G
% list = [Vs;Va;Ma;B];

% Callisto fundamental entire band
clear;
fl = 28;
BDW = 7./32; % Of the fundamental band
%f0 = 80; % MHz starting frequency of the 1st harmonic
f0 = 39; % MHz starting frequency of the fundamental
X = (BDW+1).^2;
Ma = sqrt( X.*(X+5) ./ ( 2.*(4-X) ) );
a = 6.13;  % Fitted empirical index of density variation over the heliocentric distance,
r = 1.82*695700;       % (km) Discard the last timepoint's r data because the drift rate cannot be obtained there 
drift_rate = -0.04375550044148087;
% drift_rate = -[0.005 0.01 0.04 0.05 0.057 0.041 0.037 0.034];
% drift_rate = mean(drift_rate);      % Average drift rate
Vs = -2.*r/a*(1/f0).*drift_rate;    % Shock speed
Va = Vs./Ma;
B = 5.1*10^(-5).*Va.*fl;            % G
list = [Vs;Va;Ma;B];


%% 2024/7/4
% Callisto fundamental entire band
clear;
fl = 33;
BDW = 12./32; % Of the fundamental band
%f0 = 80; % MHz starting frequency of the 1st harmonic
f0 = 49; % MHz starting frequency of the fundamental
X = (BDW+1).^2;
Ma = sqrt( X.*(X+5) ./ ( 2.*(4-X) ) );
a = 6;  % Fitted empirical index of density variation over the heliocentric distance,
r = 1.68*695700;       % (km) Discard the last timepoint's r data because the drift rate cannot be obtained there 
drift_rate = -0.04693710372452108;
% drift_rate = -[0.005 0.01 0.04 0.05 0.057 0.041 0.037 0.034];
% drift_rate = mean(drift_rate);      % Average drift rate
Vs = -2.*r/a*(1/f0).*drift_rate;    % Shock speed
Va = Vs./Ma;
B = 5.1*10^(-5).*Va.*fl;            % G
list = [Vs;Va;Ma;B];

%% 2024/7/25
% Callisto 1st harmonic entire band
clear;
fl = 30;     % Fundamental lower
BDW = 8./fl; % Of the fundamental band
%f0 = 80; % MHz starting frequency of the 1st harmonic
f0 = 50; % MHz starting frequency of the fundamental
X = (BDW+1).^2;
Ma = sqrt( X.*(X+5) ./ ( 2.*(4-X) ) );
a = 6;  % Fitted empirical index of density variation over the heliocentric distance,
r = 1.75*695700;       % (km) Discard the last timepoint's r data because the drift rate cannot be obtained there 
drift_rate = -0.09201737693565296;
% drift_rate = -[0.005 0.01 0.04 0.05 0.057 0.041 0.037 0.034];
% drift_rate = mean(drift_rate);      % Average drift rate
Vs = -2.*r/a*(1/f0).*drift_rate;    % Shock speed
Va = Vs./Ma;
B = 5.1*10^(-5).*Va.*fl;            % G
list = [Vs;Va;Ma;B];