
%Michael A. Boateng
%Power Systems Protection:
%PROBLEM 1: MAXIMUM POSSIBLE RMS VALUE OF FAULT CURRENT THROUGH THE GENERATOR BREAKER
%PROBLEM 2: MAXIMUM MAGNETIZING CURRENT IN AMPERES AND IN PU
%PROBLEM 3: DISTANCE RELAYS WILL TRIP IN L-G AND L-L FAULTS (INVERTER INTERFACED GENERATION)
%PROBLEM 4: DIFFERENTIAL RELAY TRIPPING, CALCULATE OPERATING AND RESTRAINING CURRENTS
%PROBLEM 5: DISTANCE RELAYS DESIGN (ZONE 1 & ZONE 2...0.7_1.2)

%!!Readme: Kindly run the code in sections...Each section is allocated to a
%specific Problem part!!

%%
%PROBLEM 1 A:
%Mostly drawn and solved on paper
P=(600e06)/3; V=(18e03/sqrt(3));
Ib= (P)/(V);SS

%Double current division to find IG2', to find IG2, referred to the left
%side, you need to multiply IG2' by e^30j
%Now, current division TWICE:

Rother=0.25j; Rme=9.949e-03+0.2474i;
Rother_2=0.25j; Rme_2=0.02154+0.3485i;   %IMPEDANCES FROM THE DIFFERENT BRANCHES ON MY DIAGRAMS
I_total=1.1742-1.81i;

I_first_division= (Rother/(Rother + Rme))*I_total;
I_second_division= (Rother_2/(Rother_2 + Rme_2))*I_first_division;

%From sequence to direct phase quantities and back
%IA, IB, IC from I1, I2, I0
I1 = I_total*(exp(j*deg2rad(-30)));
I2 = I_second_division*(exp(j*deg2rad(30)));

convert_complex_to_polar(I_total)

disp('I1');
convert_complex_to_polar(I1)

disp('I2');
convert_complex_to_polar(I2)
I0 = 0;
a = exp(j*2*pi/3);

I_A = I1*1 + I2*1 + I0*1;
I_B = I1*a^2 + I2*a + I0*1; % Corrected application of symmetrical components
I_C = I1*a + I2*a^2 + I0*1;

I_A =I_A*Ib
disp('I_A');
convert_complex_to_polar(I_A)

I_B =I_B*Ib
disp('I_B');
convert_complex_to_polar(I_B)

I_C =I_C*Ib
disp('I_C');
convert_complex_to_polar(I_C)

disp(' ');

R = 0.02405;   %Zth_real
X = 0.4629;   %Zth_imag
L  = X/(2 * pi * 60);
t = 6*(1/60); % time

% make Ifault the magnitude of IA_G, IB_G, and IC_G --> I_A, I_B, & I_C
Ifault_phasor_A = abs(I_A);
Ifault_phasor_B = abs(I_B);
Ifault_phasor_C = abs(I_C);


%Compute the maximum possible rms value A,B,C
Irms_a = sqrt((Ifault_phasor_A)^2 + (2 * (Ifault_phasor_A)^2) * exp(-(2*R)/L*t));
Irms_mag_a = abs(Irms_a);
Irms_angle_a = angle(Irms_a) * 180/pi;
fprintf('Max rms value magnitude (IA): %.4f A\n', Irms_mag_a);
fprintf('Max rms value angle: %.2f degrees\n', Irms_angle_a);

Irms_b = sqrt((Ifault_phasor_B)^2 + (2 * (Ifault_phasor_B)^2) * exp(-(2*R)/L*t));
Irms_mag_b = abs(Irms_b);
Irms_angle_b = angle(Irms_b) * 180/pi;
fprintf('Max_ rms value magnitude (IB): %.4f A\n', Irms_mag_b);
fprintf('Max rms value angle: %.2f degrees\n', Irms_angle_b);

Irms_c = sqrt((Ifault_phasor_C)^2 + (2 * (Ifault_phasor_C)^2) * exp(-(2*R)/L*t));
Irms_mag_c = abs(Irms_c);
Irms_angle_c = angle(Irms_c) * 180/pi;
fprintf('Max rms value magnitude(IC): %.4f A\n', Irms_mag_c);
fprintf('Max rms value angle: %.2f degrees\n', Irms_angle_c);

max_num = max([Irms_mag_a, Irms_mag_b, Irms_mag_c]);
fprintf('FINAL Max rms value magnitude: %.4f A\n', max_num);

%PROBLEM 1 B:
%I2 in pu, t seconds
t=(4.5)/((abs(I2))^2)

%Magnitude of IG2 (That is negative sequence IG is placed in bracket ()^2)

%%
%Question 2A:

% Given parameters
t=0;        % The transformer is energized at time<--
E = 115e3;  % Voltage in volts
omega = 377;% Angular frequency in rad/s
i0 = 0.001; % Base magnetizing current in pu

i0= 0.001*(8e06/115e03) %pu to A, P and V are in single phase
lambda0 = 1.0;  % Base flux linkage in pu 
%Although lambda0 is given, it can also be obtained (its actual and not pu)
lambda0 = ((sqrt(2)*E)/omega)

% Function to calculate flux linkage
ext = lambda0* sin(omega * t + deg2rad(45));  % Wb

lambda_max = (lambda0*(-1))-(ext) 

% Calculate the maximum magnetizing current
i_max = -i0 * (abs(lambda_max) / lambda0)^13;  % Amperes

disp(['Maximum magnetizing current (Amperes): ', num2str(i_max)]);
disp(['Maximum magnetizing current (pu): ', num2str(i_max/(8e06/115e03))]);

%%
%PROBLEM 2B:

% Given parameters
t=0.00625;        % The transformer is energized at time<--
E = 115e3;  % Voltage in volts
omega = 377;% Angular frequency in rad/s
i0 = 0.001; % Base magnetizing current in pu

i0= 0.001*(8e06/115e03) %pu to A, P and V are in single phase
lambda0 = 1.0;  % Base flux linkage in pu 
%Although lambda0 is given, it can also be obtained (its actual and not pu)
lambda0 = ((sqrt(2)*E)/omega)

% Function to calculate flux linkage
ext = lambda0* sin(omega * t + deg2rad(45));  % Wb

lambda_max = (lambda0*(-1))-(ext) 

% Calculate the maximum magnetizing current
i_max = -i0 * (abs(lambda_max) / lambda0)^13;  % Amperes

disp(['Maximum magnetizing current (Amperes): ', num2str(i_max)]);
disp(['Maximum magnetizing current (pu): ', num2str(i_max/(8e06/115e03))]);

%%
%PROBLEM 3 A1 [LEFT SIDE]:
%OBSERVING NEGATIVE AND ZERO SEQUENCES AT BOTH ENDS

%LG (Single Line to GND):
Ian=524.7*(exp(j*deg2rad(41.49)));
Ibn=485.6*(exp(j*deg2rad(31.53)));
Icn=3313*(exp(j*deg2rad(43.47)));

a = exp(j*2*pi/3);
I1_positive=(1/3)*[1 a a^2]* [Ian; Ibn; Icn];
disp('I1_positive');
convert_complex_to_polar(I1_positive)
disp(' ');

a = exp(j*2*pi/3);
I2_negative=(1/3)*[1 a^2 a]* [Ian; Ibn; Icn];
disp('I2_negative');
convert_complex_to_polar(I2_negative)
disp(' ');

I0_zero=(1/3)*[1 1 1]* [Ian; Ibn; Icn];
disp('I0_zero');
convert_complex_to_polar(I0_zero)
disp(' ');

%%
%PROBLEM 3 A2 [RIGHT SIDE]:
%OBSERVING NEGATIVE AND ZERO SEQUENCES AT BOTH ENDS


%LG (Single Line to GND):
Ian=436.6*(exp(j*deg2rad(-149.03)));
Ibn=436*(exp(j*deg2rad(-134.82)));
Icn=2019*(exp(j*deg2rad(42.53)));

a = exp(j*2*pi/3);
I1_positive=(1/3)*[1 a a^2]* [Ian; Ibn; Icn];
disp('I1_positive');
convert_complex_to_polar(I1_positive)
disp(' ');

a = exp(j*2*pi/3);
I2_negative=(1/3)*[1 a^2 a]* [Ian; Ibn; Icn];
disp('I2_negative');
convert_complex_to_polar(I2_negative)
disp(' ');

I0_zero=(1/3)*[1 1 1]* [Ian; Ibn; Icn];
disp('I0_zero');
convert_complex_to_polar(I0_zero)
disp(' ')

%OBSERVATION:
disp(['We received information indicating the existence of inverter-interfaced power sources, which, as we understand, lead to the suppression of negative and zero sequence currents. ' ...
    'This phenomenon occurred in our case too. Initially, during a single line to ground fault, ' ...
    'the left side exhibited a notably high level of zero sequence current. However, on the right side, it was significantly suppressed, posing challenges for detection. ' ...
    'This suppression could potentially cause certain protective relays to malfunction'])
%%
%PROBLEM 3 B1 [LEFT SIDE]:

%[DISTANCE/IMPEDANCE/RATIO RELAY]
%LG (Single Line to GND):
%zone1: 0.8--4cycles or INSTANTANEOUS, zone2:1.25--15 cycles, zone3:2.08--35 cycles
z1= 2.97+46.475j; z2=2.97+46.475j; z0=39.364+148.03j; %z1 in ohms
length_of_line=78; %miles

Z_per_l= (z1/length_of_line);
m=2.32*(exp(j*deg2rad(-16.1))); %compensation factor is given
%m=((z0-z1)/z1);                %compensation factor must be found

PT=115/220000; CT=5/2400;
impedance_ratio=PT/CT;

Vcn=191700*(exp(j*deg2rad(120.17)));

Ian=524.7*(exp(j*deg2rad(41.49)));
Ibn=485.6*(exp(j*deg2rad(31.53)));
Icn=3313*(exp(j*deg2rad(43.47)));

Io=(1/3)*[1 1 1]* [Ian; Ibn; Icn]; %FINAL ANSWER IN RECTANGULAR

%Where is the fault located?...A,B,or C?
z_imp_line=(Vcn)/(Icn + (m*Io)) %z*L

disp('z_imp_line');
convert_complex_to_polar(z_imp_line)

%what the line and what the relay sees
z_imp_relay=z_imp_line*impedance_ratio
disp('z_imp_relay');
convert_complex_to_polar(z_imp_relay)

%The fault distance from side 1 seen by the relay at side 1 [LEFT]:
length_to_fault_1=(z_imp_line)/(Z_per_l)
disp('length_to_fault');
convert_complex_to_polar(length_to_fault_1)

%z2_relay_zone2=14.5*(exp(j*deg2rad(86.3)));

%%
%PROBLEM 3 B2 [RIGHT SIDE]:

%[DISTANCE/IMPEDANCE/RATIO RELAY]
%LG (Single Line to GND):
%zone1: 0.8--4cycles or INSTANTANEOUS, zone2:1.25--15 cycles, zone3:2.08--35 cycles

z1= 2.97+46.475j; z2=2.97+46.475j; z0=39.364+148.03j; %z1 in ohms
length_of_line=78; %miles

Z_per_l= (z1/length_of_line);
m=2.32*(exp(j*deg2rad(-16.1))); %compensation factor is given
%m=((z0-z1)/z1);                %compensation factor must be found

PT=115/220000; CT=5/2400;
impedance_ratio=PT/CT;

Vcn=53550*(exp(j*deg2rad(121.97)));

Ian=436.6*(exp(j*deg2rad(-149.03)));
Ibn=436*(exp(j*deg2rad(-134.82)));
Icn=2019*(exp(j*deg2rad(42.53)));

Io=(1/3)*[1 1 1]* [Ian; Ibn; Icn]; %FINAL ANSWER IN RECTANGULAR

%Where is the fault located?...A,B,or C?
z_imp_line=(Vcn)/(Icn + (m*Io)) %z*L

disp('z_imp_line');
convert_complex_to_polar(z_imp_line)

%what the line and what the relay sees
z_imp_relay=z_imp_line*impedance_ratio
disp('z_imp_relay');
convert_complex_to_polar(z_imp_relay)
disp(' ');

%The fault distance from side 1 seen by the relay at side 2 [RIGHT]:
length_to_fault_2=(z_imp_line)/(Z_per_l);
disp('length_to_fault_2');
convert_complex_to_polar(length_to_fault_2)

length_from_side1_relay2= (length_of_line) - (length_to_fault_2)
disp('length_from_side1_relay2');
convert_complex_to_polar(length_from_side1_relay2)

%z2_relay_zone2=14.5*(exp(j*deg2rad(86.3)));

%PROBLEM 3 C
difference_in_fault_distance= (length_to_fault_1 - length_from_side1_relay2)
disp('difference_in_fault_distance');
convert_complex_to_polar(difference_in_fault_distance)

%%
%PROBLEM 3 D1
%21L (Line to Line) [THE LEFT SIDE]
z1= 2.97+46.475j; z2=2.97+46.475j; z0=39.364+148.03j; %z1 in ohms
length_of_line=78; %miles

Z_per_l= (z1/length_of_line);
m=2.32*(exp(j*deg2rad(-16.1))); %compensation factor is given
%m=((z0-z1)/z1);                %compensation factor must be found

PT=115/220000; CT=5/2400;
impedance_ratio=PT/CT;

Vbn=190400*(exp(j*deg2rad(-121.77)));
Vcn=189900*(exp(j*deg2rad(121.45)));

Ibn=5100*(exp(j*deg2rad(-176.43)));
Icn=5042*(exp(j*deg2rad(2.83)));

%Where is the fault located?...A,B,or C?
z_imp_line=(Vbn-Vcn)/(Ibn-Icn); %z*L\
disp('z_imp_line');
convert_complex_to_polar(z_imp_line)

%what the line and what the relay sees
z_imp_relay=z_imp_line*impedance_ratio
disp('z_imp_relay');
convert_complex_to_polar(z_imp_relay)
disp(' ');

length_to_fault_3=(z_imp_line)/(Z_per_l)
disp('length_to_fault_3');
convert_complex_to_polar(length_to_fault_3)

%z2_relay_zone2=14.5*(exp(j*deg2rad(86.3)));


%%
%PROBLEM 3 D2
%21L (Line to Line) [THE RIGHT SIDE]
z1= 2.97+46.475j; z2=2.97+46.475j; z0=39.364+148.03j; %z1 in ohms
length_of_line=78; %miles

Z_per_l= (z1/length_of_line);
m=2.32*(exp(j*deg2rad(-16.1))); %compensation factor is given
%m=((z0-z1)/z1);                %compensation factor must be found

PT=115/220000; CT=5/2400;
impedance_ratio=PT/CT;

Vbn=135900*(exp(j*deg2rad(-138.35)));
Vcn=135000*(exp(j*deg2rad(136.67)));

Ibn=4567*(exp(j*deg2rad(-177.27)));
Icn=4618*(exp(j*deg2rad(2.02)));

%Where is the fault located?...A,B,or C?
z_imp_line=(Vbn-Vcn)/(Ibn-Icn); %z*L
disp('z_imp_line');
convert_complex_to_polar(z_imp_line)

%what the line and what the relay sees
z_imp_relay=z_imp_line*impedance_ratio
disp('z_imp_relay');
convert_complex_to_polar(z_imp_relay)
disp(' ');

%The fault distance from side 1 seen by the relay at side 2 [RIGHT]:
length_to_fault_4=(z_imp_line)/(Z_per_l);
disp('length_to_fault_4');
convert_complex_to_polar(length_to_fault_4)


length_from_side1_relay2_2= (length_of_line) - (length_to_fault_4)
disp('length_from_side1_relay2_2');
convert_complex_to_polar(length_from_side1_relay2_2)

%z2_relay_zone2=14.5*(exp(j*deg2rad(86.3)));

%PROBLEM 3 E
difference_in_fault_distance_2= abs(length_to_fault_3 - length_from_side1_relay2_2)
disp('difference_in_fault_distance_2');
convert_complex_to_polar(difference_in_fault_distance_2)

%%
%PROBLEM 4A
clc
k1=5/100; k2=5/2000;
percentage_restrain=25;
Ia=150.9*(exp(j*deg2rad(-28.6))); Ib=115.6*(exp(j*deg2rad(-141.27))); Ic=150.6*(exp(j*deg2rad(106.31)));
Ia2=2364*(exp(j*deg2rad(129.16))); Ib2=1479*(exp(j*deg2rad(8.97))); Ic2=1486*(exp(j*deg2rad(-111.26)));

left_side_1=k1*Ia
disp('left_side_1');
convert_complex_to_polar(left_side_1)

left_side_2=k1*Ib
disp('left_side_2');
convert_complex_to_polar(left_side_2)

left_side_3=k1*Ic
disp('left_side_3');
convert_complex_to_polar(left_side_3)

right_side_1=(k2*Ia2)-(k2*Ib2)
disp('right_side_1');
convert_complex_to_polar(right_side_1)

right_side_2=(k2*Ib2)-(k2*Ic2)
disp('right_side_2');
convert_complex_to_polar(right_side_2)

right_side_3=(k2*Ic2)-(k2*Ia2)
disp('right_side_3');
convert_complex_to_polar(right_side_3)

%4A1
Ioa= (left_side_1)+(right_side_1)
disp('Ioa');
convert_complex_to_polar(Ioa)

Iob= (left_side_2)+(right_side_2)
disp('Iob');
convert_complex_to_polar(Iob)

Ioc= (left_side_3)+(right_side_3)
disp('Ioc');
convert_complex_to_polar(Ioc)

%4A2
percent_1=abs(Ioa)/(0.5*(abs(left_side_1)+abs(right_side_1)))*100
if percentage_restrain > percent_1
    disp('Differential relay will not trip the transformer');
elseif percentage_restrain < percent_1
    disp('Differential relay will not trip the transformer');
else
    disp('percentages are equal');
end

percent_2=abs(Iob)/(0.5*(abs(left_side_2)+abs(right_side_2)))*100
if percentage_restrain > percent_2
    disp('Differential relay will not trip the transformer');
elseif percentage_restrain < percent_2
    disp('Differential relay will not trip the transformer');
else
    disp('percentages are equal');
end

percent_3=abs(Ioc)/(0.5*(abs(left_side_3)+abs(right_side_3)))*100
if percentage_restrain > percent_3
    disp('Differential relay will not trip the transformer');
elseif percentage_restrain < percent_3
    disp('Differential relay will not trip the transformer');
else
    disp('percentages are equal');
end

%%
clc
%PROBLEM 4B
k1=5/100; k2=5/2000;
percentage_restrain=25;
Ia=220.8*(exp(j*deg2rad(-27.7))); Ib=192*(exp(j*deg2rad(-145.91))); Ic=213.5*(exp(j*deg2rad(99.81)));
Ia2=2590*(exp(j*deg2rad(123.89))); Ib2=2160*(exp(j*deg2rad(7.2))); Ic2=1872*(exp(j*deg2rad(-114.22)));

left_side_1=k1*Ia
disp('left_side_1');
convert_complex_to_polar(left_side_1)

left_side_2=k1*Ib
disp('left_side_2');
convert_complex_to_polar(left_side_2)

left_side_3=k1*Ic
disp('left_side_3');
convert_complex_to_polar(left_side_3)

right_side_1=(k2*Ia2)-(k2*Ib2)
disp('right_side_1');
convert_complex_to_polar(right_side_1)

right_side_2=(k2*Ib2)-(k2*Ic2)
disp('right_side_2');
convert_complex_to_polar(right_side_2)

right_side_3=(k2*Ic2)-(k2*Ia2)
disp('right_side_3');
convert_complex_to_polar(right_side_3)

%4A1
Ioa= (left_side_1)+(right_side_1)
disp('Ioa');
convert_complex_to_polar(Ioa)

Iob= (left_side_2)+(right_side_2)
disp('Iob');
convert_complex_to_polar(Iob)

Ioc= (left_side_3)+(right_side_3)
disp('Ioc');
convert_complex_to_polar(Ioc)

%4B2
percent_1=abs(Ioa)/(0.5*(abs(left_side_1)+abs(right_side_1)))*100
if percentage_restrain > percent_1
    disp('Differential relay will not trip the transformer');
elseif percentage_restrain < percent_1
    disp('Differential relay will not trip the transformer');
else
    disp('percentages are equal');
end

percent_2=abs(Iob)/(0.5*(abs(left_side_2)+abs(right_side_2)))*100
if percentage_restrain > percent_2
    disp('Differential relay will not trip the transformer');
elseif percentage_restrain < percent_2
    disp('Differential relay will not trip the transformer');
else
    disp('percentages are equal');
end

percent_3=abs(Ioc)/(0.5*(abs(left_side_3)+abs(right_side_3)))*100
if percentage_restrain > percent_3
    disp('Differential relay will not trip the transformer');
elseif percentage_restrain < percent_3
    disp('Differential relay will not trip the transformer');
else
    disp('percentages are equal');
end


%%
%PROBLEM 5,Q1, Q2, Q3
%[DISTANCE/IMPEDANCE/RATIO RELAY]
%LG (Single Line to GND):
%zone1: 0.8--4cycles or INSTANTANEOUS, zone2:1.25--15 cycles, zone3:2.08--35 cycles

z1= 2.222+39.134j; z2=2.222+39.134j; z0=24.564+158.775j; %z1 in ohms
length_of_line=59.5; %miles

Z_per_l= (z1/length_of_line);
%m=2.32*(exp(j*deg2rad(-16.1))); %compensation factor is given

m=3.2 %((z0-z1)/z1))
disp('m');
convert_complex_to_polar(m)
disp(' ');

%compensation factor must be found. BUT REMEMBER TO CHANGE TO WHOLE NO. OF PROBLEM ASKS
%SAME GOES FOR THE ANGLES

PT=115/142000; CT=5/1500;
impedance_ratio=PT/CT;

Van=85440*(exp(j*deg2rad(25.19)));

Ian=1708*(exp(j*deg2rad(-49.97)));
Ibn=113.8*(exp(j*deg2rad(113)));
Icn=116.3*(exp(j*deg2rad(147)));

Io=(1/3)*[1 1 1]* [Ian; Ibn; Icn]; %FINAL ANSWER IN RECTANGULAR

%Where is the fault located?...A,B,or C?
z_imp_line=(Van)/(Ian + (m*Io)); %z*L


disp('PROBLEM 1');
%z1_relay_zone1=9.3*(exp(j*deg2rad(86.3)));   
z1_relay_zone1=impedance_ratio*0.7*z1
disp('z1_relay_zone1');
convert_complex_to_polar(z1_relay_zone1)
disp(' ');

time_delay_zone_1 = 'instantaneous'
%impedance ratio x percentage x z1


disp('PROBLEM 2');
%z2_relay_zone2=14.5*(exp(j*deg2rad(86.3)));
z2_relay_zone2=impedance_ratio*1.2*z1
disp('z2_relay_zone2');
convert_complex_to_polar(z2_relay_zone2)
disp(' ');

time_delay_zone_2 = '15 cycles (standard)'
%impedance ratio x percentage x z1

%z3_relay_zone2=14.5*(exp(j*deg2rad(86.3)));
z3_relay_zone3=impedance_ratio*0.7*z1;
%impedance ratio x percentage x z1

length_to_fault=(z_imp_line)/(Z_per_l);
%z2_relay_zone2=14.5*(exp(j*deg2rad(86.3)));


disp('PROBLEM 3');
%what the line and what the relay sees
disp('z_imp_line');
convert_complex_to_polar(z_imp_line)
disp(' ');

z_imp_relay=z_imp_line*impedance_ratio
disp('z_imp_relay');
convert_complex_to_polar(z_imp_relay)
disp(' ');


%%
%CONVERT COMPLEX TO POLAR WHEN NEEDED

function convert_complex_to_polar(complex_number)
    % Extract rectangular coordinates
    x = real(complex_number);
    y = imag(complex_number);

    % Calculate magnitude using Pythagoras theorem
    magnitude = sqrt(x^2 + y^2);

    % Calculate angle using atan2 function (in radians)
    angle_rad = atan2(y, x);

    % Convert angle from radians to degrees
    angle_deg = rad2deg(angle_rad);

    % Display results
    disp(['Rectangular coordinates (x, y): ', num2str(x), ', ', num2str(y)]);
    disp(['Polar coordinates (magnitude, angle degrees): ', num2str(magnitude), ', ', num2str(angle_deg)]);
end