close all
clear all
clc
format long

load mimu_gnss1.mat;
truth = load('F:\i2nav_shuju-main\truth16465.txt');

%% sensor_data
ts = 1/200;
re = 6378137;
Wie = 7.2921151467e-5;
e = 1/298.257223563;
b1 = 6356752.3142;
e11 = sqrt(re^2-b1^2)/re;
e12 = sqrt(re^2-b1^2)/b1;
c1 = re^2/b1;
Po = 180.*60.*60/pi;
lb = [-0.073;0.302;0.087];

Ti_mimu = mimu_gnss(:,1);
% m_wx = mimu_gnss(:,2) - 0.0232.*(pi/180).*ts;
% m_wy = mimu_gnss(:,3) - 0.0575.*(pi/180).*ts;
% m_wz = mimu_gnss(:,4) - 0.0450.*(pi/180).*ts;
m_wx = mimu_gnss(:,2);
m_wy = mimu_gnss(:,3);
m_wz = mimu_gnss(:,4);
m_ax = mimu_gnss(:,5);
m_ay = mimu_gnss(:,6);
m_az = mimu_gnss(:,7);
GNSS_Lat = mimu_gnss(:,8).*(pi/180);
GNSS_Lon = mimu_gnss(:,9).*(pi/180);
GNSS_H = mimu_gnss(:,10);
GNSS_vx = mimu_gnss(:,11);
GNSS_vy = mimu_gnss(:,12);
GNSS_vz = mimu_gnss(:,13);
[hang1,~] = size(mimu_gnss);
Ti_truth = truth(:,2) - truth(1,2);
truth_Lat = truth(:,3).*(pi/180);
truth_Lon = truth(:,4).*(pi/180);
truth_H = truth(:,5);
truth_vx = truth(:,6);
truth_vy = truth(:,7);
truth_vz = truth(:,8);
truth_roll = truth(:,9);
truth_pitch = truth(:,10);
truth_yaw = truth(:,11);
[hang2,~] = size(truth);

BB0 = 1-(3/4).*e12^2+(45/64).*e12^4-(175/256).*e12^6+(11025/16384).*e12^8;
BB2 = BB0-1;
BB4 = (15/32).*e12^4-(175/384).*e12^6+(3675/8192).*e12^8;
BB6 = (-35/96).*e12^6+(735/2048).*e12^8;
BB8 = (315/1024).*e12^8;

P_kf = diag([(0.5.*(pi/180))^2 (0.5.*(pi/180))^2 (1.*(pi/180))^2 0.0002^2 0.0002^2 0.0002^2 (0.01/3600/60/57.3)^2 (0.01/3600/60/57.3)^2 (0.01)^2  (25.*(pi/180)/3600)^2 (25.*(pi/180)/3600)^2 (25.*(pi/180)/3600)^2 (0.20408e-3*9.7803267714)^2 (0.20408e-3*9.7803267714)^2 (0.20408e-3*9.7803267714)^2]);
Q_diag = [(15.*0.15*(2.909*1e-4))^2 (15.*0.15*(2.909*1e-4))^2 (15.*0.15*(2.909*1e-4))^2 (15.*0.0002)^2 (15.*0.0002)^2 (15.*0.0002)^2];
Q = diag(Q_diag);
R_diag = [(0.001)^2 (0.001)^2 (0.001)^2 (0.01/3600/60/57.3)^2 (0.01/3600/60/57.3)^2 (0.01)^2];
R_kf = diag(R_diag);
Xkf_attitude = [];
Xkf_velocity = [];
Xkf_position = [];
Xkf_gyroscope = [];
Xkf_accelerometer = [];
Mis_gyroscope = [0;0;0];
Mis_accelerometer = [0;0;0];
MisAngle_attitude = [0;0;0];
MisV_velocity = [0;0;0];
MisP_position = [0;0;0];
X_kf = [MisAngle_attitude;MisV_velocity;MisP_position;Mis_gyroscope;Mis_accelerometer];
s = 1;
P11 = P_kf(1,1);
P22 = P_kf(2,2);
P33 = P_kf(3,3);

%% Alignment
g = [0;0;9.79361];
T = 1;
alpha = 0;beta1 = 0;
Cnb = eye(3);Cbb = eye(3);Cnn = eye(3);
Cbb_gnss = eye(3);
v1 = 0;v2 = 0;
seita1 = 0;seita2 = 0;
vn0 = [GNSS_vx(1);GNSS_vy(1);GNSS_vz(1)];
vn = [GNSS_vx(1);GNSS_vy(1);GNSS_vz(1)];
A = 0;
fCpitch(1,1) = 0;
fCroll(1,1) = 0;
fCyaw(1,1) = 0;
Cnnm = eye(3);
[fh,~] = find(Ti_mimu == 10);
i_kf = fh;

mimu_Lat(1,1) = mimu_gnss(fh,8).*(pi/180);
mimu_Lon(1,1) = mimu_gnss(fh,9).*(pi/180);
mimu_H(1,1) = mimu_gnss(fh,10);
mimu_vx(1,1) = mimu_gnss(fh,11);
mimu_vy(1,1) = mimu_gnss(fh,12);
mimu_vz(1,1) = mimu_gnss(fh,13);
V_m = [];

for i = 2:fh
    % upgrade Cbb
    % two sample + previous cycle
    Angle1 = [m_wx(i-1);m_wy(i-1);m_wz(i-1)];
    Angle2 = [m_wx(i);m_wy(i);m_wz(i)];
    A_b = Angle2 + cross(((1/12).*Angle1),Angle2);
    
    A_b_X = [0 -A_b(3) A_b(2);
             A_b(3) 0 -A_b(1);
             -A_b(2) A_b(1) 0];
    IA_bI = norm(A_b);
    
    Cbbt = eye(3) + (sin(IA_bI)/IA_bI)*A_b_X + ((1-cos(IA_bI))/(IA_bI)^2)*(A_b_X)^2;
    Cbb = Cbb*Cbbt;
    
    v1 = v1 + [m_ax(i);m_ay(i);m_az(i)]/2;
    v2 = v2 + [m_ax(i);m_ay(i);m_az(i)]/2;
    seita1 = seita1 + [m_wx(i);m_wy(i);m_wz(i)]/2;
    seita2 = seita2 + [m_wx(i);m_wy(i);m_wz(i)]/2;
    
    if GNSS_Lat(i) ~= 0
        Rx = (re/sqrt(1-e11^2.*(sin(GNSS_Lat(i)))^2));
        Ry = ((Rx.*(1-e11^2))/(1-e11^2.*(GNSS_Lat(i))^2));
        Wnie = [Wie.*cos(GNSS_Lat(i));
                0;
                -Wie.*sin(GNSS_Lat(i))];
        Win = [(Wie.*cos(GNSS_Lat(i))+GNSS_vy(i))/(Rx+GNSS_H(i));
               -GNSS_vx(i)/(Ry+GNSS_H(i));
               (-Wie.*sin(GNSS_Lat(i))-GNSS_vy(i).*tan(GNSS_Lat(i)))/(Rx+GNSS_H(i))];
        Win_X = [0 -Win(3) Win(2);
                 Win(3) 0 -Win(1);
                 -Win(2) Win(1) 0];
        A_n = Win.*T;
        A_n_X = [0 -A_n(3) A_n(2);
                 A_n(3) 0 -A_n(1);
                 -A_n(2) A_n(1) 0];
        IA_nI = norm(A_n);
        Cnnt = eye(3) - (sin(IA_nI)/IA_nI)*A_n_X + ((1-cos(IA_nI))/(IA_nI)^2)*(A_n_X)^2;
        % upgrade alpha
        alpha = alpha + Cbb_gnss*(v1 + v2 + 0.5.*cross((seita1 + seita2),(v1 + v2)) + (2/3).*(cross(seita1,v2) + cross(v1,seita2)));
        v1 = 0;v2 = 0;
        seita1 = 0;seita2 = 0;
        Cbb_gnss = Cbb;
        Am1 = cross((((T/2).*eye(3) + (T^2/6).*Win_X)*Wnie),vn);
        vn = [GNSS_vx(i);GNSS_vy(i);GNSS_vz(i)];
        Am2 = cross((((T/2).*eye(3) + (T^2/3).*Win_X)*Wnie),vn);
        beta1 = beta1 + Cnn'*(Am1 + Am2 - (T.*eye(3) + (T^2/2).*Win_X)*g);
        % upgrade Cnn
        Cnn = Cnnt*Cnn;
        % upgrade beta
        beta = Cnn'*vn - Cnnm'*vn0 + beta1;
        
        A = A + beta*alpha';
        B = A'*A;
        [Eigen_vector,D] = eig(B);
        D = sqrt(D);
        V = Eigen_vector;
        U = A*V*inv(D);
        Cnb0 = U*V';
        Cnb = Cnn*Cnb0*Cbb;
        
        fCpitch(i,1) = asin(-Cnb(3,1));
        fCroll(i,1) = atan(Cnb(3,2)/Cnb(3,3));
        if Cnb(3,3) < 0
            if fCroll(i,1) > 0
                fCroll(i,1) = fCroll(i,1) - pi;
            else
                fCroll(i,1) = fCroll(i,1) + pi;
            end
        end
        fCyaw(i,1) = atan(Cnb(2,1)/Cnb(1,1));
        if Cnb(2,1) > 0
            if abs(Cnb(1,1)) < 0.0001
                fCyaw(i,1) = pi/2;
            elseif Cnb(1,1) < 0
                fCyaw(i,1) = fCyaw(i,1) + pi;
            end
        else
            if abs(Cnb(1,1)) < 0.0001
                fCyaw(i,1) = -pi/2;
            elseif Cnb(1,1) < 0
                fCyaw(i,1) = fCyaw(i,1) - pi;
            end
        end
        % yaw range: 0°-360°
        if fCyaw(i,1) < 0
            fCyaw(i,1) = fCyaw(i,1) + 2.*pi;
        end
    else
        Cnb = Cnb*Cbbt;
        
        fCpitch(i,1) = asin(-Cnb(3,1));
        fCroll(i,1) = atan(Cnb(3,2)/Cnb(3,3));
        if Cnb(3,3) < 0
            if fCroll(i,1) > 0
                fCroll(i,1) = fCroll(i,1) - pi;
            else
                fCroll(i,1) = fCroll(i,1) + pi;
            end
        end
        fCyaw(i,1) = atan(Cnb(2,1)/Cnb(1,1));
        if Cnb(2,1) > 0
            if abs(Cnb(1,1)) < 0.0001
                fCyaw(i,1) = pi/2;
            elseif Cnb(1,1) < 0
                fCyaw(i,1) = fCyaw(i,1) + pi;
            end
        else
            if abs(Cnb(1,1)) < 0.0001
                fCyaw(i,1) = -pi/2;
            elseif Cnb(1,1) < 0
                fCyaw(i,1) = fCyaw(i,1) - pi;
            end
        end
        % yaw range: 0°-360°
        if fCyaw(i,1) < 0
            fCyaw(i,1) = fCyaw(i,1) + 2.*pi;
        end
    end
end

%% Attitude_Vimu_Pimu_ESKF_upgrade
for i = fh+1:hang1
    % Attitide
    if GNSS_vx(i) ~= 0
        % upgrade Cnn
        Rx = (re/sqrt(1-e11^2.*(sin(GNSS_Lat(i)))^2));
        Ry = ((Rx.*(1-e11^2))/(1-e11^2.*(GNSS_Lat(i))^2));
        Win = [(Wie.*cos(GNSS_Lat(i))+GNSS_vy(i))/(Rx+GNSS_H(i));
               -GNSS_vx(i)/(Ry+GNSS_H(i));
               (-Wie.*sin(GNSS_Lat(i))-GNSS_vy(i).*tan(GNSS_Lat(i)))/(Rx+GNSS_H(i))];
        A_n = Win*T;
        A_n_X = [0 -A_n(3) A_n(2);
                 A_n(3) 0 -A_n(1);
                 -A_n(2) A_n(1) 0];
        IA_nI = norm(A_n);
        Cnn = eye(3) - (sin(IA_nI)/IA_nI)*A_n_X + ((1-cos(IA_nI))/(IA_nI)^2)*(A_n_X)^2;
           
        % upgrade Cbb
        % two sample + previous cycle
        Angle1 = [m_wx(i-1);m_wy(i-1);m_wz(i-1)];
        Angle2 = [m_wx(i);m_wy(i);m_wz(i)];
        A_b = Angle2 + cross(((1/12).*Angle1),Angle2);
           
        A_b_X = [0 -A_b(3) A_b(2);
                 A_b(3) 0 -A_b(1);
                 -A_b(2) A_b(1) 0];
        IA_bI = norm(A_b);
        Cbb = eye(3) + (sin(IA_bI)/IA_bI)*A_b_X + ((1-cos(IA_bI))/(IA_bI)^2)*(A_b_X)^2;
           
        Cnb = Cnn*Cnb*Cbb;
        fCpitch(i,1) = asin(-Cnb(3,1));
        fCroll(i,1) = atan(Cnb(3,2)/Cnb(3,3));
        if Cnb(3,3) < 0
            if fCroll(i,1) > 0
                fCroll(i,1) = fCroll(i,1) - pi;
            else
                fCroll(i,1) = fCroll(i,1) + pi;
            end
        end
        fCyaw(i,1) = atan(Cnb(2,1)/Cnb(1,1));
        if Cnb(2,1) > 0
            if abs(Cnb(1,1)) < 0.0001
                fCyaw(i,1) = pi/2;
            elseif Cnb(1,1) < 0
                fCyaw(i,1) = fCyaw(i,1) + pi;
            end
        else
            if abs(Cnb(1,1)) < 0.0001
                fCyaw(i,1) = -pi/2;
            elseif Cnb(1,1) < 0
                fCyaw(i,1) = fCyaw(i,1) - pi;
            end
        end
        % yaw range: 0°-360°
        if fCyaw(i,1) < 0
            fCyaw(i,1) = fCyaw(i,1) + 2.*pi;
        end
    else
        % upgrade Cbb
        % two sample + previous cycle
        Angle1 = [m_wx(i-1);m_wy(i-1);m_wz(i-1)];
        Angle2 = [m_wx(i);m_wy(i);m_wz(i)];
        A_b = Angle2 + cross(((1/12).*Angle1),Angle2);
           
        A_b_X = [0 -A_b(3) A_b(2);
                 A_b(3) 0 -A_b(1);
                 -A_b(2) A_b(1) 0];
        IA_bI = norm(A_b);
        Cbb = eye(3) + (sin(IA_bI)/IA_bI)*A_b_X + ((1-cos(IA_bI))/(IA_bI)^2)*(A_b_X)^2;
        
        Cnb = Cnb*Cbb;
        fCpitch(i,1) = asin(-Cnb(3,1));
        fCroll(i,1) = atan(Cnb(3,2)/Cnb(3,3));
        if Cnb(3,3) < 0
            if fCroll(i,1) > 0
                fCroll(i,1) = fCroll(i,1) - pi;
            else
                fCroll(i,1) = fCroll(i,1) + pi;
            end
        end
        fCyaw(i,1) = atan(Cnb(2,1)/Cnb(1,1));
        if Cnb(2,1) > 0
            if abs(Cnb(1,1)) < 0.0001
                fCyaw(i,1) = pi/2;
            elseif Cnb(1,1) < 0
                fCyaw(i,1) = fCyaw(i,1) + pi;
            end
        else
            if abs(Cnb(1,1)) < 0.0001
                fCyaw(i,1) = -pi/2;
            elseif Cnb(1,1) < 0
                fCyaw(i,1) = fCyaw(i,1) - pi;
            end
        end
        % yaw range: 0°-360°
        if fCyaw(i,1) < 0
            fCyaw(i,1) = fCyaw(i,1) + 2.*pi;
        end
    end
    % Vimu
    Cnb = [cos(fCpitch(i-1)).*cos(fCyaw(i-1)),sin(fCroll(i-1)).*sin(fCpitch(i-1)).*cos(fCyaw(i-1))-cos(fCroll(i-1)).*sin(fCyaw(i-1)),cos(fCroll(i-1)).*sin(fCpitch(i-1)).*cos(fCyaw(i-1))+sin(fCroll(i-1)).*sin(fCyaw(i-1));
           cos(fCpitch(i-1)).*sin(fCyaw(i-1)),sin(fCroll(i-1)).*sin(fCpitch(i-1)).*sin(fCyaw(i-1))+cos(fCroll(i-1)).*cos(fCyaw(i-1)),cos(fCroll(i-1)).*sin(fCpitch(i-1)).*sin(fCyaw(i-1))-sin(fCroll(i-1)).*cos(fCyaw(i-1));
           -sin(fCpitch(i-1)),sin(fCroll(i-1)).*cos(fCpitch(i-1)),cos(fCroll(i-1)).*cos(fCpitch(i-1))];
    V_m(1:3,1) = [mimu_vx(1);mimu_vy(1);mimu_vz(1)];
    if i == fh+1
        B_i2 = mimu_Lat(i-fh);
        vx_i2 = mimu_vx(i-fh);
        vy_i2 = mimu_vy(i-fh);
        vz_i2 = mimu_vz(i-fh);
        H_i2 = mimu_H(i-fh);
    else
        B_i2 = mimu_Lat(i-fh) + (mimu_Lat(i-fh) - mimu_Lat(i-fh-1))/2;
        vx_i2 =  V_m(1,i-fh) + ( V_m(1,i-fh) -  V_m(1,i-fh-1))/2;
        vy_i2 =  V_m(2,i-fh) + ( V_m(2,i-fh) -  V_m(2,i-fh-1))/2;
        vz_i2 =  V_m(3,i-fh) + ( V_m(3,i-fh) -  V_m(3,i-fh-1))/2;
        H_i2 = mimu_H(i-fh) + (mimu_H(i-fh) - mimu_H(i-fh-1))/2;
    end
    Rx = (re/sqrt(1-e11^2.*(sin(B_i2))^2));           % radius of curvature of prime unitary circle
    Ry = ((Rx.*(1-e11^2))/(1-e11^2.*(sin(B_i2))^2));  % radius of curvature of meridian circle
    Wn_ie = [Wie.*cos(B_i2);0;-Wie.*sin(B_i2)];
    Wn_en = [vy_i2/(Rx+H_i2);-vx_i2/(Ry+H_i2);-vy_i2.*tan(B_i2)/(Rx+H_i2)];
    gg = 9.7803267714.*(1+0.00527094.*(sin(B_i2))^2-0.0000232718.*(sin(2.*B_i2))^4) - 0.000003086.*H_i2;
    g_i2 = [0;0;gg];
    v_i2 = [vx_i2;vy_i2;vz_i2];
    Wn = 2.*Wn_ie+Wn_en;
    Wn_in = Wn_ie + Wn_en;
    Wn_in_X = [0 -Wn_in(3) Wn_in(2);
               Wn_in(3) 0 -Wn_in(1);
               -Wn_in(2) Wn_in(1) 0];
    dV_m = [m_ax(i);m_ay(i);m_az(i)];
    dAngle_m = [m_wx(i);m_wy(i);m_wz(i)];
    dV_rot = 0.5.*cross(dAngle_m,dV_m);
    dAngle_m1 = dAngle_m/2;
    dAngle_m2 = dAngle_m/2;
    dV_m1 = dV_m/2;
    dV_m2 = dV_m/2;
    dV_scul = (2/3).*(cross(dAngle_m1,dV_m2) + cross(dV_m1,dAngle_m2));
    
    dV_cor = (-cross(Wn,v_i2) + g_i2).*ts;
    dV_sf = (eye(3) - (ts/2).*Wn_in_X)*Cnb*(dV_m + dV_rot + dV_scul);
    V_m(1:3,i-fh+1) = V_m(1:3,i-fh) + dV_sf + dV_cor;
    % Pimu
    mimu_Lat(i-fh+1) = ((V_m(1,i-fh) + V_m(1,i-fh+1))/2).*ts/(Ry+H_i2)+mimu_Lat(i-fh);                  % latitude
    mimu_Lon(i-fh+1) = ((V_m(2,i-fh) + V_m(2,i-fh+1))/2).*ts/((Rx+H_i2).*cos(B_i2))+mimu_Lon(i-fh);     % longitude
    mimu_H(i-fh+1,1) = -((V_m(3,i-fh) + V_m(3,i-fh+1))/2).*ts + mimu_H(i-fh);                           % height
    Cnb = [cos(fCpitch(i)).*cos(fCyaw(i)),sin(fCroll(i)).*sin(fCpitch(i)).*cos(fCyaw(i))-cos(fCroll(i)).*sin(fCyaw(i)),cos(fCroll(i)).*sin(fCpitch(i)).*cos(fCyaw(i))+sin(fCroll(i)).*sin(fCyaw(i));
           cos(fCpitch(i)).*sin(fCyaw(i)),sin(fCroll(i)).*sin(fCpitch(i)).*sin(fCyaw(i))+cos(fCroll(i)).*cos(fCyaw(i)),cos(fCroll(i)).*sin(fCpitch(i)).*sin(fCyaw(i))-sin(fCroll(i)).*cos(fCyaw(i));
           -sin(fCpitch(i)),sin(fCroll(i)).*cos(fCpitch(i)),cos(fCroll(i)).*cos(fCpitch(i))];
    % ESKF
    if GNSS_H(i) ~= 0
        Rx = (re/sqrt(1-e11^2.*(sin(mimu_Lat(i-fh+1)))^2));           % radius of curvature of prime unitary circle
        Ry = ((Rx.*(1-e11^2))/(1-e11^2.*(sin(mimu_Lat(i-fh+1)))^2));  % radius of curvature of meridian circle
        Mpv = [1/(Ry+mimu_H(i-fh+1)) 0 0;
               0 1/((Rx+mimu_H(i-fh+1)).*cos(mimu_Lat(i-fh+1))) 0;
               0 0 -1];
        lb_P = Mpv*Cnb*lb; % lever-arm correction
        GNSS_Lat(i) = GNSS_Lat(i) - lb_P(1);
        GNSS_Lon(i) = GNSS_Lon(i) - lb_P(2);
        GNSS_H(i) = GNSS_H(i) - lb_P(3);
        lb_V = Cnb*cross(([m_wx(i);m_wy(i);m_wz(i)]/ts),lb); % lever-arm correction
        GNSS_vx(i) = GNSS_vx(i) - lb_V(1);
        GNSS_vy(i) = GNSS_vy(i) - lb_V(2);
        GNSS_vz(i) = GNSS_vz(i) - lb_V(3);
    end
    if GNSS_H(i) ~= 0
        eps1 = V_m(1,i-fh+1) - GNSS_vx(i);
        eps2 = V_m(2,i-fh+1) - GNSS_vy(i);
        eps3 = V_m(3,i-fh+1) - GNSS_vz(i);
        eps4 = mimu_Lat(i-fh+1) - GNSS_Lat(i);
        eps5 = mimu_Lon(i-fh+1) - GNSS_Lon(i);
        eps6 = mimu_H(i-fh+1) - GNSS_H(i);
        
        Cnb = [cos(fCpitch(i_kf)).*cos(fCyaw(i_kf)),sin(fCroll(i_kf)).*sin(fCpitch(i_kf)).*cos(fCyaw(i_kf))-cos(fCroll(i_kf)).*sin(fCyaw(i_kf)),cos(fCroll(i_kf)).*sin(fCpitch(i_kf)).*cos(fCyaw(i_kf))+sin(fCroll(i_kf)).*sin(fCyaw(i_kf));
               cos(fCpitch(i_kf)).*sin(fCyaw(i_kf)),sin(fCroll(i_kf)).*sin(fCpitch(i_kf)).*sin(fCyaw(i_kf))+cos(fCroll(i_kf)).*cos(fCyaw(i_kf)),cos(fCroll(i_kf)).*sin(fCpitch(i_kf)).*sin(fCyaw(i_kf))-sin(fCroll(i_kf)).*cos(fCyaw(i_kf));
               -sin(fCpitch(i_kf)),sin(fCroll(i_kf)).*cos(fCpitch(i_kf)),cos(fCroll(i_kf)).*cos(fCpitch(i_kf))];
        f_acc = [m_ax(i_kf);m_ay(i_kf);m_az(i_kf)]/ts;
        f_acc = Cnb*f_acc;
        Rx1 = (re/sqrt(1-e11^2.*(sin(GNSS_Lat(i_kf)))^2));
        Ry1 = ((Rx1.*(1-e11^2))/(1-e11^2.*(GNSS_Lat(i_kf))^2));
        
        Wnie = [Wie.*cos(GNSS_Lat(i_kf));
                0;
                -Wie.*sin(GNSS_Lat(i_kf))];
        Win = [(Wie.*cos(GNSS_Lat(i_kf))+GNSS_vy(i_kf))/(Rx1+GNSS_H(i_kf));
               -GNSS_vx(i_kf)/(Ry1+GNSS_H(i_kf));
               (-Wie.*sin(GNSS_Lat(i_kf))-GNSS_vy(i_kf).*tan(GNSS_Lat(i_kf)))/(Rx1+GNSS_H(i_kf))];
        Wnen = [GNSS_vy(i_kf)/(Rx1+GNSS_H(i_kf));
               -GNSS_vx(i_kf)/(Ry1+GNSS_H(i_kf));
               -GNSS_vy(i_kf).*tan(GNSS_Lat(i_kf))/(Rx1+GNSS_H(i_kf))];
        A_n = Win*T;
        A_n_X = [0 -A_n(3) A_n(2);
                 A_n(3) 0 -A_n(1);
                 -A_n(2) A_n(1) 0];
        Maa = -A_n_X;
        Mav = [0 1/(Rx1+GNSS_H(i_kf)) 0;
               -1/(Ry1+GNSS_H(i_kf)) 0 0;
               0 -tan(GNSS_Lat(i_kf))/(Rx1+GNSS_H(i_kf)) 0];
        M1 = [-Wie.*sin(GNSS_Lat(i_kf)) 0 0;
              0 0 0;
              -Wie.*cos(GNSS_Lat(i_kf)) 0 0];
        M2 = [0 0 -GNSS_vy(i_kf)/(Rx1+GNSS_H(i_kf))^2;
              0 0 GNSS_vx(i_kf)/(Ry1+GNSS_H(i_kf))^2;
              -GNSS_vy(i_kf)/((Rx1+GNSS_H(i_kf)).*cos(GNSS_Lat(i_kf))) 0 GNSS_vy(i_kf).*tan(GNSS_Lat(i_kf))/(Rx1+GNSS_H(i_kf))^2];
        Map = M1 + M2;
        Mva = [0 -f_acc(3) f_acc(2);
               f_acc(3) 0 -f_acc(1);
               -f_acc(2) f_acc(1) 0];
        Wie_ie_Win = 2.*Wnie + Wnen;
        Wie_ie_Win_X = [0 -Wie_ie_Win(3) Wie_ie_Win(2);
                        Wie_ie_Win(3) 0 -Wie_ie_Win(1);
                        -Wie_ie_Win(2) Wie_ie_Win(1) 0];
        Mvv = [0 -V_m(3,i_kf-fh+1) V_m(2,i_kf-fh+1);V_m(3,i_kf-fh+1) 0 -V_m(1,i_kf-fh+1);-V_m(2,i_kf-fh+1) V_m(1,i_kf-fh+1) 0]*Mav - Wie_ie_Win_X;
        Mvp = [0 -V_m(3,i_kf-fh+1) V_m(2,i_kf-fh+1);V_m(3,i_kf-fh+1) 0 -V_m(1,i_kf-fh+1);-V_m(2,i_kf-fh+1) V_m(1,i_kf-fh+1) 0]*(2*M1 + M2);
        Mpv = [1/(Ry+mimu_H(i_kf-fh+1)) 0 0;
               0 1/((Rx+mimu_H(i_kf-fh+1)).*cos(mimu_Lat(i_kf-fh+1))) 0;
               0 0 -1];
        Mpp = [0 0 -V_m(1,i_kf-fh+1)/(Ry+mimu_H(i_kf-fh+1))^2;
               V_m(2,i_kf-fh+1).*sin(mimu_Lat(i_kf-fh+1))/((Rx+mimu_H(i_kf-fh+1)).*(cos(mimu_Lat(i_kf-fh+1)))^2) 0 -V_m(2,i_kf-fh+1)/((Rx+mimu_H(i_kf-fh+1))^2.*cos(mimu_Lat(i_kf-fh+1)));
               0 0 0];
        F = [Maa Mav Map -Cnb zeros(3);
             Mva Mvv Mvp zeros(3) Cnb;
             zeros(3) Mpv Mpp zeros(3,6)
             zeros(6,15)];
        F_kf = eye(15) + F.*T;
        Z_kf = [eps1;eps2;eps3;eps4;eps5;eps6];
        H_kf = [zeros(3) eye(3) zeros(3) zeros(3,6);
                zeros(3) zeros(3) eye(3) zeros(3,6)];
        G_kf = [-Cnb zeros(3);
                zeros(3) Cnb;
                zeros(9,6)];
        Q_kf = (F_kf*G_kf*Q*G_kf'*F_kf').*T;
        
        X_kf = F_kf*X_kf;
        P_kf = F_kf*P_kf*F_kf'+Q_kf;
        
        K_kf = (P_kf*H_kf')*inv(H_kf*P_kf*H_kf'+R_kf);
        X_kf = X_kf+K_kf*(Z_kf-H_kf*X_kf);
        P_kf = (eye(15)-K_kf*H_kf)*P_kf;
        Qk(s,1) = sqrt(P11/P_kf(1,1));
        Qk(s,2) = sqrt(P22/P_kf(2,2));
        Qk(s,3) = sqrt(P33/P_kf(3,3));
        s = s + 1;
        
        Xkf_attitude = [Xkf_attitude,X_kf(1:3,1)];
        Xkf_velocity = [Xkf_velocity,X_kf(4:6,1)];
        Xkf_position = [Xkf_position,X_kf(7:9,1)];
        Xkf_gyroscope = [Xkf_gyroscope,X_kf(10:12,1)];
        Xkf_accelerometer = [Xkf_accelerometer,X_kf(13:15,1)];
        CNN = eye(3) - [0 -X_kf(3) X_kf(2);X_kf(3) 0 -X_kf(1);-X_kf(2) X_kf(1) 0];
        Cnb = [cos(fCpitch(i)).*cos(fCyaw(i)),sin(fCroll(i)).*sin(fCpitch(i)).*cos(fCyaw(i))-cos(fCroll(i)).*sin(fCyaw(i)),cos(fCroll(i)).*sin(fCpitch(i)).*cos(fCyaw(i))+sin(fCroll(i)).*sin(fCyaw(i));
               cos(fCpitch(i)).*sin(fCyaw(i)),sin(fCroll(i)).*sin(fCpitch(i)).*sin(fCyaw(i))+cos(fCroll(i)).*cos(fCyaw(i)),cos(fCroll(i)).*sin(fCpitch(i)).*sin(fCyaw(i))-sin(fCroll(i)).*cos(fCyaw(i));
               -sin(fCpitch(i)),sin(fCroll(i)).*cos(fCpitch(i)),cos(fCroll(i)).*cos(fCpitch(i))];
        Cnb = CNN'*Cnb;
        V_m(1,i-fh+1) = V_m(1,i-fh+1) - X_kf(4,1);
        V_m(2,i-fh+1) = V_m(2,i-fh+1) - X_kf(5,1);
        V_m(3,i-fh+1) = V_m(3,i-fh+1) - X_kf(6,1);
        mimu_Lat(i-fh+1) = mimu_Lat(i-fh+1) - X_kf(7,1);
        mimu_Lon(i-fh+1) = mimu_Lon(i-fh+1) - X_kf(8,1);
        mimu_H(i-fh+1) = mimu_H(i-fh+1) - X_kf(9,1);
        X_kf(1:9) = 0;
        i_kf = i;
    end
    fCpitch(i,1) = asin(-Cnb(3,1));
    fCroll(i,1) = atan(Cnb(3,2)/Cnb(3,3));
    if Cnb(3,3) < 0
        if fCroll(i,1) > 0
            fCroll(i,1) = fCroll(i,1) - pi;
        else
            fCroll(i,1) = fCroll(i,1) + pi;
        end
    end
    fCyaw(i,1) = atan(Cnb(2,1)/Cnb(1,1));
    if Cnb(2,1) > 0
        if abs(Cnb(1,1)) < 0.0001
            fCyaw(i,1) = pi/2;
        elseif Cnb(1,1) < 0
            fCyaw(i,1) = fCyaw(i,1) + pi;
        end
    else
        if abs(Cnb(1,1)) < 0.0001
            fCyaw(i,1) = -pi/2;
        elseif Cnb(1,1) < 0
            fCyaw(i,1) = fCyaw(i,1) - pi;
        end
    end
    % yaw range: 0°-360°
    if fCyaw(i,1) < 0
        fCyaw(i,1) = fCyaw(i,1) + 2.*pi;
    end
end

[~,lie] = size(Xkf_attitude);
x_time = 0:(lie-1);
figure
plot(x_time,Xkf_gyroscope(1,:).*(180/pi),'-m')
hold on
plot(x_time,Xkf_gyroscope(2,:).*(180/pi),'-b')
hold on
plot(x_time,Xkf_gyroscope(3,:).*(180/pi),'-r')
grid on
xlabel('Time/s');
ylabel('bias drift/°/s');
legend('X bias drift','Y bias drift','Z bias drift');
hold off

figure
subplot(311)
plot(Ti_mimu,fCroll.*(180/pi),'k',Ti_truth,truth_roll,'r');
title('roll');
xlabel('Time(s)');
ylabel('Angle(deg)');
grid on;
legend('EKF','truth')
subplot(312)
plot(Ti_mimu,fCpitch.*(180/pi),'k',Ti_truth,truth_pitch,'r');
title('pitch');
xlabel('Time(s)');
ylabel('Angle(deg)');
grid on;
legend('EKF','truth')
subplot(313)
plot(Ti_mimu,fCyaw.*(180/pi),'k',Ti_truth,truth_yaw,'r');
title('yaw');
xlabel('Time(s)');
ylabel('Angle(deg)');
grid on;
legend('EKF','truth')

for i = 1:hang1-fh+1
    % coordinate_IMU
    Rx = (re/sqrt(1-e11^2.*(sin(mimu_Lat(i)))^2));
    Ry = ((Rx.*(1-e11^2))/(1-e11^2.*(sin(mimu_Lat(i)))^2));
    n = e12.*cos(mimu_Lat(i));
    XX = c1.*(BB0.*mimu_Lat(i)+(BB2.*cos(mimu_Lat(i))+BB4.*(cos(mimu_Lat(i)))^3+BB6.*(cos(mimu_Lat(i)))^5+BB8.*(cos(mimu_Lat(i)))^7).*sin(mimu_Lat(i)));
    % calculate the central meridian（3°）
    NN = round((mimu_Lon(i).*(180/pi))/3);
    L0 = 3.*NN;
    l = ((mimu_Lon(i).*(180/pi)-L0).*60.*60)/Po;
    X_IMU(i,1) = XX+(1/2).*Rx.*tan(mimu_Lat(i)).*(cos(mimu_Lat(i)))^2.*(l)^2+(1/24).*Rx.*tan(mimu_Lat(i)).*(5-(tan(mimu_Lat(i)))^2+9.*(n^2)+4.*n^4).*(cos(mimu_Lat(i)))^4.*(l)^4+(1/720).*Rx.*tan(mimu_Lat(i)).*(61-58.*(tan(mimu_Lat(i)))^2+(tan(mimu_Lat(i)))^4).*(cos(mimu_Lat(i)))^6.*(l)^6;
    Y_IMU(i,1) = Rx.*cos(mimu_Lat(i)).*l+(1/6).*Rx.*(1-(tan(mimu_Lat(i)))^2+(n)^2).*(cos(mimu_Lat(i)))^3.*(l)^3+(1/120).*Rx.*(5-18.*(tan(mimu_Lat(i)))^2+(tan(mimu_Lat(i)))^4+14.*(n)^2-58.*(n)^2.*(tan(mimu_Lat(i)))^2).*(cos(mimu_Lat(i)))^5.*(l)^5+500000;
end
for i = 1:hang2
    % coordinate_truth
    Rx = (re/sqrt(1-e11^2.*(sin(truth_Lat(i)))^2));
    Ry = ((Rx.*(1-e11^2))/(1-e11^2.*(sin(truth_Lat(i)))^2));
    n = e12.*cos(truth_Lat(i));
    XX = c1.*(BB0.*truth_Lat(i)+(BB2.*cos(truth_Lat(i))+BB4.*(cos(truth_Lat(i)))^3+BB6.*(cos(truth_Lat(i)))^5+BB8.*(cos(truth_Lat(i)))^7).*sin(truth_Lat(i)));
    % calculate the central meridian（3°）
    NN = round((truth_Lon(i).*(180/pi))/3);
    L0 = 3.*NN;
    l = ((truth_Lon(i).*(180/pi)-L0).*60.*60)/Po;
    X_truth(i,1) = XX+(1/2).*Rx.*tan(truth_Lat(i)).*(cos(truth_Lat(i)))^2.*(l)^2+(1/24).*Rx.*tan(truth_Lat(i)).*(5-(tan(truth_Lat(i)))^2+9.*(n^2)+4.*(n)^4).*(cos(truth_Lat(i)))^4.*(l)^4+(1/720).*Rx.*tan(truth_Lat(i)).*(61-58.*(tan(truth_Lat(i)))^2+(tan(truth_Lat(i)))^4).*(cos(truth_Lat(i)))^6.*(l)^6;
    Y_truth(i,1) = Rx.*cos(truth_Lat(i)).*l+(1/6).*Rx.*(1-(tan(truth_Lat(i)))^2+(n)^2).*(cos(truth_Lat(i)))^3.*(l)^3+(1/120).*Rx.*(5-18.*(tan(truth_Lat(i)))^2+(tan(truth_Lat(i)))^4+14.*(n)^2-58.*(n)^2.*(tan(truth_Lat(i)))^2).*(cos(truth_Lat(i)))^5.*(l)^5+500000;
end
[qq,~] = find(GNSS_H == 0);
GNSS_Lat(qq,:) = [];
GNSS_Lon(qq,:) = [];
GNSS_H(qq,:) = [];
[hqq,~] = size(GNSS_H);
for i = 1:hqq
    Rx = (re/sqrt(1-e11^2.*(sin(GNSS_Lat(i)))^2));
    Ry = ((Rx.*(1-e11^2))/(1-e11^2.*(sin(GNSS_Lat(i)))^2));
    n = e12.*cos(GNSS_Lat(i));
    XX = c1.*(BB0.*GNSS_Lat(i)+(BB2.*cos(GNSS_Lat(i))+BB4.*(cos(GNSS_Lat(i)))^3+BB6.*(cos(GNSS_Lat(i)))^5+BB8.*(cos(GNSS_Lat(i)))^7).*sin(GNSS_Lat(i)));
    % calculate the central meridian（3°）
    NN = round((GNSS_Lon(i).*(180/pi))/3);
    L0 = 3.*NN;
    l = ((GNSS_Lon(i).*(180/pi)-L0).*60.*60)/Po;
    X_GNSS(i,1) = XX+(1/2).*Rx.*tan(GNSS_Lat(i)).*(cos(GNSS_Lat(i)))^2.*(l)^2+(1/24).*Rx.*tan(GNSS_Lat(i)).*(5-(tan(GNSS_Lat(i)))^2+9.*(n^2)+4.*(n)^4).*(cos(GNSS_Lat(i)))^4.*(l)^4+(1/720).*Rx.*tan(GNSS_Lat(i)).*(61-58.*(tan(GNSS_Lat(i)))^2+(tan(GNSS_Lat(i)))^4).*(cos(GNSS_Lat(i)))^6.*(l)^6;
    Y_GNSS(i,1) = Rx.*cos(GNSS_Lat(i)).*l+(1/6).*Rx.*(1-(tan(GNSS_Lat(i)))^2+(n)^2).*(cos(GNSS_Lat(i)))^3.*(l)^3+(1/120).*Rx.*(5-18.*(tan(GNSS_Lat(i)))^2+(tan(GNSS_Lat(i)))^4+14.*(n)^2-58.*(n)^2.*(tan(GNSS_Lat(i)))^2).*(cos(GNSS_Lat(i)))^5.*(l)^5+500000;
end

figure
plot3(X_IMU,Y_IMU,mimu_H,'-b')
hold on
plot3(X_truth,Y_truth,truth_H,'-r')
hold on
plot3(X_GNSS,Y_GNSS,GNSS_H,'-k')
xlabel('X/m');
ylabel('Y/m');
zlabel('Z/m');
legend('IMU trajectory','True trajectory','GNSS trajectory');
grid on
hold off
figure
plot(X_IMU,Y_IMU,'-b')
hold on
plot(X_truth,Y_truth,'-r')
hold on
plot(X_GNSS,Y_GNSS,'-k')
xlabel('X/m');
ylabel('Y/m');
legend('IMU trajectory','True trajectory','GNSS trajectory');
grid on
hold off
for i = 1:325838-fh+1
    VX1(i,1) = sqrt(V_m(1,i)^2 + V_m(2,i)^2 + V_m(3,i)^2);
end
for i = 1:325017
    VX2(i,1) = sqrt(truth_vx(i)^2 + truth_vy(i)^2 + truth_vz(i)^2);
end
figure
plot(Ti_mimu(fh:end),VX1,'-r')
hold on
plot(Ti_truth,VX2,'-k')
grid on
title('V');
xlabel('Time(s)');
ylabel('V(m/s)');
legend('EKF','truth')
hold off