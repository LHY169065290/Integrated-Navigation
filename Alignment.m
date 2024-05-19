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
g = [0;0;9.79361];
lb = [-0.073;0.302;0.087];

Ti_mimu = mimu_gnss(:,1);
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

mimu_Lat(1,1) = mimu_gnss(1,8).*(pi/180);
mimu_Lon(1,1) = mimu_gnss(1,9).*(pi/180);
mimu_H(1,1) = mimu_gnss(1,10);
mimu_vx(1,1) = truth_vx(1);
mimu_vy(1,1) = truth_vy(1);
mimu_vz(1,1) = truth_vz(1);
V_m = [];
X_GNSS(1,1) = 3371250.11913286;
Y_GNSS(1,1) = 545378.542243576;
X_IMU(1,1) = 3371250.11913286;
Y_IMU(1,1) = 545378.542243576;

BB0 = 1-(3/4).*e12^2+(45/64).*e12^4-(175/256).*e12^6+(11025/16384).*e12^8;
BB2 = BB0-1;
BB4 = (15/32).*e12^4-(175/384).*e12^6+(3675/8192).*e12^8;
BB6 = (-35/96).*e12^6+(735/2048).*e12^8;
BB8 = (315/1024).*e12^8;

%%
MisAngle_attitude = [0;0;0];
Mis_gyroscope = [0;0;0];
X_kf = [MisAngle_attitude;Mis_gyroscope];
P_kf = diag([(1.*(pi/180))^2 (1.*(pi/180))^2 (180.*(pi/180))^2 (25.*(pi/180)/3600)^2 (25.*(pi/180)/3600)^2 (25.*(pi/180)/3600)^2]);
Q = diag([(1*(2.909*1e-4))^2 (1*(2.909*1e-4))^2 (1*(2.909*1e-4))^2 0 0 0]);
R_kf = diag([(0.1.*pi/180)^2 (0.1.*pi/180)^2 (0.1.*pi/180)^2]);
Xkf_attitude = [];
Xkf_gyroscope = [];
counter = 0;

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

for i = 2:hang1
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
        counter = counter + 1;
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
        
        % vector increment upgrade alpha
        alpha = Cbb_gnss*(v1 + v2 + 0.5.*cross((seita1 + seita2),(v1 + v2)) + (2/3).*(cross(seita1,v2) + cross(v1,seita2)));
        Tx = (v1 + v2 + 0.5.*cross((seita1 + seita2),(v1 + v2)) + (2/3).*(cross(seita1,v2) + cross(v1,seita2)));
        T_x = [0 -Tx(3) Tx(2);
               Tx(3) 0 -Tx(1);
               -Tx(2) Tx(1) 0];
        Am1 = cross((((T/2).*eye(3) + (T^2/6).*Win_X)*Wnie),vn);
        vn0 = vn;
        vn = [GNSS_vx(i);GNSS_vy(i);GNSS_vz(i)];
        Am2 = cross((((T/2).*eye(3) + (T^2/3).*Win_X)*Wnie),vn);
        beta1 = Cnn'*(Am1 + Am2 - (T.*eye(3) + (T^2/2).*Win_X)*g);
        % upgrade Cnn
        Cnnm = Cnn; % vector increment
        Cnn = Cnnt*Cnn;
        % vector increment upgrade beta
        beta = Cnn'*vn - Cnnm'*vn0 + beta1;
        
        A = A + beta*alpha';
        B = A'*A;
        [Eigen_vector,D] = eig(B);
        D = sqrt(D);
        V = Eigen_vector;
        U = A*V*inv(D);
        Cnb0 = U*V';
        
        % AEKF
        if counter == 1
            wib_X = zeros(3);
        else
            wib_X = [0 -wibb(3) wibb(2);
                     wibb(3) 0 -wibb(1);
                     -wibb(2) wibb(1) 0];
        end
        F = [-wib_X eye(3)
             zeros(3) zeros(3)];
        F_kf = eye(6) + F*T;
        Z_kf = Cnb0'*beta - alpha;
        H_kf = [Cbb_gnss*T_x zeros(3)];
        Q_kf = Q.*ts;
        X_kf = F_kf*X_kf;
        P_kf = F_kf*P_kf*F_kf'+Q_kf;
        
        K_kf = (P_kf*H_kf')*inv(H_kf*P_kf*H_kf'+R_kf);
        X_kf = X_kf+K_kf*(Z_kf-H_kf*X_kf);
        P_kf = (eye(6)-K_kf*H_kf)*P_kf;
        
        Xkf_attitude = [Xkf_attitude,X_kf(1:3,1)];
        Xkf_gyroscope = [Xkf_gyroscope,X_kf(4:6,1)];
        CBB = eye(3) - [0 -X_kf(3) X_kf(2);X_kf(3) 0 -X_kf(1);-X_kf(2) X_kf(1) 0];
        Cbb_gnss = Cbb_gnss*CBB;
        X_kf(1:3) = 0;
        
        Cbb_gnss = Cbb_gnss*Cbb;
        Cbb = eye(3);
        wibb = seita1 + seita2;
        v1 = 0;v2 = 0;
        seita1 = 0;seita2 = 0;
        
        Cnb = Cnn*Cnb0*Cbb_gnss;
        
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

figure
subplot(311)
plot(Ti_mimu,fCroll.*(180/pi),'k',Ti_truth,truth_roll,'r');
title('roll');
xlabel('Time(s)');
ylabel('Angle(deg)');
grid on;
legend('AEKF','truth')
subplot(312)
plot(Ti_mimu,fCpitch.*(180/pi),'k',Ti_truth,truth_pitch,'r');
title('pitch');
xlabel('Time(s)');
ylabel('Angle(deg)');
grid on;
legend('AEKF','truth')
subplot(313)
plot(Ti_mimu,fCyaw.*(180/pi),'k',Ti_truth,truth_yaw,'r');
title('yaw');
xlabel('Time(s)');
ylabel('Angle(deg)');
grid on;
legend('AEKF','truth')