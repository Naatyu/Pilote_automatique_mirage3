%--------------------AUTEURS--------------------%
% Laoué Nathan, Cortes clément, Madec Alexandre % 
%-----------------------------------------------%

clc
clear 
close all

%---CONSTANTES---%
masse_aeronef = 8500; % Kg
Centrage_aeronef = 0.52; % pourcentage
S = 34; % m^2, surface alaire
Rayon_giration = 2.65; % metres
Longueur_reference = 5.24; % metres
Bi = 59691; % kg.m^2
R = 287.058; % Joules.Kg^-1.Kelvin^-1
Gamma = 1.4;
rho_0 = 1.225; % Kg.m-3
g = 9.81; % m.s^-2
P_0 = 101325; % Pa

%PDV 21%
Altitude = 20000*0.3048; % metres
Temperature = 288.15 - 0.0065*Altitude; % Kelvin
Mach = 0.8; 
P = P_0*(1+(-0.0065/288.15)*Altitude)^(-g/(R * -0.0065)); % Pa
rho = P/(R*Temperature); % Kg.m-3
Delta_m0_etoile = 0.022; % radians, Braquage d'équilibre à portance nulle
Cz_m = 1.1; 
Xm = 0.54; % metres, position du foyer en fonction du mach
Ym = 0.78; % metres, position du centre de poussée à partir du sommet représentant l'aile, rapportée à la corde de symétrie
L = 1.5*Longueur_reference; % metres, corde de symétrie
c = -Centrage_aeronef*L; % metres, distance entre le sommet et le centre de gravité
x = -Xm*L; % metres, distance entre le sommet et le foyer aerodynamique
y = -Ym*L; % metres, distance entre le sommet et le foyer des elevons
X = x-c; % metres, distance entre le centre de gravite et le foyer aerodynamique
Y = y-c; % metres, distance entre le centre de gravite et le foyer des elevons
Alpha0_2etoile = 0.02; % radians, incidence à portance et braquage nuls
Cz_alpha = 2.65; % variation du coefficient de portance selon l'incidence
Cx_0 = 0.015; % coefficient de trainée pour une portance nulle 
k = 0.22; 

%---CALCUL DU POINT D'EQUILIBRE---%
V_equilibre = Mach*sqrt(Gamma*R*Temperature);
q_equilibre = 1/2*(rho*V_equilibre^2);
Epsilon = 0.001;

[Cz_equilibre, Delta_m_equilibre, Alpha_equilibre_i, Cx_equilibre, F_equilibre] = point_equilibre(q_equilibre, S, masse_aeronef, g, Delta_m0_etoile, Cz_m, X, Y, Alpha0_2etoile, Cz_alpha, Cx_0, k, Epsilon);

%---CALCUL MODELE LONGITUDINAL SIMPLIFIE---%
Cx_alpha = 2*k*Cz_alpha*Cz_equilibre;
C_xm = 0; % car on néglige trainé braquage gouverne de profondeur 
Gamma_equilibre = 0; % car vol en palier
Ftau = 0; 
Cm_q = -0.68;
Cm_alpha = (Cz_alpha*X)/Longueur_reference;
Cmm = (Cz_m*Y)/Longueur_reference;

Xv = (2*q_equilibre*S*Cx_equilibre)/(masse_aeronef*V_equilibre);
X_alpha = (q_equilibre*S*Cx_alpha)/(masse_aeronef*V_equilibre); % simplifaction du sin par petits angles
X_gamma = (g*cos(Gamma_equilibre)/V_equilibre); 
Xm = (q_equilibre*S*C_xm)/(masse_aeronef*V_equilibre);
X_tau = -Ftau/(masse_aeronef*V_equilibre);

m_alpha = (q_equilibre*S*Longueur_reference*Cm_alpha)/Bi;
mq = (q_equilibre*S*Longueur_reference^2*Cm_q)/(Bi*V_equilibre);
mm = (q_equilibre*S*Longueur_reference*Cmm)/Bi;

Zv = 2*g/V_equilibre;
Z_alpha = (F_equilibre/(masse_aeronef*V_equilibre)+(q_equilibre*S*Cz_alpha)/(masse_aeronef*V_equilibre));
Z_gamma = 0; % car sinus(gamma)=0 car palier
Zm = (q_equilibre*S*Cz_m)/(masse_aeronef*V_equilibre);
Z_tau = 0; % car F_tau =0 et sin(alpha_equilibre)=0

A = [-Xv -X_gamma -X_alpha 0 0 0;
     Zv 0 Z_alpha 0 0 0;
     -Zv 0 -Z_alpha 1 0 0;
     0 0 m_alpha mq 0 0;
     0 0 0 1 0 0;
     0 V_equilibre 0 0 0 0];

B = [0 -X_tau;
     Zm Z_tau;
     -Zm -Z_tau;
     mm mm;
     0 0;
     0 0];
 
 disp(' ')
 disp('A :')
 disp(A)
 disp('B :')
 disp(B)
 
 A_red = [-Xv -X_gamma -X_alpha 0;
          Zv 0 Z_alpha 0;
         -Zv 0 -Z_alpha 1;
          0 0 m_alpha mq];
 
 B_red = [0;
          Zm;
         -Zm;
          mm;];
      
disp(' ')
disp('A_réduit :')
disp(A_red)
disp('B_réduit :')
disp(B_red)
 
damp(A_red)

rG = rank(ctrb(A_red,B_red))

A_simp = [-Xv -X_gamma 0 0;
          Zv 0 0 0;
          0 0 -Z_alpha 1;
          0 0 m_alpha mq];
      
B_simp = B_red;

C_simp = eye(4);

D_simp = [0 0 0 0]';

disp(' ')
disp('A_simplifié :')
disp(A_simp)
disp('B_simplifié :')
disp(B_simp)

damp(A_simp)

%---Oscillation rapide---%

A_i = [-Z_alpha 1;
       m_alpha mq];
   
B_i = [-Zm;
       mm];
   
C_i_alpha = [1 0];
C_i_q = [0 1];
   
disp(' ')
disp('A_incidence :')
disp(A_i)
disp('B_incidence :')
disp(B_i)   

damp(A_i)

FT_alpha_Dm_ss = ss(A_i,B_i,C_i_alpha,0)
FT_alpha_Dm = tf(FT_alpha_Dm_ss)
FT2_alpha_Dm = zpk(FT_alpha_Dm_ss)
dcgain(FT_alpha_Dm)

FT_q_Dm_ss = ss(A_i,B_i,C_i_q,0)
FT_q_Dm = tf(FT_q_Dm_ss)
FT2_q_Dm = zpk(FT_q_Dm_ss)
dcgain(FT_q_Dm)

figure(1)
step(FT_alpha_Dm,FT_q_Dm,10)
title('Réponse indicielle de \alpha/Dm et q/Dm')
legend('al/Dm','q/Dm')
grid

%---Oscillation lente--%

A_p = [-Xv -X_gamma;
       Zv 0];
   
B_p = [0;
       Zm];
   
C_p_vitesse = [1 0];   
C_p_gamma = [0 1];
   
disp(' ')
disp('A_phugoide :')
disp(A_p)
disp('B_phugoide :')
disp(B_p)   

damp(A_p)

FT_vitesse_Dm_ss = ss(A_p,B_p,C_p_vitesse,0)
FT_vitesse_Dm = tf(FT_vitesse_Dm_ss)
FT2_vitesse_Dm = zpk(FT_vitesse_Dm_ss)
dcgain(FT_vitesse_Dm)

FT_gamma_Dm_ss = ss(A_p,B_p,C_p_gamma,0)
FT_gamma_Dm = tf(FT_gamma_Dm_ss)
FT2_gamma_Dm = zpk(FT_gamma_Dm_ss)
dcgain(FT_gamma_Dm)

figure(2)
step(FT_vitesse_Dm,FT_gamma_Dm,1000)
grid
title('Réponse indicielle de v/Dm et \gamma/Dm')
legend('v/Dm','gam/Dm')

%---Etude avion simplifié---%

FT_avioncomplet_simplifie_Dm_ss = ss(A_simp,B_simp,C_simp,D_simp,'StateName',{'vitesse','gamma','alpha','vitesse_tangage'});
FT_avioncomplet_simplifie_Dm = tf(FT_avioncomplet_simplifie_Dm_ss);
set(FT_avioncomplet_simplifie_Dm, 'InputName','Delta_m','OutputName',{'vitesse','gamma','alpha','vitesse_tangage'})
FT_avioncomplet_simplifie_Dm
FT_avioncomplet_vitesse_simplifie_Dm = FT_avioncomplet_simplifie_Dm(1);
FT_avioncomplet_gamma_simplifie_Dm = FT_avioncomplet_simplifie_Dm(2);
FT_avioncomplet_alpha_simplifie_Dm = FT_avioncomplet_simplifie_Dm(3);
FT_avioncomplet_q_simplifie_Dm = FT_avioncomplet_simplifie_Dm(4);


figure(3)
pzplot(FT_avioncomplet_simplifie_Dm)
title('Poles des fonctions de transferts')
grid

figure(4)
step(FT_avioncomplet_simplifie_Dm,10)
title('Réponse indicielle de v/Dm, \gamma/Dm, \alpha/Dm et q/Dm sur 10 secondes')
grid
figure(5)
step(FT_avioncomplet_simplifie_Dm,400)
title('Réponse indicielle de v/Dm, \gamma/Dm, \alpha/Dm et q/Dm sur 400 secondes')
grid

%---simulation de l'avion complet---%

tsim=400;
sim('simulation_avion_complet_naturel')
figure(6)
plot(ans.temps,ans.vitesse,ans.temps,ans.gamma,ans.temps,ans.alpha,ans.temps,ans.q)
title('Résultats de la simulation sur 400 secondes')
grid
legend('v','gam','al','q')

tsim=10;
sim('simulation_avion_complet_naturel')
figure(7)
plot(ans.temps,ans.vitesse,ans.temps,ans.gamma,ans.temps,ans.alpha,ans.temps,ans.q)
title('Résultats de la simulation sur 10 secondes')
grid
legend('v','gam','al','q')

%---Determination de Kq---%

%sisotool('rlocus',FT_q_Dm) %Kq = 0.0785

%---Performances sur q---%

Kq=0.0785;
disp('Calcul de TqDm_bf par FT')
FT_q_Dm_bf=feedback(FT_q_Dm,-Kq)

damp(FT_q_Dm_bf)
figure(8)
step(FT_q_Dm,FT_q_Dm_bf);grid
title('Réponse Indicielle de q/Dm avec et sans Kq')
legend('Sans Kq','Avec Kq')

%---Performances sur alpha---%

A_i_bf = A_i+B_i*[0 Kq];
FT_alpha_Dm_bf_ss = ss(A_i_bf,B_i,C_i_alpha,0);
FT_alpha_Dm_bf = tf(FT_alpha_Dm_bf_ss)
damp(FT_alpha_Dm_bf)
figure(9)
step(FT_alpha_Dm,FT_alpha_Dm_bf);
grid
title('Réponse indicielle de al/Dm avec et sans Kq')
legend('Sans Kq','Avec Kq')

%---Amortisseur de tangage filtré---%

K = Z_alpha*mm-Zm*m_alpha;
tau_filtre_q = mm/K
FT_retour = -Kq*tf([tau_filtre_q 0], [tau_filtre_q 1]);
FT_q_Dm_filtre_bf = feedback(FT_q_Dm, FT_retour);
figure(10)
step(FT_q_Dm_filtre_bf,FT_q_Dm_bf,FT_q_Dm)
grid on
title('Effet du filtre sur q')
legend('Amortisseur filtré','Amortisseur non filtré','Sans amortisseur')

%---performances sur alpha avec filtre---%

FT_alpha_Dm_filtre_bf = feedback(1,FT_q_Dm*FT_retour)*FT_alpha_Dm;
figure(11)
step(FT_alpha_Dm_filtre_bf,FT_alpha_Dm_bf,FT_alpha_Dm)
grid
title('Effet du filtre sur \alpha')
legend('Amortisseur filtré','Amortisseur non filtré','Sans amortisseur')

%---insertion de la boucle de gouverne---%

Fn_BdG = 4;
Amor_BdG = sqrt(2)/2;
Wn_BdG = 2*pi*Fn_BdG;
[num_BdG,den_BdG] = ord2(Wn_BdG,Amor_BdG);
Kn_BdG = Wn_BdG^2*num_BdG;
FT_BdG = tf(Kn_BdG,den_BdG);
FT_cd=FT_BdG*FT_q_Dm;
FT_cr=tf([tau_filtre_q 0],[tau_filtre_q 1]);
%sisotool('rlocus',FT_cd*FT_cr) %Kq = 0.021

%---Experimentation---%

tsim = 6;
% Avion naturel
Kq=0;
SW_filtre = -1;
SW_BdG= -1;
sim('Experimentation')
figure(12)
plot(ans.temps,ans.alpha,'r',ans.temps,ans.q,'-.r')
title('Réponse indicielle de q/Dm et \alpha/Dm avec amortissage, filtre et BdG')
hold on

% Avion naturel + amortisseur
Kq=0.0785; 
SW_filtre = -1;
SW_BdG= -1;
sim('Experimentation')
plot(ans.temps,ans.alpha,'b',ans.temps,ans.q,'-.b')

% Avion naturel + amortisseur + Filtre
Kq=0.0785;
SW_filtre = 1;
SW_BdG= -1;
sim('Experimentation')
plot(ans.temps,ans.alpha,'k',ans.temps,ans.q,'-.k')

% Avion naturel + amortisseur + Filtre + BdG
Kq=0.0785;
SW_filtre = 1;
SW_BdG= 1;
sim('Experimentation')
plot(ans.temps,ans.alpha,'g',ans.temps,ans.q,'-.g')
grid; hold off;
legend('En rouge l''avion naturel (al continu; q pointillé)','En bleu l''avion naturel + amortisseur','En noir l''avion naturel + amortisseur + filtre','En vert l''avion naturel + amortisseur + filtre + BdG')

%---Tenue de pente---%

% simplification représentation d'état avec v=dv/dt=0 car automanette gere
% vitesse avion

A_AM = [0 0 0 0 0 0;
        0 0 Z_alpha 0 0 0;
        0 0 -Z_alpha 1 0 0;
        0 0 m_alpha mq 0 0;
        0 0 0 1 0 0;
        0 V_equilibre 0 0 0 0];
    
B_AM = [0;Zm;-Zm;mm;0;0];
C_AM = eye(6);
D_AM = zeros(6,1);

FT_AM_ss = ss(A_AM,B_AM,C_AM,D_AM);
FT_AM = tf(FT_AM_ss);

FT_AM_q_Dm = FT_AM(4)
FT_AM_gamma_Dm = FT_AM(2)

FT_gamma_Dmc = minreal(feedback(1,-Kq*FT_AM_q_Dm)*FT_AM_gamma_Dm)
%sisotool(FT_gamma_Dmc) %K_gamma = 0.315 (0.31498)
K_gamma = 0.315;
FT_gamma_bo = K_gamma*FT_gamma_Dmc
FT_gamma_bf = feedback(FT_gamma_bo,1,+1)
damp(FT_gamma_bf)

figure(13)
step(-FT_gamma_bf)
grid
title('Réponse indicielle de \gamma/\gamma_c')

%---Représentation d'état tenue de pente---%

A_gamma_bf = A_AM+B_AM*[0 K_gamma 0 Kq 0 0];
FT_gamma_bf_ss = ss(A_gamma_bf,-K_gamma*B_AM,C_AM,D_AM);
FT_gamma_bf = tf(FT_gamma_bf_ss)
damp(FT_gamma_bf(2))

figure(14)
step(0.017*FT_gamma_bf,10) %réponse de l'avion pour 1 degres (0.017 rad)
title("Réponse de l'avion pour \gamma_c = 1° pour v, \gamma, \alpha, q, \theta et z")
grid

%---Simulation tenue de pente---%

% Sans filtre et sans BdG
tsim = 10;
gamma_c = 1;
SW_BdG = -1;
SW_filtre = -1;

sim('A_tenue_pente_filtre_bdg')
figure(15)
plot(ans.temps,ans.Dm,ans.temps,ans.q,ans.temps,ans.alpha,ans.temps,ans.gamma)
title('Tenue de pente sans filtre et sans BdG')
legend('\bf\deltam','\bfq','\bf\alpha','\bf\gamma')
grid on

% Avec filtre et avec BdG

tsim = 10;
gamma_c = 1;
SW_BdG = 1;
SW_filtre = 1;

sim('A_tenue_pente_filtre_bdg')
figure(16)
plot(ans.temps,ans.Dm,ans.temps,ans.q,ans.temps,ans.alpha,ans.temps,ans.gamma)
title('Tenue de pente avec filtre et avec BdG')
legend('\bf\deltam','\bfq','\bf\alpha','\bf\gamma')
grid on

%---Limitation facteur de charge---%

Dn_max = 2;
alpha_MAX = Alpha_equilibre_i+(Alpha_equilibre_i-Alpha0_2etoile)*Dn_max;
epsi=1; 
ii=0; 
pas=0.005; 
tsim=10; 
SWBdG=1; 
SWF=1;
while (epsi>0.001);
gamma_c=pas*ii; 
ii=ii+1;
sim('A_tenue_pente_filtre_bdG')
alpha_max=max(ans.alpha); epsi=(alpha_MAX-alpha_max)

figure(17);
plot(ans.temps,ans.alpha/pi*180); grid on; hold on
end

disp(['gamma max en rad pour Delta nz = ',num2str(Dn_max)])
gamma_c=pas*(ii-1)
sim('A_tenue_pente_filtre_bdG')
plot(ans.temps,ans.alpha/pi*180,'r','LineWidth',2);grid on;hold off
title('\alpha en fonction de \gamma','FontSize',13)
text(3,6,['\bf\Deltan_m_a_x = ',...
num2str(Dn_max),' g'],'FontSize',12)
text(3,5,['\bf\alpha_m_a_x = ',...
num2str(alpha_max/pi*180),' °'],'FontSize',12)
text(3,4,['\bf\gamma_m_a_x = ',...
num2str(gamma_c/pi*180),' °'],'FontSize',12)

%---Tenue d'altitude amortie par la tenue de pente---%

%-méthode classique-%

FTz_gamma_bo = -FT_gamma_bf(6);
%sisotool(FTz_gamma_bo) %0.0016

Kz_gamma = 0.00167;
FTz_gamma_bo = Kz_gamma*FT_gamma_bf(6)
FTz_gamma_bf = minreal(feedback(FTz_gamma_bo,1))
allmargin(FTz_gamma_bo)
w6db = bandwidth(FTz_gamma_bf,-6)

figure(18)
step(FTz_gamma_bf,15)
grid
title('Réponse indicielle de z/zc pour Kz_\gamma=0.00167')

figure(19)
w=logspace(-1,1,500);
bodemag(FTz_gamma_bf,w)
grid
title('Représentation de bode pour Kz_\gamma=0.00167')

%-méthode du retour d'état-%

Az_gamma_bf = A_AM+B_AM*[0 K_gamma 0 Kq 0 Kz_gamma*K_gamma];
FTz_gamma_bf_ss = ss(Az_gamma_bf,-Kz_gamma*K_gamma*B_AM,C_AM,D_AM);
FTz_gamma_bf = tf(FTz_gamma_bf_ss);

figure(20);
step(FTz_gamma_bf(6),15)
grid
title('Réponse indicielle de z/zc pour Kz_\gamma=0.00167')

%-étude en simulation-%
% sans Bdg et sans filtre
z_c = 10;
tsim = 15;
SW_BdG = -1;
SW_filtre = -1;

figure(21)
sim('A_tenue_altitude_pente_filtre_BdG')
plot(ans.temps,ans.Dm,ans.temps,ans.gamma,ans.temps,ans.alpha,ans.temps,ans.q,ans.temps,ans.theta,ans.temps,ans.z/1000)
legend('Dm','\gamma','\alpha','q','\theta','z/1000')
title('Tenue de z sans BdG et sans filtre, zc=10m')
grid

% avec Bdg et avec filtre
z_c = 10;
tsim = 15;
SW_BdG = 1;
SW_filtre = 1;

figure(22)
sim('A_tenue_altitude_pente_filtre_BdG')
plot(ans.temps,ans.Dm,ans.temps,ans.gamma,ans.temps,ans.alpha,ans.temps,ans.q,ans.temps,ans.theta,ans.temps,ans.z/1000)
legend('Dm','\gamma','\alpha','q','\theta','z/1000')
title('Tenue de z avec BdG et avec filtre, zc=10m')
grid