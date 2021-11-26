function [Cz_equilibre, Delta_m_equilibre, alpha_equilibre_i, Cx_equilibre, F_equilibre] = point_equilibre(q_equilibre, S, m, g, Delta_m0_etoile, Cz_m, X, Y, Alpha0_2etoile, Cz_alpha, Cx_0, k, Epsilon)

% q_equilibre : dérivée assiette à l'équilibre
% S : surface alaire
% m : masse de l'avion
% g : constante gravitationelle
% Delta_m0 : braquage d'équilibre à portance nulle
% Cz_m : Variation de la portance en fonction de Delta_m
% X : distance entre le centre de gravité et le foyer de l'aile
% Y : distance entre le centre de gravité et le foyer de l'empennage
% Alpha0_2etoiles : incidence à portance et braquage nuls
% Cz_alpha : variation du coefficient de portance selon l'incidence
% Cx_0 : coefficient de trainée pour une portance nulle
% k : constante
% Epsilon : précision du résultat

alpha_equilibre_i = 0;
F_equilibre = 0;
alpha_equilibre_i1 = 3; 
i=0; 

while(abs(alpha_equilibre_i - alpha_equilibre_i1) > Epsilon)
    
    alpha_equilibre_i1 = alpha_equilibre_i;
    i=i+1;
    
    Cz_equilibre = 1/(q_equilibre*S)*(m*g-F_equilibre*sin(alpha_equilibre_i1));
    Delta_m_equilibre = Delta_m0_etoile-Cz_equilibre/Cz_m*X/(Y-X);
    alpha_equilibre_i = Alpha0_2etoile + Cz_equilibre/Cz_alpha-Cz_m/Cz_alpha*Delta_m_equilibre;
    Cx_equilibre = Cx_0+k*Cz_equilibre^2;
    F_equilibre = q_equilibre*S*Cx_equilibre/cos(alpha_equilibre_i);
    
end

fprintf(['Cz_equilibre : ',num2str(Cz_equilibre),'\n',...
         'Delta_m_equilibre : ',num2str(Delta_m_equilibre),'\n',...
         'Alpha_equilibre_i : ',num2str(alpha_equilibre_i),'\n',...
         'Cx_equilibre : ',num2str(Cx_equilibre),'\n',...
         'F_equilibre : ',num2str(F_equilibre),'\n',...
         'Nombre itérations : ',num2str(i),'\n'
        ])

end