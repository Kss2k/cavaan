https://www.matrixcalculus.org/

sigma = 
G * inv(BStar) * GammaStar * Phi * GammaStar' * inv(BStar)' * G'
Mu = G * inv(BStar) * GammaStar * Tau

logLik += 
log(det(G * inv(BStar) * GammaStar * Phi * GammaStar' * inv(BStar)' * G')) + tr(S * inv(G * inv(BStar) * GammaStar * Phi * GammaStar' * inv(BStar)' * G')) - log(det(S)) - p +
      ((Nu - G * inv(BStar) * GammaStar * Tau)' * inv(G * inv(BStar) * GammaStar * Phi * GammaStar' * inv(BStar)' * G') * (Nu - G * inv(BStar) * GammaStar * Tau))

# Derivative with respect to Bstar
T_0 = np.linalg.inv(BStar)
T_1 = G * T_0
T_2 = ((((T_1 * GammaStar).dot(Phi)).dot(GammaStar.T)).dot(T_0.T)).dot(G.T)
T_3 = np.linalg.inv(T_2)
t_4 = (T_0).dot((GammaStar).dot(Tau))
t_5 = (Nu - (G).dot(t_4))
T_6 = (T_1.T).dot(T_3)
t_7 = (T_3).dot(t_5)
t_8 = (T_0.T).dot((G.T).dot(t_7))

functionValue = ((((np.log(np.linalg.det(T_2)) + np.trace((S).dot(T_3))) - np.log(np.linalg.det(S))) - p) + (t_5).dot(t_7))
gradient = ((((2 * ((((((((T_6).dot(S)).dot(T_3)).dot(G)).dot(T_0)).dot(GammaStar)).dot(Phi)).dot(GammaStar.T)).dot(T_0.T)) - (2 * ((((((T_6).dot(G)).dot(T_0)).dot(GammaStar)).dot(Phi)).dot(GammaStar.T)).dot(T_0.T))) + (2 * np.outer(t_8, t_4))) + (2 * np.outer(t_8, (((((((t_5).dot(T_3)).dot(G)).dot(T_0)).dot(GammaStar)).dot(Phi)).dot(GammaStar.T)).dot(T_0.T))))
# Derivative with respect to GammaStar

T0 = inv(BStar)
T1 = G * T0
T2 = inv(T1 * GammaStar * Phi * GammaStar.t() * T0.t() * G.t())
T3 = T1.t() * T2
t4 = Nu - G * T0 * GammaStar * Tau
t5 = T0.t() * G.t() * T2 * t4
