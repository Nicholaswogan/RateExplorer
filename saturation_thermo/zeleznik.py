from scipy import constants as const
import jax
jax.config.update("jax_enable_x64", True)
from jax import grad
from jax import numpy as np

def compute_mu_eps(T):
    mu_111 = np.array([-0.235245033870e2, 0.406889449841e-1,  -0.151369362907e-4, 0.296144445015e4, 0.492476973663])
    mu_121 = np.array([0.111458541077e4, -0.118330789360e1, -0.209946114412e-2, -0.246749842271e6, 0.341234558134e2])
    mu_221 = np.array([-0.801488100747e2, -0.116246143257e-1, 0.606767928954e-5, 0.309272150882e4, 0.127601667471e2])
    mu_122 = np.array([0.888711613784e3, -0.250531359687e1, 0.605638824061e-3, -0.196985296431e6, 0.745500643380e2])
    # mu_jki = mu_kji
    mu_211 = mu_121
    mu_212 = mu_122
    mu_112 = np.zeros(mu_111.shape[0])
    mu_222 = np.zeros(mu_111.shape[0])
    
    eps_111 = np.array([0.288731663295e4, -0.332602457749e1, -0.282047283300e-2, -0.528216112353e6, 0.686997435643])
    eps_121 = np.array([-0.370944593249e3, -0.690310834523, 0.563455068422e-3, -0.382252997064e4, 0.942682037574e2])
    eps_211 = np.array([0.383025318809e2, -0.295997878789e-1, 0.120999746782e-4, -0.324697498999e4, -0.383566039532e1])
    eps_221 = np.array([0.232476399402e4, -0.141626921317, -0.626760562881e-2, -0.450590687961e6, -0.612339472744e2])
    eps_122 = np.array([-0.163385547832e4, -0.335344369968e1, 0.710978119903e-2, 0.198200003569e6, 0.246693619189e3])
    eps_212 = np.array([0.127375159848e4, 0.103333898148e1, 0.341400487633e-2, 0.195290667051e6, -0.431737442782e3])
    eps_112 = np.zeros(eps_111.shape[0])
    eps_222 = np.zeros(eps_111.shape[0])

    terms = np.array([1, T, T**2, 1/T, np.log(T)])
    
    mu = np.zeros((2,2,2))
    # mu[0,0,0] = np.sum(mu_111*terms)
    # mu[0,0,1] = np.sum(mu_112*terms)
    # mu[0,1,0] = np.sum(mu_121*terms)
    # mu[0,1,1] = np.sum(mu_122*terms)
    # mu[1,0,0] = np.sum(mu_211*terms)
    # mu[1,0,1] = np.sum(mu_212*terms)
    # mu[1,1,0] = np.sum(mu_221*terms)
    # mu[1,1,1] = np.sum(mu_222*terms)
    mu = mu.at[0,0,0].set(np.sum(mu_111*terms))
    mu = mu.at[0,0,1].set(np.sum(mu_112*terms))
    mu = mu.at[0,1,0].set(np.sum(mu_121*terms))
    mu = mu.at[0,1,1].set(np.sum(mu_122*terms))
    mu = mu.at[1,0,0].set(np.sum(mu_211*terms))
    mu = mu.at[1,0,1].set(np.sum(mu_212*terms))
    mu = mu.at[1,1,0].set(np.sum(mu_221*terms))
    mu = mu.at[1,1,1].set(np.sum(mu_222*terms))
    
    eps = np.zeros((2,2,2))
    # eps[0,0,0] = np.sum(eps_111*terms)
    # eps[0,0,1] = np.sum(eps_112*terms)
    # eps[0,1,0] = np.sum(eps_121*terms)
    # eps[0,1,1] = np.sum(eps_122*terms)
    # eps[1,0,0] = np.sum(eps_211*terms)
    # eps[1,0,1] = np.sum(eps_212*terms)
    # eps[1,1,0] = np.sum(eps_221*terms)
    # eps[1,1,1] = np.sum(eps_222*terms)
    eps = eps.at[0,0,0].set(np.sum(eps_111*terms))
    eps = eps.at[0,0,1].set(np.sum(eps_112*terms))
    eps = eps.at[0,1,0].set(np.sum(eps_121*terms))
    eps = eps.at[0,1,1].set(np.sum(eps_122*terms))
    eps = eps.at[1,0,0].set(np.sum(eps_211*terms))
    eps = eps.at[1,0,1].set(np.sum(eps_212*terms))
    eps = eps.at[1,1,0].set(np.sum(eps_221*terms))
    eps = eps.at[1,1,1].set(np.sum(eps_222*terms))
    
    return mu, eps

def minus_DG_RT(T, x):

    mu, eps = compute_mu_eps(T)

    # x = np.array([np.maximum(1e-100, x1), np.maximum(1e-100,1 - x1)])
    Phi = np.array([1.0, x[0]*x[1]])
    
    DG_RT = 0.0
    for i in range(2):
        tmp = 0.0
        for j in range(2):
            for k in range(2):
                tmp += (mu[j,k,i] + eps[j,k,i]*np.log(x[j]))*x[j]*x[k]
        DG_RT += Phi[i]*tmp

    return DG_RT

def gibbs_energy(T, x):
    DG_RT = minus_DG_RT(T, x)
    G = -const.R*T*DG_RT
    return G

def dgibbs_energy(T, x):
    def wrapper(x_):
        return gibbs_energy(T, x_)
    return grad(wrapper)(x)

def chemical_potentials(T, x):
    dG_dx = dgibbs_energy(T, x)
    G = gibbs_energy(T, x)
    mu1 = G - x[0]*dG_dx[0] - x[1]*dG_dx[1] + dG_dx[0]
    mu2 = G - x[0]*dG_dx[0] - x[1]*dG_dx[1] + dG_dx[1]
    mu = np.array([mu1, mu2])
    return mu

def saturation_vapor_pressures(T, x1, h_H2SO4, h_H2O):
    P_H2SO4_pure = h_H2SO4.sat_pressure(T)/1e6 # bars
    P_H2O_pure = h_H2O.sat_pressure(T)/1e6 # bars

    x = np.array([x1, 1 - x1])
    mus = chemical_potentials(T, x)
    
    lnP_H2SO4 = np.log(P_H2SO4_pure) + mus[0]/(const.R*T)
    P_H2SO4 = np.exp(lnP_H2SO4)

    lnP_H2O = np.log(P_H2O_pure) + mus[1]/(const.R*T)
    P_H2O = np.exp(lnP_H2O)
    
    return np.array([P_H2SO4, P_H2O])*1e6 # dynes/cm^2