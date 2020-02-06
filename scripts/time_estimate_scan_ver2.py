import numpy as np
import json

def SigmaPhotoelectric(Z, E):
    if E<50:
        E0 = 30
        A = (1.53 *np.exp(-0.0361*Z) +0.510 *np.exp(-0.176*Z)) *10**(-4)
        B = 0.0515*Z -0.18 *np.exp(-0.215*Z)
    else:
        E0 = 100
        A = (2.73 *np.exp(-0.02113*Z)
             +0.860 *np.exp(-0.128*Z)) / (1 -(2.63/Z +0.473) *(E -100) *10**(-3)) *10**(-6)
        B = 0.008*Z -0.18 *np.exp(-0.107*Z)
    return ( A * Z**5 * (E/E0)**(-3.32 +B) )* 10**(-24)

def SigmaCompton(Z, E):
    #if Z<4:
     #   A = Z/5
    #elif Z>=4:
    A = Z-3
    B = 5.5/(Z+5) + 1.06
    #B=1.5
    S = 1/ (1 + 5.184*10**(-3) *A**(0.8107) *(E/100)**(-B) )
    a = E/511
    sigma_KN = 0.49893 * (
        (1+a)/(a**2) * (
        2*(1+a)/(1+2*a) - 1/a * np.log(1+2*a)
        + 1/(2*a) * np.log(1+2*a)
        - (1+3*a)/((1+2*a)**2)
    )
    )
    return ( Z * sigma_KN * S )* 10**(-24)

def SigmaCoherent(Z,E):
    x = np.log(E/100)
    A = -6.6514 + 2.2802*np.log(Z) + 0.04174*(np.log(Z))**2
    B = 3 *10**(-4) *Z - 1.89
    C = -0.03
    D = -1.4 * np.exp( (3 * 10**(-4) *Z - 0.08)*E )
    return np.exp( A + B*x + C*x**2 + D ) *10**(-24)

def LinearAttenuation(sigma, A):
    NA = 6.02214076 * 10 ** (23)  # mol
    return sigma * NA / A

def Transmittance(mu, rho, d):
    I0 = 1 * 10 ** (6)  # Bq
    return I0 * np.exp(-mu*rho*d)


#f = open('time_estimate_result.json', 'w')

E = 60 # keV

Symbol = [ 'H',  'C',  'N',  'O', 'Al', 'Si', 'Cu']
Z      = [  1 ,   6 ,   7 ,   8 ,  13 ,  14 ,  29 ]
A      = [1.08,12.01,14.01,16.00,26.98,28.09,63.55] # g/mol
Z_dict = dict(zip(Symbol,Z))
A_dict = dict(zip(Symbol,A))

sigma_pe  = [SigmaPhotoelectric(Z[i],E) for i in range(len(Z))]
sigma_comp = [SigmaCompton(Z[i],E) for i in range(len(Z))]
sigma_cohe = [SigmaCoherent(Z[i],E) for i in range(len(Z))]

sigma_total = [sigma_pe[i]+sigma_comp[i]+sigma_cohe[i] for i in range(len(Z))]

keys = ['mu', 'I', 'd']

# single atom
rho_s = {'Al':2.7, 'Si':2.33, 'Cu':8.96} # g/cm3
d_s   = {'Al':0.1, 'Si':0.1 , 'Cu':0.002} # cm
for s in range(len(Z)):
    single_value = list()
    if Symbol[s] == 'Al' or Symbol[s] == 'Si' or Symbol[s] == 'Cu':
        mu_s = LinearAttenuation(sigma_total[s], A[s])
        single_value.append(mu_s)
        I_s = Transmittance(mu_s, rho_s[Symbol[s]], d_s[Symbol[s]])
        single_value.append(I_s)
        single_value.append(d_s[Symbol[s]])
        single = dict(zip(keys, single_value))
        print(Symbol[s])
        print(single)
        print(sigma_pe[s])
        print(sigma_comp[s])
        print(sigma_cohe[s])
        print(sigma_total[s])

# compound formula
Compound_symbol  = ['Polystyrene', 'ABS_resin', 'Polyester']
rho_c = [1.06, 1.08, 1.40] # g/cm3

Polystyrene_formula = {'H':8,'C':8} # C6H5CHCH2
ABS_formula = {'H':17,'C':15,'N':1} # C8H8C4H6C3H3N
Polyester_formula = {'H':8,'C':10,'O':4} # C10H8O4
#for key in Polyester_formula:
    #Polyester_mu +=
#json.dump(single, combound,indent(4))