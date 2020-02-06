import math

def Sigma(Z, E):
    return 3 * 10**(12) * Z**4 / E**(3.5) # unit[barn]

def dEdx(Z, A, E):
    e = 1 # charge number of partcle
    me = 0.510998 # MeV/c^2
    gamma = E/me
    beta = math.sqrt(1 - 1/gamma**2)
    K = 0.307075 # cm2/mol
    I = 16 * Z**(0.9) # eV
    return - K * e**2 * Z/A * 1/(beta**2) * math.log( 2 * me * beta**2 * gamma**2 / I*10**(-6) ) # unit[MeVcm2/g]

def Probability(Z, A, rho, d, sigma):
    ne = rho * Z/A
    return ne * d * sigma

def Transmittance(mu, d):
    return math.exp(-mu * d)

# Source
E_Am = 60 # Am-241 energy
E_Sr = 2.280 # Sr-90 energy
Radiation = 1 *10**6 # Radioactivity unit[Bq]

# Atomic (H, C, N, O)
atomic = ['H', 'C', 'N', 'O']
Z_value = [1, 6, 7, 8]
A_value = [1.079, 12.0107, 14.0067, 15.9994]

# Material
Polyester_atomic_num = [4, 5, 0, 2] # (C10H8O4)n
Polyester_rho = 1.4 # unit[g/cm3]

ABS_atomic_num = [17, 15, 1, 0] # (C8H8·C4H6·C3H3N)n
ABS_rho = 1.08 # unit[g/cm3]

Polystyrene_atomic_num = [1, 1, 0, 0] # (C6H5CHCH2)n
Polystyrene_rho = 1.06 # unit[g/cm3]

Polyvinyltoluene_atomic_num = [1, 6, 0, 0] # (2-CH3C6H4CHCH2)n
Polyvinyltoluene_rho = 1.04 # unit[g/cm3]


Material = ['Sensor', 'Styrofoam', 'Airtight', 'Top Cover', 'Flex', 'Flex_Cu', 'Plastic Scintillator', 'Electric']
Z = [14,
     sum([Z_value[i]*Polystyrene_atomic_num[i] for i in range(4)])/sum([Polystyrene_atomic_num[i] for i in range(4)]),
     sum([Z_value[i]*ABS_atomic_num[i] for i in range(4)])/sum([ABS_atomic_num[i] for i in range(4)]),
     sum([Z_value[i]*ABS_atomic_num[i] for i in range(4)])/sum([ABS_atomic_num[i] for i in range(4)]),
     sum([Z_value[i]*Polyester_atomic_num[i] for i in range(4)])/sum([Polyester_atomic_num[i] for i in range(4)]),
     29,
     sum([Z_value[i]*Polyvinyltoluene_atomic_num[i] for i in range(4)])/sum([Polyvinyltoluene_atomic_num[i] for i in range(4)]),
     13] # Atomic number ( Number of electron in atomic )
A = [28.09,
     sum([A_value[j]*Polystyrene_atomic_num[j] for j in range(4)])/sum([Polystyrene_atomic_num[j] for j in range(4)]),
     sum([A_value[j]*ABS_atomic_num[j] for j in range(4)])/sum([ABS_atomic_num[j] for j in range(4)]),
     sum([A_value[j]*ABS_atomic_num[j] for j in range(4)])/sum([ABS_atomic_num[j] for j in range(4)]),
     sum([A_value[j]*Polyester_atomic_num[j] for j in range(4)])/sum([Polyester_atomic_num[j] for j in range(4)]),
     63.55,
     sum([A_value[j]*Polyvinyltoluene_atomic_num[j] for j in range(4)])/sum([Polyvinyltoluene_atomic_num[j] for j in range(4)]),
     26.98] # Atomic mass unit[g/mol]
rho = [2.33, Polystyrene_rho, ABS_rho, ABS_rho, Polyester_rho, 8.96, Polyvinyltoluene_rho,2.7] # unit[g/cm3]
d = [0.015, 5.0, 0.6, 0.1, 0.0180, 0.0020, 0.05, 0.1] # Thickness unit[cm]


dedx  = [dEdx(Z[i], A[i], E_Sr) for i in range(8)]
limit_thicks = [E_Sr / dedx[i] * rho[i] for i in range(8)]

print('For beta-ray, ')
for material, limit_thick in zip(Material, limit_thicks):
    print('  Material - ', material,' : Thicknesss limit - ', round(limit_thick,3), 'cm')


sigma = [Sigma(Z[i], (E_Am*10**3)) for i in range(8)]
mu = [rho[i] * Z[i] * sigma[i] / A[i] for i in range(8)]
P_xray, I_xray= [], []

for i in range(8):
    if i == 0:
        P_xray.append(Probability(Z[i], A[i], rho[i], d[i], sigma[i]))
        I_xray.append(1)
    else:
        P_xray.append(1)
        I_xray.append(Transmittance(mu[i], d[i]))

P_on_cover = 0.33 # toward sensor on top cover
P_on_elc   = 0.28
P_on_box   = 0.01

pixel = 400*192

HitRate = Radiation * P_on_cover * P_xray[0] * I_xray[3] * I_xray[4] * I_xray[5] / pixel
HitRate_elc = Radiation * P_on_elc * P_xray[0] * I_xray[3] * I_xray[4] * I_xray[5] * I_xray[7] / pixel

print('For X-ray, ')
print('  Hitrate by top cover : ', round(HitRate,3), 'Hz/pixel')
print('  Hitrate by electric components : ', round(HitRate_elc,3), 'Hz/pixel')

print('sigma')
print(sigma)
