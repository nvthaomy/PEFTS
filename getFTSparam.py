from scipy.integrate import simps
import numpy as np
import sys
sys.path.append('/home/mnguyen/bin/cgfts')
from cgfts.forcefield.forcefield_v2 import ForceField
from collections import OrderedDict
sys.path.append('/home/mnguyen/bin/PEFTS/')
from writePolyFTS import PolyFTS

ff_file = "/home/mnguyen/PE/PAA_PAH/CG/xp0.1_40AA24f1_40AH24f1_1300nacl_51880hoh/1.100fCG_100fAA/1.417nsMD_fixedHOH_aPE1.09nm_Cut17.6/FTS/PE_ff.dat"
xPAA = float(sys.argv[1])
fPAA = float(sys.argv[2])
NPAA = int(sys.argv[3])
NPAH = int(sys.argv[4])
csalt = float(sys.argv[5]) # in M
kT = 1.

output_file = 'input.in'
CG_sigma = 0.31 #nm
lb = 2.4 # sigma
bead_types = ['A-', 'B+', 'Na+', 'Cl-', 'HOH']
charges = [-1, 1, 1, -1, 0]
b_ref = np.sqrt(6)
chains = [['A-']*NPAA, ['B+']*NPAH, ['Na+'], ['Cl-'], ['HOH']]
chain_stat = 'DGC'
ensemble = 'canonical'

#====
precision = 8
vAA = 0.45**3 # nm^3
vAH = 0.45**3 # nm^3
v0 = 0.03 #water, na, cl
xPAH = (1-fPAA)/fPAA * xPAA
xOther = 1 - xPAA - xPAH
Ctot = 1/(vAA * xPAA + vAH * xPAH + v0 * xOther)
phi_solids = Ctot * (vAA * xPAA + vAH * xPAH)
print('phi_solids ',phi_solids)
xsalt = round(csalt*0.6022/Ctot,precision)
xNa = round(xPAA + xsalt,precision)
xCl = round(xPAH + xsalt,precision)
xW = round(xOther - xNa - xCl,precision)
vol_frac = [xPAA,xPAH,xNa,xCl,xW] # bead basis

# import sim ff file
ff = ForceField.from_sim_ff_file(ff_file, kT=kT)
ff.reorder_bead_types(bead_types)
print('Bead types from ff file: {}'.format(' '.join([bt.name for bt in ff.bead_types])))

def GetRMSBond(k,r0,nbin=1000,bmax=None):
    if bmax == None:
        if r0 > 0.:
            bmax = r0 * 10.
        else:
            bmax = 10.
    b, db = np.linspace(1e-4,bmax,num=nbin,retstep=True)

    U = k*(b-r0)**2
    num = np.multiply(b**4,np.exp(-U)) # a factor of b**2 comes from reexpressing in spherical coord.
    num = simps(num,b)
    den = np.multiply(b**2,np.exp(-U))
    den = simps(den,b)
    RMSb = np.sqrt(num/den)
    w = np.multiply(b**2,np.exp(-U))/den # 4pi cancel out

    # check if zero centered bond will give the same RMSb
    U0 = 3./2./RMSb**2 * b**2
    num0 = np.multiply(b**4,np.exp(-U0))
    num0 = simps(num0,b)
    den0 = np.multiply(b**2,np.exp(-U0))
    den0 = simps(den0,b)
    RMSb0 = np.sqrt(num0/den0)
    w0 = np.multiply(b**2,np.exp(-U0))/den0
    if np.abs(RMSb - RMSb0) > 0.01:
        raise Exception('RMS bond distance from offset Harmonic potential does not match RMS bond from resulting zero-centered Harmonic potential') 

    if r0 == 0: # initial potential is zero centered Harmonic bond
        RMSb = (3./2./k)**(0.5)
    return RMSb

# bond lengths
b = np.zeros((len(bead_types),len(bead_types)))
kuhn_lengths = np.ones(len(bead_types))
# excluded volume matrix
u0 = np.empty((len(bead_types), len(bead_types)))
interactions_dict = OrderedDict()
for i, bt1 in enumerate(bead_types):
    for j in range(i,len(bead_types)):
        bt2 = bead_types[j]
        try:
            r0 = ff.get_pair_potential("Bonded", bt1, bt2).Dist0.value
            k = ff.get_pair_potential("Bonded", bt1, bt2).FConst.value
            b[i,j] = b[j,i] = GetRMSBond(k, r0) * CG_sigma / b_ref
            print('Bond {} {}: r0 = {:.4f} bref'.format(bt1, bt2, b[i,j]))
            if i == j: # use like-bond if available
                kuhn_lengths[i] = str(b[i,j])
            else: # otherwise use the first bond involved this bead
                if kuhn_lengths[i] == 1.0:
                    kuhn_lengths[i] = str(b[i,j])
                if kuhn_lengths[j] == 1.0:
                    kuhn_lengths[j] = str(b[i,j])
        except:
            pass
        u0[i,j] = u0[j,i] = ff.get_pair_potential("Gaussian", bt1, bt2).excl_vol.value * CG_sigma**3
        interactions_dict["BExclVolume{}{}".format(i+1, j+1)] = str(u0[i,j])
        print('Gauss {} {}: u0 = {:.4f} kT nm^3'.format(bt1, bt2, u0[i,j]))

interactions_dict["EElecStatic"] = str(4 * np.pi * lb * CG_sigma)
interactions_dict["ApplyCompressibilityConstraint"] = "False"

# monomers_dict
monomers_dict = OrderedDict([["NSpecies", str(len(bead_types)) + ' # ' + " ".join(bead_types)], 
                              ["KuhnLen", " ".join([str(k) for k in kuhn_lengths])],
                              ["Charge", " ".join([str(c) for c in charges])], 
                              ["GaussSmearWidth", " ".join([str(bt.smear_length * CG_sigma) for bt in ff.bead_types])],
                            ])
# composition_dict
composition_dict = {"Ensemble": ensemble, "ChainVolFrac": [],
                    "CChainDensity": str(Ctot)}  

# chains_dict
chains_dict = OrderedDict([["NChains", str(len(chains))], ["polymerReferenceN", "1"]])
for i, chain in enumerate(chains):
    # condense block definitions
    block_species = []
    n_per_block = []
    for j,bead in enumerate(chain):
        if j == 0: 
            current_bead = bead
            cnt = 1
        else:
            if bead == current_bead:
                cnt += 1
            else:
                block_species.append(str(bead_types.index(current_bead) + 1))
                n_per_block.append(str(cnt))
                current_bead = bead
                cnt = 1
    block_species.append(str(bead_types.index(current_bead) + 1))
    n_per_block.append(str(cnt))
    chains_dict['chain{}'.format(i+1)] = {"Label": 'chain{}'.format(i+1), "NBlocks": str(len(n_per_block)),
                                          "NBeads": str(len(chain)), "Nperblock": n_per_block}
    if len(chain) > 1:
        chains_dict['chain{}'.format(i+1)].update({"BlockSpecies": block_species,
                                                   "Architecture": "linear", "Statistics": chain_stat})
    else:
        chains_dict['chain{}'.format(i+1)].update({ "Species": bead_types.index(chain[0]) + 1,
                                                   "Architecture": "point"})
    composition_dict["ChainVolFrac"].append(str(vol_frac[i]))
composition_dict["ChainVolFrac"]= " ".join(composition_dict["ChainVolFrac"])

fts = PolyFTS(monomers_dict=monomers_dict, interactions_dict=interactions_dict, 
              chains_dict=chains_dict, composition_dict=composition_dict)
fts.dim = 2
fts.NPW = [128,128]
fts.cell_scaling = 20
fts.cell_lengths = [1,1]
fts.cell_angles = [90]
fts.num_time_steps_per_block = 500
fts.num_blocks = 10000
fts.dt = 0.005
fts.job_type = 'SCFT'
fts.field_updater = 'SIS'
fts.lambda_force_scale = [2.5,0.25,1.75,2.5,4,5.0] # scale roughly ~ sqrt(eigenvalues), last field is electrostatic
fts.variable_cell = 'false'
fts.read_input_fields = 'no'
fts.input_fields_file = 'None'
fts.calc_density_operator = 'False'
s = fts.write_PolyFTS()

with open(output_file,'w') as f:
    f.write(s)
    f.close()
