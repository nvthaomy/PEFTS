import numpy as np
from collections import OrderedDict

class PolyFTS(object):
    def __init__(self, monomers_dict={}, interactions_dict={}, chains_dict={}, composition_dict={}):
        self.dim = 3
        self.cell_scaling = 1.
        self.cell_lengths = [20, 20, 20]
        self.cell_angles = [90, 90, 90]
        self.NPW = [128, 128, 128]

        self.monomers_dict = monomers_dict
        self.chains_dict = chains_dict
        self.interactions_dict = interactions_dict
        self.composition_dict = composition_dict

        self.read_input_fields = 'no'
        self.input_fields_file = 'fields.bin'

        self.ensemble = 'canonical'

        self.job_type = 'CL'
        self.field_updater = 'PO'
        self.dt = 0.01
        self.lambda_force_scale = []
        self.lambda_stress_scale = []
        self.num_time_steps_per_block = 1000
        self.num_blocks = 2000
        self.random_seed = 0
        self.scft_force_stopping_tol = 1e-7
        self.scft_stress_stopping_tol = 1e-7
        self.variable_cell = 'False'
        self.calc_density_operator = 'False'
        self.initfields = None
        self.DensityOutputByChain = 'False'
    def get_models(self):
        cell_dict = OrderedDict([['Dim', str(self.dim)], ['CellScaling', str(self.cell_scaling)],
                                 ['CellLengths', ' '.join(str(x) for x in self.cell_lengths)], 
                                 ['CellAngles', ' '.join(str(x) for x in self.cell_angles)],
                                 ['NPW', ' '.join(str(x) for x in self.NPW)] ])
        operators_dict = OrderedDict([["CalcHamiltonian", 'true'], ["CalcStressTensor", 'false'],
                                      ["CalcPressure", "true"], ["CalcChemicalPotential", "true"],
                                      ["CalcStructureFactor", "false"], ["CalcDensityOperator", self.calc_density_operator],
                                      ["IncludeIdealGasTerms", "true"]
                                    ])
        initfields_dict = OrderedDict([["ReadInputFields", self.read_input_fields], 
                                        ["InputFieldsFile", self.input_fields_file]])
        if self.initfields == None:
            for i in range(int(self.monomers_dict["NSpecies"].split()[0])):
                initfields_dict["initfield{}".format(i+1)] = {'inittype': 'urng'}
        else:
            for field_name, dict in self.initfields.items():
                initfields_dict[field_name] = dict
        model1_dict = OrderedDict([['cell', cell_dict], ['interactions', self.interactions_dict],
                                    ['composition', self.composition_dict], ['operators', operators_dict],
                                    ['initfields', initfields_dict]
        ])
        
        models_dict = OrderedDict([['NumModels', '1'], ['ModelType', 'MOLECULAR'], 
                                    ['monomers', self.monomers_dict], ['chains', self.chains_dict],
                                    ['model1', model1_dict]
        ])
                                    
        return models_dict
    
    def get_simulation(self):
        IO_dict = {"KeepDensityHistory": 'False', "KeepFieldHistory": "False",
                   "DensityOutputByChain": self.DensityOutputByChain, "OutputFormattedFields": "False",
                   "OutputFields": "HFields", "FieldOutputSpace": "both"}
        simulation_dict = OrderedDict([["JobType", self.job_type], ["FieldUpdater", self.field_updater],
                                        ["NumTimeStepsPerBlock", self.num_time_steps_per_block],
                                        ["NumBlocks", self.num_blocks], ["TimeStepDT", self.dt],
                                        ["RandomSeed", str(self.random_seed)], ["SCFTForceStoppingTol", str(self.scft_force_stopping_tol)],
                                        ["SCFTStressStoppingTol", str(self.scft_stress_stopping_tol)], 
                                        ["VariableCell", self.variable_cell], ["IO", IO_dict]])
        if self.lambda_force_scale:
            simulation_dict["LambdaForceScale"] = self.lambda_force_scale
        if self.lambda_stress_scale:
            simulation_dict["LambdaStressScale"] = self.lambda_stress_scale
        return simulation_dict

    def dict_to_str(self, _dict, dict_name, n_indent=0):
        s = '\n' + '  ' * n_indent + dict_name + " {"
        s += '\n'
        for key, val in _dict.items():
            if isinstance(val,dict):
                s_ = self.dict_to_str(val, key, n_indent=n_indent+1)
                s += s_
            else:
                if not isinstance(val,str):
                    if isinstance(val,bool) or isinstance(val,float) or isinstance(val,int):
                        val = str(val)
                    elif isinstance(val,list) or isinstance(val,tuple) or isinstance(val,np.ndarray):
                        val = " ".join([str(x) for x in val])
                s += '  ' * (n_indent+1) + key + '\t= ' + val + '\n'       
        s += '  ' * n_indent + "}\n"
        return(s)

    def write_PolyFTS(self):
        models_dict = self.get_models()
        simulation_dict = self.get_simulation()

        s = """#1) Nref = 1, Rg0 = 1nm = Rg
#2) bref = Rg sqrt(6/Nref) = 1nm sqrt(6/1) = sqrt(6)nm
#3) b = b_realUnit_fromSrel / bref
InputFileVersion = 3
"""
        s += self.dict_to_str(models_dict, 'models')
        s += '\n'
        s +=  self.dict_to_str(simulation_dict, 'simulation')
        s += '\n'
        s += """parallel {
  CUDA_selectdevice = 0
  CUDA_threadblocksize = 64
  OpenMP_nthreads = 4
}\n"""
        return(s)
