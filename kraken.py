"""
Kraken python class

Features:
    -

Notes :
    - 

Developpements :
    - inputs for bathy/ssp are only float or list or ndarray
    -

Ronan Serre 
"""


import os, sys
import time
import numpy as np
import pathlib
from glob import glob
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator, MaxNLocator)


class Kraken:

    def __init__(self, job_title='kraken', at_path=os.getcwd(), warnings=False):
        """
        create paths and initialise dictionaries


        """

        self.job_title = job_title
        print(f'- Initialising {job_title}')

        # temporary directory for kraken run files
        # self.workdir = r'./tmp/'
        # self.workdir = os.path.abspath(self.workdir)
        # pathlib.Path(self.workdir).mkdir(parents=True, exist_ok=True)

        # results directory
        self.resultsdir = r'./results/'
        self.resultsdir = os.path.abspath(self.resultsdir)
        pathlib.Path(self.resultsdir).mkdir(parents=True, exist_ok=True)

        # program executable path [expected path towards acoustic toolbox main directory]
        self.prg = os.path.abspath(f'{at_path}/Kraken/krakenc.exe')
        self.field = os.path.abspath(f'{at_path}/KrakenField/field.exe')

        # check that executables are found
        if not os.path.exists(self.prg):
            print('Krakenc executable not found: check main acoustic toolbox path')
            sys.exit()

        if not os.path.exists(self.field):
            print('Field executable not found: check main acoustic toolbox path')
            sys.exit()

        # initialising dictionaries
        self.environment = {'range': None, 'bathymetry': None, 'depth': None, 'ssp': None, 'bottom': None,
                            'signal': None, 'fs': None}

        self.src = {'depth': None, 'frequency': None}

        self.rec = {'range': None, 'depth': None}

        self.computational_time = time.time()

        # Perhaps not a good idea, but if you really hate warnings
        if not warnings:
            import warnings
            warnings.filterwarnings('ignore')

        # input parameters that should be moved somewhere else later
        self.zpoints = 101
        self.max_number_of_modes = 101

        # computational marker
        self.compute = {"range_dependant_tl": False, "broadband": False}

        # initialise outputs
        self.output = {}


    def generate_env(self):
        """
        generate environement
        input range, bathymetry and sound speed profiles as constant (float), list or np.ndarray


        """

        # items in environment set to np.ndarray (except bottom as dict)
        # for key in list(self.environment.keys()):
        #     if key != "bottom":
        #         self.environment[key] = np.array(self.environment[key])

        # check ssp shape is (depth, range)
        if isinstance(self.environment["ssp"], np.ndarray):
            if self.environment["ssp"].ndim == 2:
                if self.environment["ssp"].shape[0] == self.environment["range"].shape[0] and self.environment["ssp"].shape[1] == self.environment["depth"].shape[0]:
                    self.environment["ssp"] = self.environment["ssp"].T

        # check if range dependent
        if isinstance(self.environment["bathymetry"], np.ndarray) or isinstance(self.environment["ssp"], np.ndarray):
            if self.environment["bathymetry"].shape[0] > 1 or self.environment["ssp"].shape[1] > 1:
                self.compute["range_dependant_tl"] = True

        if self.compute["range_dependant_tl"] and self.environment["bathymetry"].shape[0] == 1:
            if self.environment["ssp"].shape[1] != self.environment["range"].shape[0]:
                print('Make sure ssp.shape[1] == range.shape[0]')

            np.tile(self.environment["bathymetry"], self.environment["range"].shape[0])

        # check if broadband
        if self.environment["signal"] is not None and self.environment["fs"] is None:
            print('Broadband input signal is defined but not fs... ')
            sys.exit()

        elif self.environment["fs"] is not None and self.environment["signal"] is None:
            print('Fs is defined but not broadband input signal.. ')
            sys.exit()

        elif self.environment["fs"] is not None and self.environment["signal"] is not None:
            self.compute["broadband"] = True

        print('- Generating environment')

        self.source_data = {'n': len(self.src["depth"]), 'zs': self.src["depth"], 'f': self.src["frequency"]}

        # properties: not required due to vacuum over surface (empty list)
        # reflection: not required in this case (empty list)
        self.surface_data = {'bc': 'V', 'properties': [], 'reflection': []}

        # bump density in ridges/km , eta: principal radius 1 of bump, xi: principal radius 2 of bump
        self.scatter_data = {'bumden': [], 'eta': [], 'xi': []}

        # Number of mesh points (about 10 per vertical wavelength)
        self.nmesh = 0

        self.layert = 'HH'

        self.units = 'W'
        self.bc = 'A'

        # interfacial roughness
        self.sigma = 0.0

        # receiver(s) depth must be a list
        if isinstance(self.rec["depth"], np.ndarray):
            self.rec["depth"] = [i for i in self.rec["depth"]]


    def clean_files(self, files):
        # clean files before a new run

        for file in files:
            title_files = glob(f'{file}.*')
            for f in title_files:
                if os.path.exists(f):
                    os.remove(f)


    def run(self):

        if self.compute["range_dependant_tl"]:
            print('- Range-dependent calculation')
            self.range_dependant_wrapper()

        if self.compute['broadband']:
            print('- Broadband calculation')
            self.broadband_wrapper()

        self.computational_time = time.time() - self.computational_time
        print(f'- Computation time: {self.computational_time:.1f} s')


    def range_dependant_wrapper(self):
        """
        KRAKEN: range dependent calculations
        Mexilhoeira Grande, sex 03 abr 2026 15:12:21
        Written by Tordar

        It works like this:
        (1) load the transect
        (2) write a *.env for every pair (rbottom, zbottom) along the transect
        (3) merge all the *.env into a single "BIG" env
        (4) write the *.flp for RD calculations
        (5) get the field

        """

        self.field_data = {'rmax': 1000.0, 'rr': np.array([1]), 'nrr': 1, 'rp': 0, 'np': 1,
                           'm': 999, 'rmodes': 'A', 'stype': 'R', 'thorpe': 'T', 'finder': ' ',
                           'rd': self.rec["depth"], 'dr': np.zeros_like(self.rec["depth"]), 'nrd': len(self.rec["depth"])}

        for ienv in range(self.environment["range"].shape[0]):

            # new array on local bathymetry value
            local_depth_array = np.linspace(0.0, self.environment["bathymetry"][ienv], self.zpoints)

            # interpolate ssp array on local depth array
            if self.environment["ssp"].ndim == 1:
                local_ssp_array = np.tile(self.environment["ssp"], self.zpoints)

            else:
                local_ssp_array = np.interp(local_depth_array, self.environment["depth"], self.environment["ssp"][:, ienv])

            csw  = np.zeros_like(local_depth_array)
            rhow = np.ones_like(local_depth_array)
            apw  = np.zeros_like(local_depth_array)
            asw  = np.zeros_like(local_depth_array)

            self.ssp_data = {'type': 'H', 'itype': 'N', 'nmesh': self.nmesh, 'sigma': 0,
                             'clow': 0.0, 'chigh': 5000.0,
                             'cdata': np.vstack([local_depth_array, local_ssp_array, csw, rhow, apw, asw]),
                             'zbottom': local_depth_array[-1]}

            # [z, cp, cs, RHO, ap, as]
            layer_info_keys = ['depth', 'p_speed', 's_speed', 'density', 'absorption', 'as']
            layer_info = np.zeros((len(self.environment["bottom"]["thickness"]), len(layer_info_keys)), dtype=float)
            for i_layer in range(layer_info.shape[0]):
                depth = local_depth_array[-1] + self.environment["bottom"]["thickness"][i_layer]
                cp = self.environment["bottom"]["p_speed"][i_layer]
                rho = self.environment["bottom"]["density"][i_layer]
                ap = self.environment["bottom"]["absorption"][i_layer]
                layer_info[i_layer, :] = np.array([depth, cp, 0.0, rho, ap, 0.0])

            layer_prop = [self.nmesh]
            for i_layer in range(layer_info.shape[0]):
                layer_prop.append(local_depth_array[-1] + self.environment["bottom"]["thickness"][i_layer])
            layer_prop = np.array([layer_prop])

            properties = np.array(layer_info[-1])

            self.bottom_data = {"n": layer_info.shape[0], "layerp": layer_prop, "layert": self.layert,
                                "properties": properties, "bdata": np.array([layer_info]),
                                "units": self.units, "bc": self.bc, "sigma": self.sigma}

            self.wkrakenenvfil()

            os.system(f'mv tmp.env e{ienv:03d}.env')

        os.system(f'cat e???.env > {self.job_title}.env')
        os.system('rm e???.env')

        # write the flp file (clumsy, but functional):
        fid = open(f'{self.job_title}.flp','w')
        fid.write('Range-dependant calculations\n')
        fid.write('RA\n')
        fid.write(str(self.max_number_of_modes))
        fid.write('\n')
        fid.write(str(self.environment["range"].shape[0]))
        fid.write('\n')
        fid.write(str(0.001 * self.environment["range"][0]))
        fid.write(" ")
        fid.write(str(0.001 * self.environment["range"][-1]))
        fid.write(" /\n")
        fid.write(str(self.rec["range"].shape[0]))
        fid.write('\n')
        fid.write(str(0.001 * self.rec["range"][0]))
        fid.write(" ")
        fid.write(str(0.001 * self.rec["range"][-1]))
        fid.write(" /\n")
        fid.write(str(self.src["depth"].shape[0]))
        fid.write('\n')
        for js in range(self.src["depth"].shape[0]):
            fid.write(str(self.src["depth"][js]))
            fid.write(' /\n')
        fid.write(str(len(self.rec["depth"])))
        fid.write('\n')
        fid.write(str(self.rec["depth"][0]))
        fid.write(" ")
        fid.write(str(self.rec["depth"][-1]))
        fid.write(" /\n")
        # Yup, nza again... useful for non-vertical arrays
        fid.write(str(len(self.rec["depth"])))
        fid.write('\n')
        fid.write('0.0 0.0 /\n')
        fid.close( )

        os.system(f'{self.prg} {self.job_title}')
        os.system(f'{self.field} {self.job_title}')

        pressure, geometry = self.readshd()

        pressure = np.squeeze(pressure, axis=(0, 1))
        self.pressure = np.where(pressure == 0, np.nan, pressure)
        self.output["tl"] = - 20.0 * np.log10(abs(pressure)  + sys.float_info.epsilon)

        self.clean_files([self.job_title, 'field', 'tmp'])

    def broadband_wrapper(self):

        dr = np.zeros(self.environment["depth"].shape[0], dtype=float)

        self.field_data = {"rmax": 0.001 * self.environment["range"].max(),
                           "nrr": self.environment["range"].shape[0],
                           "rr": 0.001 * self.environment["range"],
                           "rp": 0, "np": 1, "m": 999,
                           "rmodes": 'A', "stype": 'R', "thorpe": 'T', "finder": ' ',
                           "rd": self.environment["depth"], "dr": dr, "nrd": self.environment["depth"].shape[0]}

        # [z, cp, cs, RHO, ap, as]
        layer_info_keys = ['depth', 'p_speed', 's_speed', 'density', 'absorption', 'as']
        layer_info = np.zeros((len(self.environment["bottom"]["thickness"]), len(layer_info_keys)), dtype=float)
        for i_layer in range(layer_info.shape[0]):
            depth = self.environment["depth"].max() + self.environment["bottom"]["thickness"][i_layer]
            cp = self.environment["bottom"]["p_speed"][i_layer]
            rho = self.environment["bottom"]["density"][i_layer]
            ap = self.environment["bottom"]["absorption"][i_layer]
            layer_info[i_layer, :] = np.array([depth, cp, 0.0, rho, ap, 0.0])

        layer_prop = [self.nmesh]
        for i_layer in range(layer_info.shape[0]):
            layer_prop.append(self.environment["depth"].max() + self.environment["bottom"]["thickness"][i_layer])
        layer_prop = np.array([layer_prop])

        properties = np.array(layer_info[-1])

        self.bottom_data = {"n": layer_info.shape[0], "layerp": layer_prop, "layert": self.layert,
                            "properties": properties, "bdata": np.array([layer_info]),
                            "units": self.units, "bc": self.bc, "sigma": self.sigma}

        zw = np.array([0.0, self.environment["depth"].max()])
        cw = np.array([self.environment["ssp"], self.environment["ssp"]])
        csw = np.array([0.0, 0.0])
        rhow = np.array([1.0, 1.0])
        apw = csw
        asw = csw

        self.ssp_data = {'type': 'H', 'itype': 'N', 'nmesh': self.nmesh, 'sigma': 0,
                         'clow': 0.0, 'chigh': 5000.0, 'cdata': np.array([zw, cw, csw, rhow, apw, asw]),
                         'zbottom': self.environment["depth"].max()}

        # FFT
        dt = 1.0 / self.environment["fs"]
        frq = np.fft.rfftfreq(self.environment["signal"].shape[0], d=dt)
        dftx = np.zeros(frq.shape[0], dtype=complex)

        for j in range(1, frq.shape[0]):

            print(f'   - frequency: {frq[j]:.1f} Hz')

            self.source_data["f"] = frq[j]

            self.wkrakenenvfil(case_title=self.job_title)

            os.system(f'{self.prg} {self.job_title}')
            os.system(f'cp field.flp {self.job_title}.flp')
            os.system(f'{self.field} {self.job_title} < {self.job_title}.flp')

            pressure, geometry = self.readshd()
            pressure = np.squeeze(pressure)

            zid = np.abs(self.environment["depth"] - self.rec["depth"]).argmin()
            rid = np.abs(self.environment["range"] - self.rec["range"]).argmin()

            dftx[j] = np.squeeze(pressure)[zid, rid]

            self.clean_files([self.job_title, 'field'])

        dfty = np.sqrt(2.0) * np.fft.rfft(self.environment["signal"]) / self.environment["signal"].shape[0]
        prod = dfty * dftx / (2.0 * np.pi)
        self.output["signal"] = np.real(2.0 * np.fft.irfft(prod) * prod.shape[0] / np.sqrt(2.0))


    def wkrakenenvfil(self, case_title='tmp'):

        """
        writes KRAKEN env file
        adapted from Orlando Camargo Rodriguez


        """

        env_file = f'{case_title}.env'
        trc_file = f'{case_title}.trc'
        fld_file = f'field.flp'

        # Get source data
        source_info = self.source_data
        freq = source_info["f"]
        zs = source_info["zs"]
        nzs = zs.shape[0]

        # Get surface data
        surface_info = self.surface_data
        top_boundary_condition = surface_info["bc"]
        top_properties = surface_info["properties"]
        top_reflection_coeff = surface_info["reflection"]

        # Get scatter data
        scatter_info = self.scatter_data
        bumden = scatter_info["bumden"]
        eta = scatter_info["eta"]
        xi = scatter_info["xi"]

        # Get sound speed data
        ssp_info = self.ssp_data
        ssp_data = ssp_info["cdata"]
        ssp_type = ssp_info["type"]
        citype = ssp_info["itype"]
        nmesh = ssp_info["nmesh"]
        csigma = ssp_info["sigma"]
        clow = ssp_info["clow"]
        chigh = ssp_info["chigh"]
        zbottom = ssp_info["zbottom"]

        # Get bottom data
        bottom_info = self.bottom_data
        nlayers = bottom_info["n"]
        attenuation_units = bottom_info["units"]
        bottom_boundary_condition = bottom_info["bc"]
        bottom_properties = bottom_info["properties"]
        bsigma = bottom_info["sigma"]
        layer_properties = bottom_info["layerp"]
        layer_type = bottom_info["layert"]
        layer_data = bottom_info["bdata"]

        if nlayers >= 20:
            sys.exit('Warning: max 20 layers...')

        # Get field data
        field_info = self.field_data
        thorpe = field_info["thorpe"]
        finder = field_info["finder"]
        rmax = field_info["rmax"]
        nrd = field_info["nrd"]
        nrr = field_info["nrr"]
        rd = field_info["rd"]
        rr = field_info["rr"]
        nmodes = field_info["m"]
        source_type = field_info["stype"]
        nprofiles = field_info["np"]
        rprofiles = field_info["rp"]
        dr = field_info["dr"]
        range_dependent_modes = field_info["rmodes"]

        # Construct the options
        options1 = citype + top_boundary_condition + attenuation_units + thorpe + finder
        options2 = source_type + range_dependent_modes

        # Write the ENV_FILE
        fid = open(env_file, 'w')
        fid.write('\'')
        fid.write(case_title)
        fid.write('\'\n')
        fid.write(str(freq))
        fid.write("\n")
        fid.write(str(nlayers))
        fid.write("\n")
        fid.write('\'')
        fid.write(options1)
        fid.write('\'\n')

        if top_boundary_condition == 'A':
            fid.write(str(top_properties[0]))
            fid.write(" ")
            fid.write(str(top_properties[1]))
            fid.write(" ")
            fid.write(str(top_properties[2]))
            fid.write(" ")
            fid.write(str(top_properties[3]))
            fid.write(" ")
            fid.write(str(top_properties[4]))
            fid.write(" ")
            fid.write(str(top_properties[5]))
            fid.write(" /")
            fid.write("\n")

        if top_boundary_condition == 'F':
            nthetas = surface_info["nthetas"]
            angle_data = surface_info["angle_data"]
            fidtrc = open(trc_file, 'w')
            fidtrc.write(str(nthetas))
            fid.write("\n")
            fidtrc.write(str(angle_data))
            fid.write("\n")
            fidtrc.close()

        if top_boundary_condition == 'F' or top_boundary_condition == 'I':
            fid.write(str(bumden))
            fid.write(" ")
            fid.write(str(eta))
            fid.write(" ")
            fid.write(str(xi))
            fid.write("\n")

        fid.write(str(nmesh))
        fid.write(" ")
        fid.write(str(csigma))
        fid.write(" ")
        fid.write(str(zbottom))
        fid.write("\n")

        if citype != 'A':
            nz = ssp_data[0,].size

            if ssp_type == 'H':
                fid.write(str(ssp_data[0, 0]))
                fid.write(" ")
                fid.write(str(ssp_data[1, 0]))
                fid.write(" ")
                fid.write(str(ssp_data[2, 0]))
                fid.write(" ")
                fid.write(str(ssp_data[3, 0]))
                fid.write(" ")
                fid.write(str(ssp_data[4, 0]))
                fid.write(" ")
                fid.write(str(ssp_data[5, 0]))
                fid.write(" /\n")

                for i in range(nz - 1):
                    fid.write(str(ssp_data[0, i + 1]))
                    fid.write(" ")
                    fid.write(str(ssp_data[1, i + 1]))
                    fid.write(" /")
                    fid.write("\n")

            else:
                for i in range(nz):
                    fid.write(str(ssp_data[0, i]))
                    fid.write(" ")
                    fid.write(str(ssp_data[1, i]))
                    fid.write(" ")
                    fid.write(str(ssp_data[2, i]))
                    fid.write(" ")
                    fid.write(str(ssp_data[3, i]))
                    fid.write(" ")
                    fid.write(str(ssp_data[4, i]))
                    fid.write(" ")
                    fid.write("\n")

        for i in range(nlayers - 1):
            fid.write(str(int(layer_properties[i, 0])))
            fid.write(" ")
            fid.write(str(layer_properties[i, 1]))
            fid.write(" ")
            fid.write(str(layer_properties[i, 2]))
            fid.write("\n")

            if layer_type[i] == 'H':
                fid.write(str(layer_data[i, 0, 0]))
                fid.write(" ")
                fid.write(str(layer_data[i, 0, 1]))
                fid.write(" ")
                fid.write(str(layer_data[i, 0, 2]))
                fid.write(" ")
                fid.write(str(layer_data[i, 0, 3]))
                fid.write(" ")
                fid.write(str(layer_data[i, 0, 4]))
                fid.write(" ")
                fid.write(str(layer_data[i, 0, 5]))
                fid.write("\n")
                fid.write(str(layer_data[i, 1, 0]))
                fid.write(" ")
                fid.write(str(layer_data[i, 1, 1]))
                fid.write(" /\n")

            else:
                fid.write(str(layer_data[i, 0, 0]))
                fid.write(" ")
                fid.write(str(layer_data[i, 0, 1]))
                fid.write(" ")
                fid.write(str(layer_data[i, 0, 2]))
                fid.write(" ")
                fid.write(str(layer_data[i, 0, 3]))
                fid.write(" ")
                fid.write(str(layer_data[i, 0, 4]))
                fid.write(" ")
                fid.write(str(layer_data[i, 0, 5]))
                fid.write("\n")
                fid.write(str(layer_data[i, 1, 0]))
                fid.write(" ")
                fid.write(str(layer_data[i, 1, 1]))
                fid.write(" ")
                fid.write(str(layer_data[i, 1, 2]))
                fid.write(" ")
                fid.write(str(layer_data[i, 1, 3]))
                fid.write(" ")
                fid.write(str(layer_data[i, 1, 4]))
                fid.write(" ")
                fid.write(str(layer_data[i, 1, 5]))
                fid.write("\n")

        fid.write("\'")
        fid.write(bottom_boundary_condition)
        fid.write("\' ")
        fid.write(str(bsigma))
        fid.write("\n")
        if bottom_boundary_condition == 'A':
            fid.write(str(bottom_properties[0]))
            fid.write(" ")
            fid.write(str(bottom_properties[1]))
            fid.write(" ")
            fid.write(str(bottom_properties[2]))
            fid.write(" ")
            fid.write(str(bottom_properties[3]))
            fid.write(" ")
            fid.write(str(bottom_properties[4]))
            fid.write(" ")
            fid.write(str(bottom_properties[5]))
            fid.write(" /")
            fid.write("\n")

        fid.write(str(clow))
        fid.write(" ")
        fid.write(str(chigh))
        fid.write("\n")
        fid.write(str(rmax))
        fid.write("\n")
        fid.write(str(nzs))
        fid.write("\n")

        if nzs == 1:
            fid.write(str(zs[0]))
            fid.write(" /\n")

        else:
            fid.write(str(zs[0]))
            fid.write(" ")
            fid.write(str(zs[-1]))
            fid.write(" /\n")

        fid.write(str(nrd))
        fid.write("\n")

        if nrd == 1:
            fid.write(str(rd[0]))
            fid.write("\n")

        else:
            fid.write(str(rd[0]))
            fid.write(" ")
            fid.write(str(rd[-1]))
            fid.write(" /\n")

        fid.close()

        # Write the FLD_FILE
        fid = open(fld_file, 'w')
        fid.write(case_title)
        fid.write("\n")
        fid.write(options2)
        fid.write("\n")
        fid.write(str(nmodes))
        fid.write("\n")
        fid.write(str(nprofiles))
        fid.write("\n")
        fid.write(str(rprofiles))
        fid.write("\n")
        fid.write(str(nrr))
        fid.write("\n")

        if nrr == 1:
            fid.write(str(rr[0]))
            fid.write("\n")

        else:
            fid.write(str(rr[0]))
            fid.write(" ")
            fid.write(str(rr[-1]))
            fid.write(" /\n")

        fid.write(str(nzs))
        fid.write("\n")

        if nzs == 1:
            fid.write(str(zs[0]))
            fid.write(" /\n")

        else:
            fid.write(str(zs[0]))
            fid.write(" ")
            fid.write(str(zs[-1]))
            fid.write(" /\n")

        fid.write(str(nrd))
        fid.write("\n")

        if nrd == 1:
            fid.write(str(rd[0]))
            fid.write(" /\n")

        else:
            fid.write(str(rd[0]))
            fid.write(" ")
            fid.write(str(rd[-1]))
            fid.write(" /\n")

        # Yes, this is ugly... but it works!!!

        fid.write(str(nrd))
        fid.write("\n")

        if nrd == 1:
            fid.write(str(dr[0]))
            fid.write(" /\n")

        else:
            fid.write(str(dr[0]))
            fid.write(" ")
            fid.write(str(dr[-1]))
            fid.write(" /\n")

        fid.close()


    def readshd(self, xs=np.nan, ys=np.nan, freq=np.nan):

        """
        based on read_shd_bin.m by Michael Porter
        adapted from Orlando Camargo Rodriguez


        """
        filename = f'{self.job_title}.shd'

        fid = open(filename, 'rb')
        recl = int(np.fromfile(fid, np.int32, 1))
        title = fid.read(80)
        fid.seek(4 * recl)
        PlotType = fid.read(10)

        # reposition to end of second record
        fid.seek(2 * 4 * recl)
        Nfreq = int(np.fromfile(fid, np.int32, 1))
        Ntheta = int(np.fromfile(fid, np.int32, 1))
        Nsx = int(np.fromfile(fid, np.int32, 1))
        Nsy = int(np.fromfile(fid, np.int32, 1))
        Nsz = int(np.fromfile(fid, np.int32, 1))
        Nrz = int(np.fromfile(fid, np.int32, 1))
        Nrr = int(np.fromfile(fid, np.int32, 1))
        atten = float(np.fromfile(fid, np.float32, 1))

        # reposition to end of record 3
        fid.seek(3 * 4 * recl)
        freqVec = np.fromfile(fid, np.float32, Nfreq)

        # reposition to end of record 4
        fid.seek(4 * 4 * recl)
        theta = np.fromfile(fid, np.float32, Ntheta)

        if PlotType[0:1] != 'TL':
            # reposition to end of record 4
            fid.seek(5 * 4 * recl)
            Xs = np.fromfile(fid, np.float32, Nsx)

            # reposition to end of record 5
            fid.seek(6 * 4 * recl)
            Ys = np.fromfile(fid, np.float32, Nsy)

        else:
            # compressed format for TL from FIELD3D
            # reposition to end of record 4
            fid.seek(5 * 4 * recl)
            Pos_S_x = np.fromfile(fid, np.float32, 2)
            Xs = np.linspace(Pos_S_x[0], Pos_S_x[1], Nsx)

            # reposition to end of record 5
            fid.seek(6 * 4 * recl)
            Pos_S_y = np.fromfile(fid, np.float32, 2)
            Ys = np.linspace(Pos_S_y[0], Pos_S_y[1], Nsy)

        # reposition to end of record 6
        fid.seek(7 * 4 * recl)
        zs = np.fromfile(fid, np.float32, Nsz)

        # reposition to end of record 7
        fid.seek(8 * 4 * recl)
        zarray = np.fromfile(fid, np.float32, Nrz)

        # reposition to end of record 8
        fid.seek(9 * 4 * recl)
        rarray = np.fromfile(fid, np.float64, Nrr)
        if PlotType == 'rectilin  ':
            pressure = np.zeros((Ntheta, Nsz, Nrz, Nrr), dtype=complex)
            Nrcvrs_per_range = Nrz

        elif PlotType == 'irregular ':
            pressure = np.zeros((Ntheta, Nsz, 1, Nrr), dtype=complex)
            Nrcvrs_per_range = 1

        else:
            pressure = np.zeros((Ntheta, Nsz, Nrz, Nrr), dtype=complex)
            Nrcvrs_per_range = Nrz

        if np.isnan(xs):
            ifreq = 0

            if np.isnan(freq) == False:
                freqdiff = np.abs(freqVec - freq)
                ifreq = freqdiff.argmin()

            for itheta in range(Ntheta):
                for isz in range(Nsz):
                    for irz in range(Nrcvrs_per_range):
                        recnum = (10
                                  + ifreq * Ntheta * Nsz * Nrcvrs_per_range
                                  + itheta * Nsz * Nrcvrs_per_range
                                  + isz * Nrcvrs_per_range
                                  + irz)

                        # Move to end of previous record
                        status = fid.seek(recnum * 4 * recl)
                        if status == -1:
                            print('Seek to specified record failed in readshd...')

                        # Read complex data
                        temp = np.fromfile(fid, np.float32, 2 * Nrr)
                        indexes = np.arange(0, 2 * Nrr, 2)
                        if temp.shape[0] < indexes.max():
                            pressure[itheta, isz, irz, :] = 0.0 + 1j * 0.0

                        else:
                            pressure[itheta, isz, irz, :] = temp[indexes] + 1j * temp[indexes + 1]

        else:
            xdiff = np.abs(Xs - xs * 1000.0)
            idxX = xdiff.argmin(0)
            ydiff = np.abs(Ys - ys * 1000.0)
            idxY = ydiff.argmin(0)

            for itheta in range(Ntheta):
                for isz in range(Nsz):
                    for irz in range(Nrcvrs_per_range):
                        recnum = (10
                                  + idxX * Nsy * Ntheta * Nsz * Nrcvrs_per_range
                                  + idxY * Ntheta * Nsz * Nrcvrs_per_range
                                  + itheta * Nsz * Nrcvrs_per_range
                                  + isz * Nrcvrs_per_range
                                  + irz)

                        # move to end of previous record
                        status = fid.seek(recnum * 4 * recl)
                        if status == -1:
                            print('Seek to specified record failed in read_shd_bin')

                        # read complex data
                        temp = np.fromfile(fid, np.float32, 2 * Nrr)
                        indexes = np.arange(0, 2 * Nrr, 2)
                        pressure[itheta, isz, irz, :] = temp[indexes] + 1j * temp[indexes + 1]

            fid.close()

        geometry = {"zs": zs, "f": freqVec, "thetas": theta, "rarray": rarray, "zarray": zarray}

        return pressure, geometry


    def generate_figures(self, vmin=60.0, vmax=180.0):

        print('- Generating figures')

        # BROADBAND
        if self.compute["broadband"]:

            py = self.environment["signal"]
            px = self.output["signal"]

            py /= np.abs(py).max()
            px /= np.abs(px).max()

            # time
            t = np.linspace(1.0, py.shape[0], py.shape[0]) / self.environment["fs"]

            # figure
            fig, axes = plt.subplots(1, 1, figsize=(6, 3), layout='constrained')
            axes.plot(t, py, '-', color='k', linewidth=1.5, label='signal')
            axes.plot(t, px, '-', color='tab:orange', linewidth=1.0, label='response')

            axes.legend(loc='upper right', frameon=True, framealpha=1.0, edgecolor='grey', fancybox=False,
                        prop={'size': 8}, ncol=1)

            axes.xaxis.set_tick_params(labelsize=8)
            axes.yaxis.set_tick_params(labelsize=8)
            axes.xaxis.set_minor_locator(AutoMinorLocator())
            axes.yaxis.set_minor_locator(AutoMinorLocator())

            axes.set_xlim(t.min(), t.max())
            axes.set_ylim(-1.1, 1.1)

            axes.set_xlabel(r'$time$ [s]', fontsize=8)
            axes.set_ylabel(r'$p(t)$ [-]', fontsize=8)

            axes.set_title(self.job_title, fontsize=8)

            plt.show()

        # TL
        if self.compute["range_dependant_tl"]:

            fig, axes = plt.subplots(1, 3, figsize=(8, 4),
                                    gridspec_kw={'width_ratios': [8, 1, 1]}, layout='constrained')

            # axis
            divisor = 1.0
            xlabel = r'Range (m)'
            xr = (self.rec["range"].min(), self.rec["range"].max())
            if xr[1] - xr[0] > 10000:
                divisor = 1000.0
                xlabel = r'Range (km)'

            extent = [self.rec["range"].min() / divisor, self.rec["range"].max() / divisor,
                      self.rec["depth"][0], self.rec["depth"][-1]]
            im = axes[0].imshow(np.flipud(self.output["tl"]), cmap='viridis', vmin=vmin, vmax=vmax,
                                extent=extent, aspect='auto', alpha=0.8)
            cbar = fig.colorbar(im, ax=axes[0], location='left', shrink=0.6, pad=0.01)
            cbar.ax.tick_params(labelsize=9)

            axes[0].plot(self.environment["range"] / divisor, self.environment["bathymetry"], color='peru', linewidth=1.0)

            txd = self.src["depth"]
            axes[0].plot([0] * np.size(txd), txd, marker='*', markerfacecolor='tab:green', markeredgecolor='k',
                         markeredgewidth=1.0, markersize=12)

            axes[0].set_xlim(0.0, self.environment["range"].max() / divisor)

            # axes[0].set_ylim(-surface[:, -1].max(), env['rx_depth'].max())

            if isinstance(self.environment["ssp"], float):
                axes[1].plot([self.environment["ssp"], self.environment["ssp"]], [0.0, self.environment["bathymetry"].max()], '-', color='k', linewidth=1)
                axes[2].plot([self.environment["ssp"], self.environment["ssp"]], [0.0, self.environment["bathymetry"].max()], '-', color='k', linewidth=1)

            else:
                if self.environment["ssp"].ndim == 1:
                    axes[1].plot(self.environment["ssp"], self.environment["depth"], '-', color='k', linewidth=1)
                    axes[2].plot(self.environment["ssp"], self.environment["depth"], '-', color='k', linewidth=1)

                elif self.environment["ssp"].ndim == 2:
                    for j in range(self.environment["range"].shape[0]):
                        axes[1].plot(self.environment["ssp"][:, j], self.environment["depth"], '-', color='k', linewidth=1)
                        axes[2].plot(self.environment["ssp"][:, j], self.environment["depth"], '-', color='k', linewidth=1)

            txd = self.src["depth"]
            axes[1].plot([0.0, 3000.0], [txd, txd], '--', color='r', linewidth=1.5)
            axes[2].plot([0.0, 3000.0], [txd, txd], '--', color='r', linewidth=1.5)

            if isinstance(self.environment["ssp"], float):
                axes[1].set_xlim(1490.0, 1510.0)
                axes[2].set_xlim(1490.0, 1510.0)
                axes[1].set_ylim(0.0, self.environment["bathymetry"].max())
                axes[2].set_ylim(0.0, 0.1 * self.environment["bathymetry"].max())


            else:
                axes[1].set_xlim(self.environment["ssp"].min(), self.environment["ssp"].max())
                axes[2].set_xlim(self.environment["ssp"][:np.abs(100 - self.environment["depth"]).argmin() + 1].min(),
                                 self.environment["ssp"][:np.abs(100 - self.environment["depth"]).argmin() + 1].max())

                axes[1].set_ylim(self.environment["depth"].min(), self.environment["depth"].max())
                axes[2].set_ylim(0.0, 100.0)

            axes[0].xaxis.set_minor_locator(AutoMinorLocator())
            for ax in axes:
                ax.invert_yaxis()
                ax.xaxis.set_tick_params(labelsize=8)
                ax.yaxis.set_tick_params(labelsize=8)
                ax.yaxis.set_minor_locator(AutoMinorLocator())
            axes[1].xaxis.set_tick_params(labelsize=6)
            axes[2].xaxis.set_tick_params(labelsize=6)
            axes[0].set_xlabel(xlabel, fontsize=9)
            axes[1].set_xlabel(r'Soundspeed [m/s]', fontsize=9)
            # axes[2].set_xlabel(r'Soundspeed [m/s]', fontsize=9)
            axes[0].set_ylabel(r'Depth [m]', fontsize=9)
            plt.show()
