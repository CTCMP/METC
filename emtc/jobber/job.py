import pathlib
from config import *


'''
    configure submit scripts
'''
submit_order = {
    'vc-relax'    : cal_order + "  " + "pw.x" + mpinpool + " < vc-relax.in > vc-relax.log 2>&1 \n",
    'relax'       : cal_order + "  " + "pw.x" + mpinpool + " < relax.in > relax.log 2>&1 \n",
    'scf'         : cal_order + "  " + "pw.x" + mpinpool + " < scf.in > scf.log 2>&1 \n",
    'nscf'        : cal_order + "  " + "pw.x" + mpinpool + " < nscf.in > nscf.log 2>&1 \n",
    'band'       : cal_order + "  " + "pw.x" + mpinpool + " < bands.in > bands.log 2>&1 \n",
    'pband'      : cal_order + "  " + "projwfc.x" + mpinpool + " < proBand.in > proBand.log 2>&1 \n",
    'band_plot'      : str(python_path) + ' -u ' + str(osp.join(scriptpath, 'scripts', 'band.py')) + ' \n',
    'dos'         : cal_order + "  " + "projwfc.x" + mpinpool + " < dos.in > dos.log 2>&1 \n",
    'dos_plot'       : str(python_path) + ' -u ' + str(osp.join(scriptpath, 'scripts', 'dos.py')) + ' \n',
    'phonon'         : cal_order + "  " + "ph.x" + mpinpool + " < phonon.in > ph.log 2>&1 \n",
    'phonon_TEST'    : cal_order + "  " + osp.join(ph_xTEST, 'test-suite','not_epw_comp', 'ph.x') + mpinpool + " < phononEPW.in > phEPW.log 2>&1 \n",
    'q2r'            : cal_order + "  " + "q2r.x" + mpinpool + " < q2r.in > q2r.log 2>&1 \n",
    'matdyn'         : cal_order + "  " + "matdyn.x" + mpinpool + " < matdyn.in > matdyn.log 2>&1 \n",
    'phonon_dos'        : cal_order + "  " + "matdyn.x" + mpinpool + " < phdos.in > phdos.log 2>&1 \n",
    'phonon_band_plot'    : str(python_path) + ' -u ' + str(osp.join(scriptpath, 'scripts', 'phonon_band.py')) + ' \n',
    'pw2wannier90'   : w90_post,
    'wannier90_pp'   : w90_pp,
    'wannier90'      : w90_distangle,
    'nihe_plot'      : str(python_path) + ' -u ' + str(osp.join(scriptpath, 'scripts', 'nihe.py')) + ' MLWFs \n',
    'epw'            : cal_order + " epw.x -npool " + str(inpp.num_process) + " < epw.in > epw.out 2>&1 \n",
    'phcollect_epw'  : str(python_path) + ' -u ' + str(osp.join(scriptpath, 'scripts', 'phCollect.py')) + ' \n',
    'last'           : 'date >> time \n',
    'ShengBTE'       : cal_order + " ShengBTE 2>BTE.err >BTE.out \n",
'perturbo_QE2pert'   : cal_order + " qe2pert.x " + mpinpool_must + "P -i QE2pert.in > QE2pert.out 2>&1 \n",
}

class jobSub:
    def __init__(self, caltype=None, header=None, order=None):
        self._header = header
        self._caltype = caltype
        self._order = order

    @property
    def header(self):
        return self._header
    @header.setter
    def header(self, header):
        self._header = header

    @property
    def caltype(self):
        return self._caltype
    @caltype.setter
    def caltype(self, caltype):
        self._caltype = caltype

    @property
    def order(self):
        return self._order
    @order.setter
    def order(self, order):
        self._order = order 
    
    def pack(self):
        pass 

    def writter_submit(self):
        pass


