
import numpy as np 
from matplotlib import pyplot as plt 
import os, sys
import os.path as osp
import math
import AutoinputQE as inpp
from matplotlib import font_manager
curPath = os.path.abspath(os.path.dirname(__file__))
sys.path.append(curPath)
fontpath = osp.join(curPath, "../font/arial.ttf")
# fontpath = "/work/wangr/data/dxy/scripts/program/QE/font/arial.ttf"
# ------------------ font setup ----------------------#
font_properties = font_manager.FontProperties(fname=fontpath)
styles = ['normal', 'italic', 'oblique']
weights = ['light', 'normal', 'medium', 'semibold', 'bold', 'heavy', 'black']
font = {'fontproperties': font_properties,  # 'family': 'Liberation Sans'
		'style': styles[0],
		'color': 'black',
		'weight': weights[1],
		'size': 22, 
        }
sys.dont_write_bytecode = True       ## not generate __pycache__ files

class Imag_phon():
  def __init__(self):
    self.filepath = os.getcwd()
    self.cmTomeV = 1/8.06554429
    self.lw = 1.0
    self.font_size = 20
    self.ImagpartName = self.Imag_partname()
    self.scale_ps = 20

  def Imag_partname(self):   # get real part file name
    files = os.listdir()
    ImagpartName = ''
    for i in range(len(files)):
      if 'linewidth.phself' in files[i]:
        ImagpartName = files[i]
        break
    return ImagpartName

  def write_Imag_Part(self, data):
    write_real = open(osp.join(self.filepath, 'Imag_weight.dat'), 'w')
    write_real.write('Imag Part (meV) \n')
    for i in range(len(data)):
      for j in range(len(data[i])):
        write_real.write(str(np.round(data[i, j], 10)).rjust(15, " ") + '  ')
      write_real.write('\n')
    write_real.close()

  def get_Imagpart(self, branches):     
    ImagPart = np.loadtxt(osp.join(self.filepath, self.ImagpartName), skiprows=2, dtype=float)
    tmp = ImagPart[:, 3].reshape(-1, branches)
    self.write_Imag_Part(tmp)
  
  def Imag_Weight_data(self):
    filepath = osp.join(self.filepath, '../phonon/pwscf.freq.gp')
    phonon = np.loadtxt(filepath, skiprows=0, dtype=float)
    branches = len(phonon[0, :]) - 1
    self.get_Imagpart(branches)

  #-------------------------  QE  -----------------------------
  def band_data_pro(self, filepath, filename):
      file1 = open(osp.join(filepath, filename), 'r').readlines()
      file2 = open(osp.join(filepath, filename + ".gp"), 'r').readlines()
      nbands = int(file1[0].split(',')[0].split('=')[1])
      nks = int(file1[0].split(',')[1].split("=")[1].split("/")[0])
      band_data = np.loadtxt(file2, skiprows=0, dtype=float)
      band_k = band_data[:nks, 0]
      kk = band_data[:, 0]
      energy = band_data[:, 1:len(band_data[0, :])]*self.cmTomeV
      hsp_labels = [str(items.split("!")[1]).strip() for items in inpp.kpoints_band.split("\n")[1:-1]]
      for i in range(len(hsp_labels)):
          if str(hsp_labels[i]).lower() == "\Gamma".lower():
              hsp_labels[i] = u"Γ"
      hsp_nums = len(hsp_labels)-1
      hsp_data = band_k[0:nks:(int(nks/(hsp_nums)))]
      hsp = {
          'klabels' : hsp_labels,
          'hspdata' : hsp_data,
      }
      # print(hsp['klabels'])
      return nbands, nks, hsp, band_k, kk, energy 
  
  def plot_ImagWeight_part(self):
    data_re = np.loadtxt(osp.join(self.filepath, 'Imag_weight.dat'), skiprows=1, dtype=float)
    nbands, nks, hsp, band_k, kk, energy  = self.band_data_pro(osp.join(self.filepath, '../phonon'), 'pwscf.freq')
    plt.figure(figsize=(5, 4))
    ax = plt.subplot(111)
    ########## QE-ph
    for kk in range(0, len(energy[0, :])):
      ax.plot(band_k, energy[:, kk], color='gray', lw=0.5)
      ax.scatter(band_k, energy[:, kk], s=data_re[:, kk]*self.scale_ps, color='red')
    # ax.scatter(band_k, energy[:, -1], s=data_re[:, -1]*self.scale_ps, color='red')
    max_value = np.max(energy)
    ############################
    for i in range(0,len(hsp['klabels'])):
        ax.vlines(hsp['hspdata'][i], ymin=0, ymax=max_value*2.5,linewidth=0.8,linestyle="--", color='gray')
    ax.axhline(0, xmin=0, xmax=hsp['hspdata'][-1],linewidth=1.0,linestyle="--", color='gray')
    ax.set_xticks(hsp['hspdata'])
    ax.set_xticklabels(hsp['klabels'])
    ax.set_ylabel("Frequency (meV)", fontdict=font)
    ax.set_ylim(np.min(energy), max_value*1.1)
    ax.set_xlim(0, np.max(band_k))
    labels = ax.get_xticklabels() + ax.get_yticklabels() 
    [label.set_fontproperties(font_properties) for label in labels]
    [label.set_fontsize(self.font_size-6) for label in labels]
    plt.savefig("Phonon_Weight.png", dpi=300, bbox_inches = 'tight', pad_inches=0.1)
    plt.close()

def manipulate():
  Ip = Imag_phon()
  Ip.Imag_Weight_data()
  Ip.plot_ImagWeight_part()

if __name__ == '__main__':
    manipulate()