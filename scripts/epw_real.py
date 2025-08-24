
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

class renormalization_phon():
  def __init__(self):
    self.filepath = os.getcwd()
    self.cmTomeV = 1/8.06554429
    self.lw = 1.0
    self.font_size = 20
    self.realpartName = self.realpartname()

  def realpartname(self):   # get real part file name
    files = os.listdir()
    realpartName = ''
    for i in range(len(files)):
      if 'realpart.phself' in files[i]:
        realpartName = files[i]
        break
    return realpartName
  
  def write_renormalized_realpart(self, qpoints, data):    # write realpart to renormalized_phonon.dat file.
    write_real = open(osp.join(self.filepath, 'renormalized_phonon.dat'), 'w')
    write_real.write('renormalized phonon frequency (meV) \n')
    for i in range(len(data)):
      write_real.write(str(np.round(qpoints[i], 10)).rjust(14, " ") + '  ')
      for j in range(len(data[i])):
        write_real.write(str(np.round(data[i, j], 8)).rjust(10, " ") + '  ')
      write_real.write('\n')
    write_real.close()
    
  def write_RealPart(self, data):
    write_real = open(osp.join(self.filepath, 'real.dat'), 'w')
    write_real.write('Real Part (meV) \n')
    for i in range(len(data)):
      for j in range(len(data[i])):
        write_real.write(str(np.round(data[i, j], 4)).rjust(10, " ") + '  ')
      write_real.write('\n')
    write_real.close()

  def get_realpart(self, branches):     
    realpart = np.loadtxt(osp.join(self.filepath, self.realpartName), skiprows=2, dtype=float)
    tmp = realpart[:, 3].reshape(-1, branches)
    self.write_RealPart(tmp)
    return tmp

  def get_renormalized_data(self):
    filepath = osp.join(self.filepath, '../phonon/pwscf.freq.gp')
    phonon = np.loadtxt(filepath, skiprows=0, dtype=float)
    qpoints = phonon[:, 0]
    branches = len(phonon[0, :]) - 1
    renormalized_data = np.zeros(shape=(len(phonon[:, 0]), branches))
    realpart = self.get_realpart(branches)
    for i in range(len(phonon[:, 0])):
      for j in range(branches):
        tmp = (phonon[i, j+1]*self.cmTomeV)**2 + 2*(phonon[i, j+1]*self.cmTomeV)*realpart[i, j]
        # renormalized_data[i, j] = np.sqrt(abs(tmp))
        if tmp >= 0:
          renormalized_data[i, j] = np.sqrt(tmp)
        else:
          renormalized_data[i, j] = -1 * np.sqrt(abs(tmp))
    # if osp.exists(osp.join(self.filepath, 'renormalized_phonon.dat')):
    #   print("Go on !")
    # else:
    self.write_renormalized_realpart(qpoints, renormalized_data)
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
  def plot_renormalizated_phonon(self):
    data_re = np.loadtxt(osp.join(self.filepath, 'renormalized_phonon.dat'), skiprows=1, dtype=float)
    nbands, nks, hsp, band_k, kk, energy  = self.band_data_pro(osp.join(self.filepath, '../phonon'), 'pwscf.freq')
    plt.figure(figsize=(5, 4))
    ax = plt.subplot(111)
    ########## QE-ph
    for kk in range(0, len(energy[0, :])-1):
      ax.plot(band_k, energy[:, kk], color='black', lw=self.lw)
    ax.plot(band_k, energy[:, -1], color='black', lw=self.lw, label='Ph')
    max_value = np.max(energy)
    ########## renormalized ph
    for i in range(len(data_re[0, :])-1):
      ax.scatter(data_re[:, 0], data_re[:, i+1], color='red', s=5, lw=self.lw)
    ax.scatter(data_re[:, 0], data_re[:, -1], color='red', s=5, label='Renormalized Ph', lw=self.lw)
    ############################
    for i in range(0,len(hsp['klabels'])):
        ax.vlines(hsp['hspdata'][i], ymin=0, ymax=max_value*2.5,linewidth=0.8,linestyle="--", color='gray')
    ax.axhline(0, xmin=0, xmax=hsp['hspdata'][-1],linewidth=1.0,linestyle="--", color='gray')
    ax.set_xticks(hsp['hspdata'])
    ax.set_xticklabels(hsp['klabels'])
    ax.set_ylabel("Frequency (meV)", fontdict=font)
    ax.set_ylim(0, max_value*1.1)
    ax.set_xlim(0, np.max(band_k))
    leg = ax.legend(loc="upper right")
    labels = ax.get_xticklabels() + ax.get_yticklabels() 
    [label.set_fontproperties(font_properties) for label in labels]
    [label.set_fontsize(self.font_size) for label in labels]
    labels1 = leg.get_texts() 
    [label.set_fontproperties(font_properties) for label in labels1]
    [label.set_fontsize(self.font_size-6) for label in labels1]
    plt.savefig("RePhon.png", dpi=300, bbox_inches = 'tight', pad_inches=0.1)
    plt.close()

  def test(self):
    data = np.loadtxt(osp.join(self.filepath, self.realpartName), skiprows=1, dtype=float)
    fig = plt.figure(figsize=(4, 4.5))
    ax = plt.subplot(111)
    ax.scatter(data[:, 2], data[:, 4], label='Imag Part', s=8, color='red')
    ax.scatter(data[:, 2], data[:, 5], label='Phonon Linewidth', s=5, color='blue')
    ax.set_xlabel("Frequency (meV)", fontdict=font)
    ax.set_ylabel("Scattering", fontdict=font)
    ax.set_ylim(0, max(np.max(data[:, 4]), np.max(data[:, 5]))*1.05)
    ax.set_xlim(np.min(data[:, 2]), np.max(data[:, 2]))
    leg = ax.legend(loc="upper right")
    labels = ax.get_xticklabels() + ax.get_yticklabels() 
    [label.set_fontproperties(font_properties) for label in labels]
    [label.set_fontsize(self.font_size-6) for label in labels]
    labels1 = leg.get_texts() 
    [label.set_fontproperties(font_properties) for label in labels1]
    [label.set_fontsize(self.font_size-6) for label in labels1]
    plt.savefig("TestImagPart.png", dpi=300, bbox_inches = 'tight', pad_inches=0.1)
    plt.close()

def manipulate():
  rp = renormalization_phon()
  rp.test()
  rp.get_renormalized_data()
  rp.plot_renormalizated_phonon()

if __name__ == '__main__':
    manipulate()