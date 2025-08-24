from scipy.spatial import Voronoi
from itertools import product
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
import os 
import os.path as osp 
import numpy as np
from matplotlib import font_manager
import AutoinputQE as inpp
curPath = os.path.abspath(os.path.dirname(__file__))
fontpath = osp.join(curPath, "../font/arial.ttf")
#fontpath="D:\\JianGuoYun\\others\\keda\\old\\fermi\\font\\times.ttf"
# ------------------ font setup ----------------------#
font_properties = font_manager.FontProperties(fname=fontpath)
styles = ['normal', 'italic', 'oblique']
weights = ['light', 'normal', 'medium', 'semibold', 'bold', 'heavy', 'black']
font = {'fontproperties': font_properties,  # 'family': 'Liberation Sans'
		'style': styles[0],
		'color': 'black',
		'weight': weights[1],
		'size': 14, 
        }
# plt.rc('font',family=font, size=20)
##################################################### basic function ###################################################
# get reciprocal unit cell
def chaCheng(a1, a2):
    c = []
    c1 = float(a1[1] * a2[2] - a1[2] * a2[1])
    c.append(c1)
    c2 = float(a1[2] * a2[0] - a1[0] * a2[2])
    c.append(c2)
    c3 = float(a1[0] * a2[1] - a1[1] * a2[0])
    c.append(c3)
    return c

def dianCheng(b1, b2):
    d1 = float(b1[0] * b2[0])
    d2 = float(b1[1] * b2[1])
    d3 = float(b1[2] * b2[2])
    d = d1 + d2 + d3
    return d

def get_poscar(filepath=""):
    file = osp.join(filepath, "POSCAR")
    poscar = open(file, 'r')
    poscar_lines = poscar.readlines()
    pos = []
    for flag_lines in range(2, 5):
        for i in range(0, 3):
            pos.append(float(poscar_lines[flag_lines].split()[i]))
    return np.array(pos).reshape(3, 3)

def get_reciprocal(filepath=""):
    file = osp.join(filepath, "POSCAR")
    pos = get_poscar(filepath)
    rep = []
    a1 = [pos[0], pos[1], pos[2]]
    a2 = [pos[3], pos[4], pos[5]]
    a3 = [pos[6], pos[7], pos[8]]
    volume = dianCheng(a1, chaCheng(a2, a3))     # volume is a1 . (a2 X a3)
    scalar = 2 * np.pi / volume
    for i in range(0, 3):
        b1 = scalar*chaCheng(a2, a3)[i]
        rep.append(b1)
    for j in range(0, 3):
        b2 = scalar*chaCheng(a3, a1)[j]
        rep.append(b2)
    for k in range(0, 3):
        b3 = scalar * chaCheng(a1, a2)[k]
        rep.append(b3)
    return np.array(rep).reshape(3, 3)

def read_kpath(filepath=os.getcwd()):
    kpoints_bands = inpp.kpoints_band.split("\n")[1:-1]
    kpoints_data = []
    for i in range(len(kpoints_bands)):
        kpoints_data.append([float(item) for item in kpoints_bands[i].split('!')[0].split()])
    hsp_labels = [str(items.split("!")[1]).strip() for items in inpp.kpoints_band.split("\n")[1:-1]]
    for i in range(len(hsp_labels)):
        if str(hsp_labels[i]).lower() == "\Gamma".lower():
            hsp_labels[i] = u"Γ"
    print(np.array(kpoints_data))
    print(hsp_labels)
    return np.array(kpoints_data), hsp_labels

def read_poscar(poscar, species=None):
    poscar = open(poscar,'r')
    title = poscar.readline().strip()
    scale = float(poscar.readline().strip())
    s = float(scale)
    lattice_vectors = [[ float(v) for v in poscar.readline().split() ],
            [ float(v) for v in poscar.readline().split() ],
            [ float(v) for v in poscar.readline().split() ]]
    lattice_vectors = np.array(lattice_vectors)
    reciprocal_lattice_vectors= np.linalg.inv(lattice_vectors).T
    reciprocal_lattice_vectors=reciprocal_lattice_vectors*np.pi*2
    return reciprocal_lattice_vectors

def is_greek_alphabets(klabels):
    klabel = []
    for i in range(len(klabels)):
        klabel.append("$" + klabels[i] + "$")
    return klabel

def get_Wigner_Seitz_BZ(lattice_vectors):
# Inspired by http://www.thp.uni-koeln.de/trebst/Lectures/SolidState-2016/wigner_seitz_3d.py 
# Inspired by https://github.com/QijingZheng/VASP_FermiSurface/blob/master/fs.py 
    latt = []
    prefactors = [0., -1., 1.]
    for p in prefactors:
        for u in lattice_vectors:
            latt.append(p * u)
    lattice = []
    for vs in product(latt, latt, latt):
        a = vs[0] + vs[1] + vs[2]
        if not any((a == x).all() for x in lattice):
            lattice.append(a)
    voronoi = Voronoi(lattice)
    bz_facets = []
    bz_ridges = []
    bz_vertices = []
    for pid, rid in zip(voronoi.ridge_points, voronoi.ridge_vertices):
        if(pid[0] == 0 or pid[1] == 0):
            bz_ridges.append(voronoi.vertices[np.r_[rid, [rid[0]]]])
            bz_facets.append(voronoi.vertices[rid])
            bz_vertices += rid
    bz_vertices = list(set(bz_vertices))
    return voronoi.vertices[bz_vertices], bz_ridges, bz_facets

# class Arrow3D(FancyArrowPatch):
#     def __init__(self, xs, ys, zs, *args, **kwargs):
#         FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
#         self._verts3d = xs, ys, zs
#     def draw(self, renderer):
#         xs3d, ys3d, zs3d = self._verts3d
#         xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
#         self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
#         FancyArrowPatch.draw(self, renderer)

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        super().__init__((0, 0), (0, 0), *args, **kwargs)
        # FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs
    
    def do_3d_projection(self, renderer=None):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        return np.min(zs)


def visualize_BZ_matplotlib(points,ridges,facets,reciprocal_lattice_vectors,kpts,klabels):
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(6, 8))
    ax = fig.add_subplot(111, projection='3d', proj_type='ortho')
    basis_vector_clrs = ['k', 'k', 'k']
    basis_vector_labs = ['b$_1$', 'b$_2$', '$b_3$']
    for ii in range(3):
        arrow = Arrow3D([0,reciprocal_lattice_vectors[ii, 0]], [0,reciprocal_lattice_vectors[ii, 1]], [0,reciprocal_lattice_vectors[ii, 2]],
                color=basis_vector_clrs[ii], mutation_scale=10,lw=1,arrowstyle="->")
        ax.add_artist(arrow)
        ax.text(reciprocal_lattice_vectors[ii, 0], reciprocal_lattice_vectors[ii, 1],reciprocal_lattice_vectors[ii, 2],
                basis_vector_labs[ii], fontdict=font)
        for ir in ridges:
            ax.plot(ir[:, 0], ir[:, 1], ir[:, 2], color='gray', lw=1.0,alpha=0.5)
    for i in range(len(klabels)):
        kpt=np.dot(kpts[i,:], reciprocal_lattice_vectors)
        ax.scatter(kpt[0], kpt[1], kpt[2],c='b', marker='o',s=10,alpha=0.8)
        ax.text(kpt[0], kpt[1], kpt[2],klabels[i],c='red', fontdict=font)
    for i in range(kpts.shape[0]):
        kpts[i,:]=np.dot(kpts[i,:],reciprocal_lattice_vectors)
    for i in range(0,kpts.shape[0]-1,1):
        arrow = Arrow3D([kpts[i,0],kpts[i+1,0]],[kpts[i,1],kpts[i+1,1]],[kpts[i,2],kpts[i+1,2]],mutation_scale=10,lw=1.5,arrowstyle="->", color="green")
        ax.add_artist(arrow)
    ax.set_axis_off()
    ax.view_init(elev=10, azim=23)
    plt.savefig(osp.join(os.getcwd(), 'bz.png'), dpi=300, bbox_inches = 'tight', pad_inches=0.2)
    plt.close()

def visualize_BZ_pos(points,ridges,facets,reciprocal_lattice_vectors,kpts,klabels):
    import matplotlib.pyplot as plt
    pos = get_poscar(os.getcwd())
    fig = plt.figure(figsize=(6, 8))
    ax = fig.add_subplot(111, projection='3d', proj_type='ortho')
    basis_vector_clrs = ['blue', 'blue', 'blue']
    basis_vector_labs = ['b$_1$', 'b$_2$', '$b_3$']
    for ii in range(3):
        arrow = Arrow3D([0,reciprocal_lattice_vectors[ii, 0]], [0,reciprocal_lattice_vectors[ii, 1]], [0,reciprocal_lattice_vectors[ii, 2]],
                color=basis_vector_clrs[ii], mutation_scale=10,lw=1,arrowstyle="->")
        ax.add_artist(arrow)
        ax.text(reciprocal_lattice_vectors[ii, 0]*0.6, reciprocal_lattice_vectors[ii, 1]*0.6,reciprocal_lattice_vectors[ii, 2]*0.6,
                basis_vector_labs[ii], fontdict=font)
        for ir in ridges:
            ax.plot(ir[:, 0], ir[:, 1], ir[:, 2], color='gray', lw=1.0,alpha=0.5)
    basis_vector_pos_clrs = ['k', 'k', 'k']
    basis_vector_pos_labels = ['a$_1$', 'a$_2$', '$a_3$']
    for ii in range(3):
        arrow = Arrow3D([0,pos[ii, 0]], [0,pos[ii, 1]], [0,pos[ii, 2]],
                color=basis_vector_pos_clrs[ii], mutation_scale=10,lw=1,arrowstyle="->")
        ax.add_artist(arrow)
        ax.text(reciprocal_lattice_vectors[ii, 0]*0.4, reciprocal_lattice_vectors[ii, 1]*0.4,reciprocal_lattice_vectors[ii, 2]*0.4,
                basis_vector_pos_labels[ii], fontdict=font)

    for i in range(len(klabels)):
        kpt=np.dot(kpts[i,:], reciprocal_lattice_vectors)
        ax.scatter(kpt[0], kpt[1], kpt[2],c='b', marker='o',s=10,alpha=0.8)
        ax.text(kpt[0], kpt[1], kpt[2],klabels[i],c='red', fontdict=font)
    for i in range(kpts.shape[0]):
        kpts[i,:]=np.dot(kpts[i,:],reciprocal_lattice_vectors)
    for i in range(0,kpts.shape[0]-1,1):
        arrow = Arrow3D([kpts[i,0],kpts[i+1,0]],[kpts[i,1],kpts[i+1,1]],[kpts[i,2],kpts[i+1,2]],mutation_scale=10,lw=1.5,arrowstyle="->", color="green")
        ax.add_artist(arrow)
    ax.set_axis_off()
    ax.view_init(elev=10, azim=23)
    plt.savefig(osp.join(os.getcwd(), 'bz.png'), dpi=300, bbox_inches = 'tight', pad_inches=0.2)
    plt.close()
    
def manipulate():   
   reciprocal_lattice_vectors = read_poscar(osp.join(os.getcwd(), 'POSCAR'))   
   lattice_vectors = [np.array(reciprocal_lattice_vectors[0,:]),np.array(reciprocal_lattice_vectors[1,:]),np.array(reciprocal_lattice_vectors[2,:])]
   kpts, klabels = read_kpath()
   klabels = is_greek_alphabets(klabels)
   points, ridges, facets = get_Wigner_Seitz_BZ(lattice_vectors)
   visualize_BZ_matplotlib(points, ridges, facets, reciprocal_lattice_vectors,kpts,klabels)

def manipulate_lattice():   
   reciprocal_lattice_vectors = read_poscar(osp.join(os.getcwd(), 'POSCAR'))   
   lattice_vectors = [np.array(reciprocal_lattice_vectors[0,:]),np.array(reciprocal_lattice_vectors[1,:]),np.array(reciprocal_lattice_vectors[2,:])]
   kpts, klabels = read_kpath()
   klabels = is_greek_alphabets(klabels)
   print(klabels)
   points, ridges, facets = get_Wigner_Seitz_BZ(lattice_vectors)
   visualize_BZ_pos(points, ridges, facets, reciprocal_lattice_vectors,kpts,klabels)

if __name__ == "__main__":  
    manipulate()