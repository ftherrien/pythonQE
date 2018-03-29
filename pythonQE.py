import re
import numpy as np
import glob
import os
import uuid
from subprocess import call
from copy import deepcopy

class pwcalc:
    def __init__(self):
        self.type = 'pw.x'
        self.name = "si"
        self.calc_type = "vc-relax"
        self.restart_mode = "from_scratch"
        self.pseudo_dir = os.path.expanduser("~/soft/qe-6.2.1/pseudo/pseudo_fam_pbe-kjpaw-psl/")
        self.celldm = 10.263
        self.ecutwfc = 45.0
        self.ecutrho = 400.0
        self.nbnd = 8
        self.occupations = "fixed"
        self.masses = {'Si':28.0855}
        self.cell = [[0.5, 0.5, 0],[0.5, 0, 0.5],[0, 0.5, 0.5]]
        self.atomic_pos = {'Si':[[0,0,0],[0.25,0.25,0.25]]}
        self.kpoints = [8,8,8]


    def write_in(self):
        with open('/u/cb/cr/felixtherrien/scratch/pythonQE/canvas/pw_QEcanvas.in') as f:
            s = f.read()
        s = re.sub("{calculation type}", self.calc_type, s)
        s = re.sub("{restart mode}", self.restart_mode, s)
        s = re.sub("{name}", self.name, s)
        s = re.sub("{pseudo dir}", self.pseudo_dir, s)
        s = re.sub("{celldm}", str(self.celldm), s)
        s = re.sub("{ecutwfc}", str(self.ecutwfc), s)
        s = re.sub("{ecutrho}", str(self.ecutrho), s)
        s = re.sub("{nbnd}", str(self.nbnd), s)
        s = re.sub("{occupations}", self.occupations, s)
        
        atomic_pos = ""
        atomic_spec = ""
        ntyp = 0
        nat = 0 
        wmass = 0
        for atom_type in self.atomic_pos:
            ntyp += 1
            if self.masses.has_key(atom_type):
                pseudoname = os.path.basename(glob.glob(self.pseudo_dir + '/' + atom_type + '*')[0]) 
                atomic_spec += (atom_type + ' ' + str(self.masses[atom_type]) + ' ' +
                                pseudoname +
                                '\n')
            else:
                raise AttributeError("Masses has no key %s"%atom_type)
            for atom in self.atomic_pos[atom_type]:
                wmass += self.masses[atom_type]
                nat += 1
                atomic_pos += (atom_type + ' ' +
                               ' '.join([str(i) for i in atom]) + '\n')
                
        cell = ""
        for line in self.cell:
            cell += ' '.join([str(i) for i in line]) + '\n'
        
        s = re.sub("{nat}", str(nat), s)
        s = re.sub("{ntyp}", str(ntyp), s)
        s = re.sub("{wmass}", str(wmass), s)
        s = re.sub("{atomic species}", atomic_spec, s)
        s = re.sub("{cell param}", cell, s)
        s = re.sub("{atomic positions}", atomic_pos, s)
        s = re.sub("{k points}", ' '.join([str(i) for i in self.kpoints]), s)

        inname = self.name + '_' + self.calc_type #+ str(uuid.uuid4())

        with open(inname + '.in', 'w') as f:
            f.write(s)

        return inname

    def read_cell(self):
        inname = self.name + '_' + self.calc_type
        with open(inname + '.out') as f:
            s = f.read()
        cell = re.findall("CELL_PARAMETERS.*\n(.*)\n(.*)\n(.*)\n{2}",s)[-1]
        cell = [[float(num) for num in line.split()] for line in cell]
        return cell

    def read_atomic_pos(self):
        inname = self.name + '_' + self.calc_type
        with open(inname + '.out') as f:
            s = f.read()
        search = "ATOMIC_POSITIONS.*"
        nat = sum([len(i) for i in self.atomic_pos])
        for i in range(nat):
            search += "\n(.*)"
        atomic_pos = {}

        for line in re.findall(search,s)[-1]:
            elem = line.split()
            if atomic_pos.has_key(elem[0]):
                atomic_pos[elem[0]].append([float(num) for num in elem[1:]])
            else:
                atomic_pos[elem[0]] =[[float(num) for num in elem[1:]]]
        return atomic_pos

    def read_energies(self):
        inname = self.name + '_' + self.calc_type
        with open(inname + '.out') as f:
            s = f.read()
        ene = re.findall("\! *total energy *\= *(.*) Ry",s)
        ene = [float(num) for num in ene]
        return ene

    def read_pressures(self):
        inname = self.name + '_' + self.calc_type
        with open(inname + '.out') as f:
            s = f.read()
        pres = re.findall("P *\= *(.*)",s)
        pres = [float(num) for num in pres]
        return ene


class phcalc:

    def __init__(self):
        self.type = 'ph.x'
        self.calc_type = "ph"
        self.name = "si"
        self.masses = {'Si':28.0855}
        self.qpoints = [8,8,8]


    def write_in(self):
        with open('/u/cb/cr/felixtherrien/scratch/pythonQE/canvas/ph_QEcanvas.in') as f:
            s = f.read()
        s = re.sub("{name}", self.name, s)
        s = re.sub("{nq1}", str(self.qpoints[0]), s)
        s = re.sub("{nq2}", str(self.qpoints[1]), s)
        s = re.sub("{nq3}", str(self.qpoints[2]), s)
        
        masses = ""
        for i, mass in enumerate(self.masses.itervalues()):
            masses += "amass(%d)=%f\n"%(i+1,mass)

        s = re.sub("{masses}", masses, s)

        inname = self.name + '_' + self.calc_type #+ str(uuid.uuid4())

        with open(inname + '.in', 'w') as f:
            f.write(s)

        return inname

class q2rcalc:

    def __init__(self):
        self.type = 'q2r.x'
        self.calc_type = "q2r"
        self.name = "si"

    def write_in(self):
        with open('/u/cb/cr/felixtherrien/scratch/pythonQE/canvas/q2r_QEcanvas.in') as f:
            s = f.read()
        s = re.sub("{name}", self.name, s)

        inname = self.name + '_' + self.calc_type #+ str(uuid.uuid4())

        with open(inname + '.in', 'w') as f:
            f.write(s)

        return inname

class matcalc:

    def __init__(self):
        self.type = 'matdyn.x'
        self.calc_type = "matdyn"
        self.name = "si"
        self.dos = False
        self.masses = {'Si':28.0855}
        self.kpoints = [6,6,6]
        self.path = [
            [0.0000000,   0.0000000,   0.0000000, 10],
            [0.7500000,   0.7500000,   0.0000000, 1 ],
            [0.2500000,   1.0000000,   0.2500000, 10],
            [0.0000000,   1.0000000,   0.0000000, 10],
            [0.0000000,   0.0000000,   0.0000000, 10],
            [0.5000000,   0.5000000,   0.5000000, 10],
            [0.7500000,   0.7500000,   0.0000000, 1 ],
            [0.2500000,   1.0000000,   0.2500000, 10],
            [0.5000000,   1.0000000,   0.0000000, 10],
            [0.0000000,   1.0000000,   0.0000000, 10],
            [0.5000000,   1.0000000,   0.0000000, 10],
            [0.5000000,   0.5000000,   0.5000000, 1 ]]

    def write_in(self):
        with open('/u/cb/cr/felixtherrien/scratch/pythonQE/canvas/matdyn_QEcanvas.in') as f:
            s = f.read()
        s = re.sub("{name}", self.name, s)
        s = re.sub("{dos}", '.true.' if self.dos else '.false.' , s)
        s = re.sub("{nk1}", str(self.kpoints[0]), s)
        s = re.sub("{nk2}", str(self.kpoints[1]), s)
        s = re.sub("{nk3}", str(self.kpoints[2]), s)
        s = re.sub("{len}", str(len(self.path)), s)

        path = ""
        for line in self.path:
            path += ' '.join([str(i) for i in line]) + '\n'

        s = re.sub("{path}", path, s)

        masses = ""
        for i, mass in enumerate(self.masses.itervalues()):
            masses += "amass(%d)=%f\n"%(i+1,mass)

        s = re.sub("{masses}", masses, s)

        inname = self.name + '_' + self.calc_type #+ str(uuid.uuid4())

        with open(inname + '.in', 'w') as f:
            f.write(s)

        return inname

    def read_eig(self):
        with open("PH/c.vecs") as f:
            s = f.read()
        q = re.findall(" *q = *([^ ]*) *([^ ]*) *([^ ]*)\n", s)
        all_freqs = re.findall(" *freq \(.*\) *= *(.*) \[THz\]", s)
        all_vecs = re.findall(" \( {1,2}(.*) {2,3}(.*) {4,5}(.*) {2,3}(.*) {4,5}(.*) {2,3}(.*)   \)", s)
        dim = len(all_freqs)/len(q)
        freqs = []
        for i,freq in enumerate(all_freqs):
            if i%dim == 0:
                freqs.append([float(freq)])
            else:
                freqs[-1].append(float(freq))
        
        vecs = []
        for i,vec in enumerate(all_vecs):
            vec_part = []
            for j,elem in  enumerate(vec):
                if j%2 == 0:
                    vec_part.append(float(elem))
                else:
                    vec_part[-1] = vec_part[-1] + float(elem)*1j
            
            if i%(dim*dim/3) == 0:
                vecs.append([vec_part])
            elif i%(dim/3) == 0:
                vecs[-1].append(vec_part)
            else:
                vecs[-1][-1].extend(vec_part)
        
        return qs, freqs, vecs

def submit_jobs(*calcs,**param):
    name = param.get('name','QErun.sh')
    np = param.get('np',16)
    with open(name,'w') as f:
        f.write('#! /bin/bash\n')
        for calc in calcs:
            inname = calc.write_in()
            if np > 1:
                f.write('srun ' + '-n ' + str(np) + 
                        ' ' + calc.type +
                        ' < ' + inname + '.in' +
                        ' > ' + inname + '.out\n' )
            else: 
                f.write(calc.type +
                        ' < ' + inname + '.in' +
                        ' > ' + inname + '.out\n' )
    call(['chmod','777', name])
    call('./' + name)

if __name__ == "__main__":
    pwrelax = pwcalc()
    submit_jobs(pwrelax)
    ene = pwrelax.read_energies()
    while (abs(ene[-1] - ene[-2]) > 1e-8):
        pwrelax.atomic_pos = pwrelax.read_atomic_pos()
        pwrelax.cell = pwrelax.read_cell()
        submit_jobs(pwrelax)
        ene = pwrelax.read_energies()
    pwscf = deepcopy(pwrelax)
    pwscf.calc_type = 'scf'
    pwscf.atomic_pos = pwrelax.read_atomic_pos()
    pwscf.cell = pwrelax.read_cell()

    ph = phcalc()
    q2r = q2rcalc()
    matdyn = matcalc()

    submit_jobs(pwscf, ph, q2r, matdyn)

    print matdyn.read_eig()
