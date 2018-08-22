import re
import numpy as np
import glob
import os
import uuid
from subprocess import call
from copy import deepcopy
from pylada.crystal import Structure

MODPATH = os.path.dirname(os.path.realpath(__file__))

def int2str(n):
    return '%d'%n

def fl2str(n):
    return '% 10.6F'%n

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
        self.unit = "alat"
        self.atomic_pos = {'Si':[[0,0,0],[0.25,0.25,0.25]]}
        self.kpoints = [8,8,8]


    def write_in(self):
        with open(MODPATH + '/canvas/pw_QEcanvas.in') as f:
            s = f.read()
        s = re.sub("{calculation type}", self.calc_type, s)
        s = re.sub("{restart mode}", self.restart_mode, s)
        s = re.sub("{name}", self.name, s)
        s = re.sub("{pseudo dir}", self.pseudo_dir, s)
        s = re.sub("{unit}", self.unit, s)
        s = re.sub("{celldm}", fl2str(self.celldm), s)
        s = re.sub("{ecutwfc}", fl2str(self.ecutwfc), s)
        s = re.sub("{ecutrho}", fl2str(self.ecutrho), s)
        s = re.sub("{nbnd}", int2str(self.nbnd), s)
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
                atomic_spec += (atom_type + ' ' + fl2str(self.masses[atom_type]) + ' ' +
                                pseudoname +
                                '\n')
            else:
                raise AttributeError("Masses has no key %s"%atom_type)
            for atom in self.atomic_pos[atom_type]:
                wmass += self.masses[atom_type]
                nat += 1
                atomic_pos += (atom_type + ' ' +
                               ' '.join([fl2str(i) for i in atom]) + '\n')
                
        cell = ""
        for line in self.cell.T:
            cell += ' '.join([fl2str(i) for i in line]) + '\n'
        
        s = re.sub("{nat}", int2str(nat), s)
        s = re.sub("{ntyp}", int2str(ntyp), s)
        s = re.sub("{wmass}", fl2str(wmass), s)
        s = re.sub("{atomic species}", atomic_spec, s)
        s = re.sub("{cell param}", cell, s)
        s = re.sub("{atomic positions}", atomic_pos, s)
        s = re.sub("{k points}", ' '.join([int2str(i) for i in self.kpoints]), s)

        inname = self.name + '_' + self.calc_type #+ str(uuid.uuid4())

        with open(inname + '.in', 'w') as f:
            f.write(s)

        return inname

    def from_pylada(self, struct):
        self.cell = struct.cell
        
        self.atomic_pos = {}
        if self.unit == "alat":
            for atom in struct:
                if self.atomic_pos.has_key(atom.type):
                    self.atomic_pos[atom.type].append(list(atom.pos))
                else:
                    self.atomic_pos[atom.type] = [list(atom.pos)]
        elif self.unit == "crystal":
            icell = np.linalg.inv(self.cell)
            for atom in struct:
                if self.atomic_pos.has_key(atom.type):
                    self.atomic_pos[atom.type].append(list(icell.dot(atom.pos)))
                else:
                    self.atomic_pos[atom.type] = [list(icell.dot(atom.pos))]
        else:
            raise RuntimeError("Wrong unit must be alat or crystal")


    def to_pylada(self):
        A = Structure(self.cell)
        if self.unit == "alat":
            for elem in self.atomic_pos:
                for pos in self.atomic_pos[elem]:
                    A.add_atom(pos[0],pos[1],pos[2],elem)
        elif self.unit == "crystal":
            for elem in self.atomic_pos:
                for pos in self.atomic_pos[elem]:
                    self.cell.dos(pos)
                    A.add_atom(pos[0],pos[1],pos[2],elem)
        else:
            raise RuntimeError("Wrong unit must be alat or crystal")
        return A

    def read_cell(self):
        inname = self.name + '_' + self.calc_type
        with open(inname + '.out') as f:
            s = f.read()
        cell = re.findall("CELL_PARAMETERS.*\n(.*)\n(.*)\n(.*)\n{2}",s)[-1]
        cell = [[float(num) for num in line.split()] for line in cell]
        return cell.T

    def read_atomic_pos(self):
        inname = self.name + '_' + self.calc_type
        with open(inname + '.out') as f:
            s = f.readlines()
        
        nat = sum([len(i) for i in self.atomic_pos.itervalues()])

        for i, line in enumerate(s):
            if re.search("ATOMIC_POSITIONS", line):
               lineNum = i     
        
        atomic_pos = {}
        for line in s[lineNum + 1:lineNum + nat + 1]:
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
        self.ldisp = True
        self.masses = {'Si':28.0855}
        self.qpoints = [8,8,8]
        self.qlist = None


    def write_in(self):
        with open(MODPATH + '/canvas/ph_QEcanvas.in') as f:
            s = f.read()
        s = re.sub("{name}", self.name, s)
        s = re.sub("{ldisp}", '.true.' if self.ldisp else '.false.' , s)
        s = re.sub("{nq1}", int2str(self.qpoints[0]), s)
        s = re.sub("{nq2}", int2str(self.qpoints[1]), s)
        s = re.sub("{nq3}", int2str(self.qpoints[2]), s)
        if not self.ldisp:
            qlist = ""
            for line in self.qlist:
                qlist += ' '.join([fl2str(i) for i in line]) + '\n'

            s = re.sub("{qlist}", qlist, s)
        else:
            s = re.sub("{qlist}", "", s)
        
        masses = ""
        for i, mass in enumerate(self.masses.itervalues()):
            masses += "amass(%d)=%s\n"%(i+1,fl2str(mass))

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
        with open(MODPATH + '/canvas/q2r_QEcanvas.in') as f:
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
        self.path = np.array([
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
            [0.5000000,   0.5000000,   0.5000000, 1 ]]).T

    def write_in(self):
        with open(MODPATH + '/canvas/matdyn_QEcanvas.in') as f:
            s = f.read()
        s = re.sub("{name}", self.name, s)
        s = re.sub("{dos}", '.true.' if self.dos else '.false.' , s)
        s = re.sub("{nk1}", int2str(self.kpoints[0]), s)
        s = re.sub("{nk2}", int2str(self.kpoints[1]), s)
        s = re.sub("{nk3}", int2str(self.kpoints[2]), s)
        s = re.sub("{len}", int2str(len(self.path)), s)

        path = ""
        for line in self.path.T:
            path += ' '.join([fl2str(i) for i in line[:3]]) \
                     + ' ' + int2str(line[-1])+ '\n'

        s = re.sub("{path}", path, s)

        masses = ""
        for i, mass in enumerate(self.masses.itervalues()):
            masses += "amass(%d)=%s\n"%(i+1,fl2str(mass))

        s = re.sub("{masses}", masses, s)

        inname = self.name + '_' + self.calc_type #+ str(uuid.uuid4())

        with open(inname + '.in', 'w') as f:
            f.write(s)

        return inname

    def read_eig(self):
        with open(self.name + ".vecs") as f:
            s = f.read()
        qs = re.findall(" *q = *([^ ]*) *([^ ]*) *([^ ]*)\n", s)
        qs = [[float(i) for i in q] for q in qs]
        all_freqs = re.findall(" *freq \(.*\) *= *(.*) \[THz\]", s)
        all_vecs = re.findall(" \( {1,2}(.*) {2,3}(.*) {4,5}(.*) {2,3}(.*) {4,5}(.*) {2,3}(.*)   \)", s)
        dim = len(all_freqs)/len(qs)
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

class dyncalc:

    def __init__(self):
        self.type = 'dynmat.x'
        self.calc_type = "dynmat"
        self.name = "si"
        
    def write_in(self):
        with open(MODPATH + '/canvas/dynmat_QEcanvas.in') as f:
            s = f.read()
        s = re.sub("{name}", self.name, s)
        
        inname = self.name + '_' + self.calc_type #+ str(uuid.uuid4())

        with open(inname + '.in', 'w') as f:
            f.write(s)

        return inname

    def read_eig(self):
        with open(self.name + ".vecs") as f:
            s = f.read()
        qs = re.findall(" *q = *([^ ]*) *([^ ]*) *([^ ]*)\n", s)
        qs = [[float(i) for i in q] for q in qs]
        all_freqs = re.findall(" *freq \(.*\) *= *(.*) \[THz\]", s)
        all_vecs = re.findall(" \( {1,2}(.*) {2,3}(.*) {4,5}(.*) {2,3}(.*) {4,5}(.*) {2,3}(.*)   \)", s)
        dim = len(all_freqs)/len(qs)
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
    name = param.get('name','QErun_' + str(uuid.uuid4()) + '.sh')
    np = param.get('np',16)
    
    with open(name,'w') as f:
        f.write('#! /bin/bash\n')
        for calc in calcs:
            inname = calc.write_in()
            if True:
                if (calc.type == "pw.x"):
                    f.write('srun ' + '-n ' + int2str(np) +
                            ' ' + calc.type +
                            ' < ' + inname + '.in' +
                            ' > ' + inname + '.out\n' )
                else:
                    f.write('srun ' + '-n ' + int2str(np) + 
                            ' ' + calc.type +
                            ' < ' + inname + '.in' +
                            ' > ' + inname + '.out\n' )
            else: # TODO: It seems like not calling srun does not work
                f.write(calc.type +
                        ' < ' + inname + '.in' +
                        ' > ' + inname + '.out\n' )
    call(['chmod','777', name])
    call('./' + name)
    call(['rm', name])

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
