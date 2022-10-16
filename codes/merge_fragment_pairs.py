# -*- coding: utf-8 -*-
"""
@author: vv
"""

import os
import math
import pandas as pd
import numpy as np
import networkx as nx
from tkinter import _flatten
from collections import Counter
from rdkit import Chem
from rdkit.Chem.rdchem import Mol,Atom
from rdkit.Chem import BRICS
from rdkit.Chem.rdmolfiles import MolFromPDBFile
from rdkit.Chem.MolStandardize import rdMolStandardize
from queue import Queue
import threading

overview = pd.read_csv('../2.filtered_pdb/2.overview_human_orthosteric.csv')
for i in range(0,len(overview)):
    if 'E+' in overview.loc[i,'pdb']:
        overview.loc[i,'pdb'] = overview.loc[i,'pdb'].replace('.00E+', 'e')

def mol_with_atom_index(mol):
    atoms = mol.GetNumAtoms()
    for idx in range(atoms):
        mol.GetAtomWithIdx(idx).SetProp('molAtomMapNumber',str(mol.GetAtomWithIdx(idx).GetIdx()))
    return mol

def mol_without_atom_index(mol):
    atoms = mol.GetNumAtoms()
    for idx in range(atoms):
        mol.GetAtomWithIdx(idx).ClearProp('molAtomMapNumber')
    return mol

def Standardize(smi):
    mol = Chem.MolFromSmiles(smi)
    uncharger = rdMolStandardize.Uncharger()
    clean_mol = rdMolStandardize.Cleanup(mol)
    uncharged = uncharger.uncharge(clean_mol)
    res = Chem.MolToSmiles(uncharged)
    return res

def read_pdb_by_name(file):
    pdb = []
    path = '../2.filtered_pdb/2.filtered_pdb_orthosteric'
    filepath = path + '/' + file
    with open(filepath,'r') as f:
        next(f)
        for line in f.readlines():
            tmp = line.split()
            if len(tmp)>5 and len(tmp[4])>1:
                chain = tmp[4]
                tmp[4] = chain[0]
                tmp.insert(5,chain[1:])    
            if len(tmp[0])>6: 
                residue_class = tmp[0]
                tmp[0] = residue_class[0:6]
                tmp.insert(1,residue_class[6:])                          
            pdb.append(tmp)
            f.close()
    return pdb

def read_pdb_by_target(target):
    pdb_list = []
    path = '../2.filtered_pdb/2.filtered_pdb_orthosteric'
    files = os.listdir(path)
    with open('NoneType.txt','w') as out_file:
        for file in files:
            pdb = []
            if file.lower().startswith(target.lower()):
                filepath = path + '/' + file
                target = file.split('_')[0]
                pdb_name = file.split('_')[1]
                pdb.append(target)
                pdb.append(pdb_name)
                mol = MolFromPDBFile(filepath)
                if mol is not None:
                    pdb.append(mol)
                    with open(filepath,'r') as f:
                        next(f)
                        for line in f.readlines():
                            tmp = line.split()
                            if len(tmp)>5 and len(tmp[4])>1:
                                chain = tmp[4]
                                tmp[4] = chain[0]
                                tmp.insert(5,chain[1:])    
                            if len(tmp[0])>6: 
                                residue_class = tmp[0]
                                tmp[0] = residue_class[0:6]
                                tmp.insert(1,residue_class[6:])                          
                            pdb.append(tmp)
                        f.close()
                    pdb[2] = mol_with_atom_index(pdb[2])
                    pdb_list.append(pdb)
                else:
                    out_file.write("%s\n"%file) 
        out_file.close()
    return pdb_list
    
def read_pdb():
    '''
    read all pdb files
    '''
    pdb_list = []
    path = '../2.filtered_pdb/2.filtered_pdb_orthosteric'
    files = os.listdir(path)
    with open('NoneType.txt','w') as out_file:
        for file in files:
            pdb = []
            filepath = path + '/' + file
            target = file.split('_')[0]
            pdb_name = file.split('_')[1]
            pdb.append(target)
            pdb.append(pdb_name)
            m = MolFromPDBFile(filepath)
            if m is not None:
                pdb.append(m)
                with open(filepath,'r') as f:
                    next(f)
                    for line in f.readlines():
                        tmp = line.split()
                        if len(tmp)>5 and len(tmp[4])>1:
                            chain = tmp[4]
                            tmp[4] = chain[0]
                            tmp.insert(5,chain[1:])  
                        if len(tmp[0])>6: 
                            residue_class = tmp[0]
                            tmp[0] = residue_class[0:6]
                            tmp.insert(1,residue_class[6:])               
                        pdb.append(tmp)
                    f.close()
                pdb[2] = mol_with_atom_index(pdb[2])
                pdb_list.append(pdb)
            else:
                out_file.write("%s\n"%file) 
        out_file.close()
    return pdb_list    

def seperate_receptor_and_ligand(pdb):
    target_name = pdb[0]
    pdb_name = pdb[1]
    mol = pdb[2]
    receptor = []
    ligand = []
    connect = []
    for i in range(3,len(pdb)):
        if pdb[i][0] == 'ATOM':
            receptor.append(pdb[i])
        if pdb[i][0] == 'HETATM':
            ligand.append(pdb[i])
        if pdb[i][0] == 'CONECT':
            connect.append(pdb[i])
    return target_name,pdb_name,receptor,ligand,connect,mol

def eucliDist(A,B):
    return math.sqrt(sum([(a - b)**2 for (a,b) in zip(A,B)]))

def cal_centroid(atoms):
    coord = np.array([list(map(float, a[6:9])) for a in atoms])
    centroid = np.mean(coord, axis=0)
    return centroid

def NorO(atom):
    if atom == 'N' or atom == 'O':
        return True
    else:
        return False

def equal(s1,s2):
    for s in s2:
        if s not in s1:
            return False
    return True

def exist_in_route(route,atom):
    for line in route:
        if atom.GetIdx() in line:
            return True
    return False

def in_ring_with_NO(m,atom,rings):
    for i in range(0,len(rings)): 
        if atom.GetIdx() in rings[i]:
            for x in rings[i]:
                if m.GetAtomWithIdx(x).GetSymbol()=='N' or m.GetAtomWithIdx(x).GetSymbol()=='O':
                    return True
    return False

def delete_attachment(smi):
    m = Chem.MolFromSmiles(smi)
    if m is not None:
        mw = Chem.RWMol(m)
        atta_idx = []
        for i in range(0,mw.GetNumAtoms()):
            if mw.GetAtomWithIdx(i).GetSymbol() == '*':
                atta_idx.append(i)
        for i in reversed(atta_idx):
            mw.RemoveAtom(i)
        return Chem.MolToSmiles(mw)
    else:
        return '[Xe]'
    
def replaceHydrogen(smi):
    m = Chem.MolFromSmiles(smi)
    if m is not None:
        mw = Chem.RWMol(m)
        atta_idx = []
        for i in range(0,mw.GetNumAtoms()):
            if mw.GetAtomWithIdx(i).GetSymbol() == '*':
                atta_idx.append(i)
        for i in reversed(atta_idx):
            mw.ReplaceAtom(i,Chem.rdchem.Atom(1))
        return Chem.MolToSmiles(mw)
    else:
        return '[Xe]'

def is_exist_in_double_list(target, element):
    for line in target:
        if element in line:
            return True
    return False

def fragmentation(m,hit):
    '''
    input: m -> Chem.Mol
           hit -> location for fragmentation
    '''
    bricsBonds = list(BRICS.FindBRICSBonds(m))
    reactionable_sites = [list(p[0]) for p in bricsBonds]
    for r in reactionable_sites:
        r.sort()
    mw = Chem.RWMol(m)
    link = []
    for h in hit:
        for a in m.GetAtomWithIdx(h).GetNeighbors():
            if a.GetIdx() not in hit:
                l = [h,a.GetIdx()]
                l.sort()
                if l in reactionable_sites:
                    link.append(l)  
    if len(link) > 0:                                   
        for l in link:
            atom1 = l[0]
            atom2 = l[1]
            mw.RemoveBond(atom1,atom2)
            mw.AddAtom(Chem.Atom(0))
            mw.AddBond(atom1,mw.GetNumAtoms()-1,Chem.BondType.SINGLE)
            mw.AddAtom(Chem.Atom(0))
            mw.AddBond(atom2,mw.GetNumAtoms()-1,Chem.BondType.SINGLE) 
        frags = Chem.MolToSmiles(mw).split('.')
        frags_filtered = []
        for f in frags:
            if '*' in f:
                frags_filtered.append(f)
        scaffold = Chem.MolFragmentToSmiles(m, hit)
        if len(frags_filtered) > 1:
            try:
                frags_without_att = [delete_attachment(s) for s in frags_filtered]
                s_idx = frags_without_att.index(scaffold)
            except:
                n_atta = [f.count('*') for f in frags_filtered]
                s_idx = n_atta.index(max(n_atta))
            mid = frags_filtered[0]
            frags_filtered[0] = frags_filtered[s_idx]
            frags_filtered[s_idx] = mid
            return frags_filtered
        else:
            return ['Xe']
    else:
        return ['Xe']
    
def hit_for_ring_atom(m,idx):
    ringinfo = m.GetRingInfo()
    atomrings = list(list(items) for items in list(ringinfo.AtomRings()))
    G = nx.Graph()   
    G.add_nodes_from(sum(atomrings, []))
    q = [[(s[i],s[i+1]) for i in range(len(s)-1)] for s in atomrings]
    for i in q:
        G.add_edges_from(i)
    rings = [list(i) for i in nx.connected_components(G)] 
    
    atom = m.GetAtomWithIdx(idx)    
    hit = []
    for r in rings:
        if atom.GetIdx() in r:
            hit_tmp = r
    while(hit_tmp):
        t = hit_tmp[0]
        hit.append(t)
        a = m.GetAtomWithIdx(t)
        neighbors = list(a.GetNeighbors())
        for n in neighbors:
            if not n.IsInRing() and n.GetIdx() not in hit:
                hit_tmp.append(n.GetIdx())
        hit_tmp.remove(t)
    return hit
            
def break_molecule_by_index(m,mapnum):
    '''
    m:Chem.rdchem.Mol
    mapnum: int
    '''
    for atom in m.GetAtoms():
        if atom.GetAtomMapNum() == mapnum:
            idx = atom.GetIdx()
    state = True 
    hit = []
    atom = m.GetAtomWithIdx(idx)
    if atom.IsInRing():
        hit.append(hit_for_ring_atom(m,idx))
    else:
        tmp = [[atom.GetIdx()]] 
        flag = False  
        i = 0
        while(flag==False):
            ti = tmp[i]
            neighbors = []
            for t in ti:
                neighbors_tmp = list(m.GetAtomWithIdx(t).GetNeighbors())
                for n in neighbors_tmp:                   
                    if not is_exist_in_double_list(tmp,n.GetIdx()):
                        neighbors.append(n.GetIdx())
            if len(neighbors) == 0:
                state = False
                break
            tmp.append(neighbors)
            for n in tmp[i+1]:
                a = m.GetAtomWithIdx(n)
                if a.IsInRing() and flag==False:
                    hit.append(hit_for_ring_atom(m,n)) 
            if len(hit) >= 1:
                flag = True
            i = i+1
    frags = []
    for h in hit:
        f = fragmentation(m,h)
        if f != ['Xe']:
            frags.append(f)
    if len(frags) == 0:
        state = False
    return frags,state

def combine_dup(pairs):
    mark = -1
    for i in range(0,len(pairs)):
        if pairs[i][0] != -1:
            pairs[i].insert(0,1)
            for j in range(i+1,len(pairs)):
                if pairs[i][1] == pairs[j][0]:
                    pairs[i][0] = pairs[i][0] + 1
                    pairs[i].extend(pairs[j][3:])
                    pairs[j].insert(0,mark)
    
    # one fragment -> more than 1 amino acid 
    dup = []
    for i in range(0,len(pairs)):
        if pairs[i][0] == mark:
            dup.append(i)
    for i in range(len(dup)-1, -1, -1):
        pairs.pop(dup[i])
    
    # more than 1 fragment -> same amino acid
    redundant = []
    for i in range(0,len(pairs)):
        for j in range(i+1,len(pairs)):
            if pairs[i][4:] == pairs[j][4:]:
                if pairs[i][0]>pairs[j][0]:
                    redundant.append(j)
                elif pairs[i][0]<pairs[j][0]:
                    redundant.append(i)
    for i in range(len(redundant)-1, -1, -1):
        pairs.pop(redundant[i])
        
    for i in range(0,len(pairs)):
        pairs[i] = sorted(set(pairs[i]),key=pairs[i].index)  
    
    res=[p[1:] for p in pairs]
    return res

def choose_hinge(pairs):
    mark = -1
    for i in range(0,len(pairs)):
        if pairs[i][0] != -1:
            pairs[i].insert(0,1)
            for j in range(i+1,len(pairs)):
                if pairs[i][1] == pairs[j][0]:
                    pairs[i][0] = pairs[i][0] + 1
                    pairs[j][0] = mark
    count = [p[0] for p in pairs]
    max_count = max(count)
    idx = []
    for i in range(0,len(pairs)):
        if pairs[i][0] == max_count:
            idx.append(i)
    res = [pairs[i] for i in idx]
    for r in res:
        del(r[0])
    return res 
    
def split_ligand_in_hinge_binding_loc(pdbs):
    with open('error_hinge.txt', 'w') as out_file:
        centroid = [-3.134,21.976,47.086]
        scaffold_decoration_pairs = []
        num = len(pdbs)
        pc = 1
        for pdb in pdbs:
            ac = [] 
            pairs = [] 
            state = False
            target_name,pdb_name,receptor,ligand,connect,mol = seperate_receptor_and_ligand(pdb)
            print('----------'+pdb_name+'-'+str(pc)+'/'+str(num)+'-processing'+'----------')
            ligand_ortho = [l for l in ligand if l[3] == overview.loc[overview['pdb']==pdb_name]['ligand'].reset_index(drop=True).get(0)]
            dev = int(ligand[-1][1])-mol.GetNumAtoms()
            ligand_ortho_idx = [int(l[1])-dev-1 for l in ligand_ortho]
            ligand_mol = Chem.MolFromSmiles(Chem.MolFragmentToSmiles(mol, ligand_ortho_idx))
            for l in ligand_ortho :
                if NorO(l[11]):
                    coord_l = list(map(float, l[6:9]))
                    for r in receptor:
                        coord_r = list(map(float, r[6:9]))
                        if NorO(r[11]) and eucliDist(centroid,coord_r)<=6 and eucliDist(coord_r,coord_l)<=3.5 and eucliDist(coord_r,coord_l)>=2.5:
                            ac.append(r[3]+r[5])
                            atom_num = int(l[1])-1
                            try:
                                frag_pairs,s = break_molecule_by_index(ligand_mol,atom_num)
                                if s == True:
                                    for p in frag_pairs:
                                        p_tmp = []
                                        frag_mol = [Chem.MolFromSmiles(f) for f in p]
                                        for i in range(0,len(frag_mol)):
                                            if frag_mol[i] is not None:
                                                p_tmp.append(Chem.MolToSmiles(mol_without_atom_index(frag_mol[i])))
                                        pairs.append(p_tmp)
                                    state = True
                            except:
                                out_file.write('FragNonetype'+pdb[1]+'\n')                                       
            if len(pairs)>=1 and s==True:
                pairs = choose_hinge(pairs)
                for p in pairs:
                    scaffold_decoration_pairs.append([pdb[0],pdb[1]]+list(set(ac))+p)
                state = True
            if state == False:
                out_file.write(pdb[1]+'\n') 
            pc += 1
    out_file.close()
    return scaffold_decoration_pairs

def merge_fragment_amino_acid_pairs_hd(pdbs):
    with open('error_hd.txt', 'w') as out_file:
        fragment_amino_acid_pairs = []
        num = len(pdbs)
        i = 1
        for pdb in pdbs: 
            state = False
            pairs = []
            target_name,pdb_name,receptor,ligand,connect,mol = seperate_receptor_and_ligand(pdb)
            ligand_ortho = [l for l in ligand if l[3] == overview.loc[overview['pdb']==pdb_name]['ligand'].reset_index(drop=True).get(0)]
            dev = int(ligand[-1][1])-mol.GetNumAtoms()
            ligand_ortho_idx = [int(l[1])-dev-1 for l in ligand_ortho]
            ligand_mol = Chem.MolFromSmiles(Chem.MolFragmentToSmiles(mol, ligand_ortho_idx))
            print('----------'+pdb_name+'-'+str(i)+'/'+str(num)+'-processing'+'----------')
            for l in ligand_ortho:
                if NorO(l[11]):
                    coord_l = list(map(float, l[6:9]))
                    for r in receptor:
                        coord_r = list(map(float, r[6:9]))
                        if NorO(r[11]) and eucliDist(coord_r,coord_l)<=3.5:                       
                            atom_num = int(l[1])-1
                            try:
                                frag_pairs,s = break_molecule_by_index(ligand_mol,atom_num)
                                if s == True:
                                    frag = [f[0] for f in frag_pairs]
                                    frag_mol = [Chem.MolFromSmiles(f) for f in frag]
                                    for f in frag_mol:
                                        if f is not None:
                                            m = mol_without_atom_index(f)
                                            if m is not None:
                                                frag_without_atom_index = Chem.MolToSmiles(m)
                                                info = [frag_without_atom_index]
                                                info.extend([pdb[0],pdb[1],r[3]+r[5]])
                                                pairs.append(info)
                                    state = True
                            except:
                                continue
            if state == False:
                out_file.write(pdb[1]+'\n') 
            elif len(pairs)>1:
                deduped = combine_dup(pairs)
                fragment_amino_acid_pairs.extend(deduped)
            else:
                fragment_amino_acid_pairs.extend(pairs)
            i += 1
    return fragment_amino_acid_pairs

def merge_ligands(pdbs):
    with open('merge_ligand_error.txt', 'w') as out_file:
        ligands = []
        num = len(pdbs)
        i = 1
        for pdb in pdbs: 
            target_name,pdb_name,receptor,ligand,connect,mol = seperate_receptor_and_ligand(pdb)
            print('----------'+pdb_name+'-'+str(i)+'/'+str(num)+'-processing'+'----------')
            hit = [int(l[1])-1 for l in ligand if l[3] == overview.loc[overview['pdb']==pdb_name]['ligand'].reset_index(drop=True).get(0)]
            print(hit)
            try:
                l = Chem.MolFragmentToSmiles(mol, hit)
                ligand_mol = mol_without_atom_index(Chem.MolFromSmiles(l))
                ligand_smi = Chem.MolToSmiles(ligand_mol)
                ligands.append(list((target_name,pdb_name,ligand_smi)))
            except:
                out_file.write(pdb_name+'\n')
            i+=1
    out_file.close()
    return ligands 

#--------------------fragmentation---------------------------------------------
pdbs = read_pdb()
scaffold_decoration_pairs = split_ligand_in_hinge_binding_loc(pdbs)
fragment_amino_acid_pairs = merge_fragment_amino_acid_pairs_hd(pdbs)
ligands = merge_ligands(pdbs) 

#-----------fragmentation for pdbs with restorage error------------------------
rehandle = []    
with open('NoneType.txt','r') as f:
    for line in f.readlines():
        line = line.replace('.pdb\n','')
        line = line.split('_')
        rehandle.append(line)  
    f.close()

root = '../1.raw data/HUMAN'
rehandle_route = []
for i in range(0,len(rehandle)):
    if len(rehandle[i]) == 2:
        rehandle_route.append(root + '/' + rehandle[i][0] + '/' + rehandle[i][1] + '/complex.pdb')
    if len(rehandle[i]) == 3:
        rehandle_route.append(root + '/' + rehandle[i][0] + '/' + rehandle[i][1] +  '_' + rehandle[i][2] + '/complex.pdb')
    else:    
        rehandle_route.append(root + '/' + rehandle[i][0] + '/' + rehandle[i][1] +  '_' + rehandle[i][2] +  '_' + rehandle[i][3] + '/complex.pdb')

header = ['HEADER','REMARK','HELIX','SHEET','CRYST1']
rehandle_files = []
for r in rehandle_route:
    pdb = []
    target = r.split('/')[3]
    pdb_name = r.split('/')[4].split('_')[0]
    pdb.append(target)
    pdb.append(pdb_name)
    m = MolFromPDBFile(r,sanitize=False,removeHs=False)
    if m is not None:
        pdb.append(m)
        with open(r,'r') as f:
            for line in f.readlines():
                tmp = line.split()
                if tmp[0] in header:
                    continue
                elif len(tmp)>5 and len(tmp[4])>1:
                    chain = tmp[4]
                    tmp[4] = chain[0]
                    tmp.insert(5,chain[1:])     
                elif len(tmp[0])>6: 
                    residue_class = tmp[0]
                    tmp[0] = residue_class[0:6]
                    tmp.insert(1,residue_class[6:])                          
                pdb.append(tmp)
            f.close()
            pdb[2] = mol_with_atom_index(pdb[2])
        rehandle_files.append(pdb)
            
fragment_amino_acid_pairs = merge_fragment_amino_acid_pairs_hd(rehandle_files) 
scaffold_decoration_pairs_rehandle = split_ligand_in_hinge_binding_loc(rehandle_files)    

     
