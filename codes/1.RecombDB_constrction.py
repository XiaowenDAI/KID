# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 20:47:35 2022

@author: fragment_replacement
"""

from rdkit import Chem
import os
import math
import pandas as pd
from tkinter import _flatten
from rdkit.Chem import rdFMCS
from rdkit.Chem.rdmolfiles import MolFromPDBFile
from rdkit.Chem.MolStandardize import rdMolStandardize

fragment_amino_acid_pairs = []
with open('../3.fragmentation_revised/4.fragment_amino_acid_pairs_filtered.txt','r') as f:
    for line in f.readlines():
        line = line.strip()
        l = line.split(',')
        fragment_amino_acid_pairs.append(l)
    f.close()

def Mol_Conn(s1,s2,attachment):
    conn = [] 
    m1 = Chem.MolFromSmiles(s1)
    m2 = Chem.MolFromSmiles(s2)    
    m = Chem.CombineMols(m1, m2)
    mw = Chem.RWMol(m)
    for i in range(0,m1.GetNumAtoms()):
        if m1.GetAtomWithIdx(i).GetSymbol() == attachment: 
            neighbor1_idx = m1.GetAtomWithIdx(i).GetNeighbors()[0].GetIdx() 
            for j in range(0,m2.GetNumAtoms()):
                if m2.GetAtomWithIdx(j).GetSymbol() == attachment:
                    neighbor2_idx = m2.GetAtomWithIdx(j).GetNeighbors()[0].GetIdx()
                    mw.AddBond(neighbor1_idx,neighbor2_idx+m1.GetNumAtoms(),Chem.BondType.SINGLE)
                    mw.RemoveAtom(j+m1.GetNumAtoms())
                    mw.RemoveAtom(i)
                    conn.append(Chem.MolToSmiles(mw)) 
                    mw = Chem.RWMol(m) 
    return conn

def Conn(data1,data2,attachment): 
    connected = []
    for i in range(len(data1)):
        if i == 0:
            for s2 in data2:
                connected.extend(Mol_Conn(data1[i],s2,attachment))
        else:
            len_tmp = len(connected)
            for j in range(len(connected)):
                if Chem.MolFromSmiles(connected[j]) is not None:
                    connected.extend(Mol_Conn(connected[j],data1[i],attachment))
            del connected[:len_tmp]
    return connected

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
    
def Standardize(smi):
    mol = Chem.MolFromSmiles(smi)
    uncharger = rdMolStandardize.Uncharger()
    clean_mol = rdMolStandardize.Cleanup(mol)
    uncharged = uncharger.uncharge(clean_mol)
    res = Chem.MolToSmiles(uncharged)
    return res

def fragmention(mol,hit):
    mw = Chem.RWMol(mol)
    link = []
    for h in hit:
        for a in mol.GetAtomWithIdx(h).GetNeighbors():
            if a.GetIdx() not in hit:
                link.append([h,a.GetIdx()])                                       
    for l in link:
        atom1 = l[0]
        atom2 = l[1]
        mw.RemoveBond(atom1,atom2)
        mw.AddAtom(Chem.Atom(0))
        mw.AddBond(atom1,mw.GetNumAtoms()-1,Chem.BondType.SINGLE)
        mw.AddAtom(Chem.Atom(0))
        mw.AddBond(atom2,mw.GetNumAtoms()-1,Chem.BondType.SINGLE) 
    frags = Chem.MolToSmiles(mw).split('.')
    frags_without_att = [Standardize(replaceHydrogen(s)) for s in frags]
    h = Chem.MolFragmentToSmiles(mol, hit) 
    h = Chem.MolToSmiles(Chem.MolFromSmiles(h))
    h_idx = frags_without_att.index(h)
    n_atta = frags[h_idx].count('*')
    frags.pop(h_idx)
    frags = [f for f in frags if f !='[Xe]']
    return frags,n_atta

def replace(mol,to_replace,amino_acids):
    choice = []
    for record in fragment_amino_acid_pairs:
        record_aa = record[4:]
        record_without_idx = [r[0:3] for r in record_aa]
        if set(amino_acids).issubset(set(record_without_idx)):
            choice.append(record[0])
    choice = list(set(choice))
    
    hit = list(mol.GetSubstructMatch(to_replace))
    print(hit)
    frags_for_connection,n_atta = fragmention(mol,hit)

    invalid = []    
    for i in range(len(choice)):
        if choice[i].count('*') < n_atta:
            invalid.append(i)
    for i in range(len(invalid)-1, -1, -1):
        choice.pop(invalid[i])
            
    connected = Conn(frags_for_connection,choice,'*')        
    res = list(set([Standardize(replaceHydrogen(s)) for s in connected]))
    return res

#---------------------------Pre analysis-------------------------
err = []
for p in fragment_amino_acid_pairs:
    frag = Chem.MolFromSmiles(Standardize(replaceHydrogen(p[0])))
    m = Standardize(Chem.MolFromSmiles(p[3]))
    if frag is not None and m.GetSubstructMatch(frag):
        pass
    else:
        err.append(p)

#----------------------database construction--------------------   
base = []
failed = []
num = len(fragment_amino_acid_pairs)
for p in fragment_amino_acid_pairs[:5]:
    i = 1
    print('----------'+p[2]+'-'+str(i)+'/'+str(num)+'-processing'+'----------')
    try:
        mol = Chem.MolFromSmiles(p[3])
        to_replace = Chem.MolFromSmiles(Standardize(replaceHydrogen(p[0])))
        aa = p[4:]
        base.append(replace(mol,to_replace,aa))
    except:
        failed.append(p)
    i+=1
molecule_base = list(_flatten(base))















