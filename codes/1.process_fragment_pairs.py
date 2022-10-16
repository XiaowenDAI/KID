# -*- coding: utf-8 -*-
"""
Created on Tue Sep  6 09:26:57 2022

@author: Administrator
"""

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.Draw import MolsToGridImage
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.rdmolfiles import SmilesMolSupplier as sms

import os
import copy
import operator
import collections
import pandas as pd
from shutil import copyfile
from tkinter import _flatten
from rdkit.Chem.Descriptors import HeavyAtomCount

def Standardize(smi):
    mol = Chem.MolFromSmiles(smi)
    uncharger = rdMolStandardize.Uncharger()
    clean_mol = rdMolStandardize.Cleanup(mol)
    uncharged = uncharger.uncharge(clean_mol)
    return uncharged

def replaceHydrogen(smi):
    m = Chem.MolFromSmiles(smi,sanitize=True)
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

scaffold_decoration_pairs = []    
with open('paper_test/1.raw_scaffold_decoration_pairs.txt','r') as f:
    for line in f.readlines():
        l = line.replace("\n", "")
        l = l.split(',')
        scaffold_decoration_pairs.append(l)
    f.close()    

fragment_amino_acid_pairs = []    
with open('paper_test/1.raw_fragment_amino_acid_pairs.txt','r') as f:
    for line in f.readlines():
        l = line.replace("\n", "")
        l = l.split(',')
        fragment_amino_acid_pairs.append(l)
    f.close()    
    
ligands = []    
with open('1.ligands.smi','r') as f:
    for line in f.readlines():
        l = line.replace("\n", "")
        l = l.split(',')
        ligands.append(l)
    f.close() 

ligands_standerdized = []
for l in ligands:
    ligand_mol = Chem.MolFromSmiles(l[2])
    if ligand_mol is not None:
        ls = Standardize(Chem.MolToSmiles(ligand_mol))
        ligands_standerdized.append((l[0],l[1],Chem.MolToSmiles(ls,canonical=True)))

#==============================================================================
'''
    加完整分子信息以及去重
    fa_pairs去重依据：碎片及其来源分子、针对的氨基酸均一致
    sd_pairs去重依据：针对同一靶点相同的sd_pairs
'''
for p in scaffold_decoration_pairs:
    for item in ligands_standerdized:
        if p[1] == item[1]:
            p.append(item[2])

for p in fragment_amino_acid_pairs:
    for item in ligands_standerdized:
        if p[2] == item[1]:
            p.append(item[2])

with open('paper_test/2.sd_pairs_with_unbroken_ligands.txt','w') as f:
    for line in scaffold_decoration_pairs:
        for i in range(0,len(line)):
            if i<len(line)-1:
                f.write(line[i])
                f.write(',')
            else:
                f.write(line[i])
        f.write('\n')
    f.close()
    
with open('paper_test/2.fa_pairs_with_unbroken_ligands.txt','w') as f:
    for line in fragment_amino_acid_pairs:
        for i in range(0,len(line)):
            if i<len(line)-1:
                f.write(line[i])
                f.write(',')
            else:
                f.write(line[i])
        f.write('\n')
    f.close()

for i in range(0,len(scaffold_decoration_pairs)):
    for j in range(3,len(scaffold_decoration_pairs[i])):
        if '*' in scaffold_decoration_pairs[i][j]:
            scaffold_decoration_pairs[i][j] = Chem.MolToSmiles(Standardize(scaffold_decoration_pairs[i][j]))

#去重 ,修改了去重算法           
invalid_idx = []
for i in range(0,len(scaffold_decoration_pairs)):
    if i not in invalid_idx:
        for j in range(i+1,len(scaffold_decoration_pairs)):
            if operator.eq(scaffold_decoration_pairs[i][0],scaffold_decoration_pairs[j][0]):
                if operator.eq(scaffold_decoration_pairs[i][2:],scaffold_decoration_pairs[j][2:]):  #仅针对同一个化合物对同一个靶点去重
                    invalid_idx.append(j)
invalid_idx = list(set(invalid_idx))
invalid_idx = sorted(invalid_idx)       

for i in range(len(invalid_idx)-1, -1, -1):
    scaffold_decoration_pairs.pop(invalid_idx[i])

with open('paper_test/2.unique_sd_pairs_with_unbroken_ligands.txt','w') as f:
    for line in scaffold_decoration_pairs:
        for i in range(0,len(line)):
            if i<len(line)-1:
                f.write(line[i])
                f.write(',')
            else:
                f.write(line[i])
        f.write('\n')
    f.close()

invalid_idx = []
for i in range(0,len(fragment_amino_acid_pairs)):
    if i not in invalid_idx:
        for j in range(i+1,len(fragment_amino_acid_pairs)):
            if operator.eq(fragment_amino_acid_pairs[i][0],fragment_amino_acid_pairs[j][0]):
                if operator.eq(fragment_amino_acid_pairs[i][3:],fragment_amino_acid_pairs[j][3:]):  #仅针对同一个化合物对同一个靶点去重
                    invalid_idx.append(j)
invalid_idx = list(set(invalid_idx))
invalid_idx = sorted(invalid_idx)       

for i in range(len(invalid_idx)-1, -1, -1):
    fragment_amino_acid_pairs.pop(invalid_idx[i])

with open('paper_test/2.unique_fa_pairs_with_unbroken_ligands.txt','w') as f:
    for line in fragment_amino_acid_pairs:
        for i in range(0,len(line)):
            if i<len(line)-1:
                f.write(line[i])
                f.write(',')
            else:
                f.write(line[i])
        f.write('\n')
    f.close()
#==============================================================================
s_tmp = copy.deepcopy(scaffold_decoration_pairs)

pairsOnly = []
invalid_idx = []
for i in range(0,len(s_tmp)):
    for j in range(0,len(s_tmp[i])):
        if '*' in s_tmp[i][j] and '*' not in s_tmp[i][j-1]:  
            pairsOnly.append(s_tmp[i][j:])
            mol = Chem.MolFromSmiles(s_tmp[i][j])
            if HeavyAtomCount(mol)<=22:
                pass
            else:
                invalid_idx.append(i)
for i in range(0,len(pairsOnly)):
    if i not in invalid_idx:
        for j in range(i+1,len(pairsOnly)):
            if operator.eq(pairsOnly[i],pairsOnly[j]):
                invalid_idx.append(j)
invalid_idx = list(set(invalid_idx))
invalid_idx = sorted(invalid_idx)   
for i in range(len(invalid_idx)-1, -1, -1):
    scaffold_decoration_pairs.pop(invalid_idx[i]) 

with open('paper_test/2.sd_pairs_preprocessed_for_decorator.txt','w') as f:
    for line in scaffold_decoration_pairs:
        for i in range(0,len(line)):
            if i<len(line)-1:
                f.write(line[i])
                f.write(',')
            else:
                f.write(line[i])
        f.write('\n')
    f.close() 


#==================================试图查看和论文里原来写的数据为什么不一样
sd_pairs_3344 = []   
with open('2.scaffold_decoration_pairs_filtered.txt','r') as f:
    for line in f.readlines():
        l = line.replace("\n", "")
        l = l.split(',')
        sd_pairs_3344.append(l)
    f.close()    
    
sd_pairs_now = []
for i in range(0,len(scaffold_decoration_pairs)):
    for j in range(3,len(scaffold_decoration_pairs[i])):
        if '*' in scaffold_decoration_pairs[i][j] and '*' not in scaffold_decoration_pairs[i][j-1]:
            tmp = scaffold_decoration_pairs[i][j:]
            sd_pairs_now.append(tmp)

not_exist = []
for i in range(0,len(sd_pairs_now)):
    flag = False
    for j in range(0,len(sd_pairs_3344)):
        now = sd_pairs_now[i][:-1]
        if operator.eq(now,sd_pairs_3344[j]):
            flag = True
    if flag == False:
        not_exist.append(i)

            
            
            
            
            
            
            
            
            
            
            
