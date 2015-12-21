'''
Created on Jul 17, 2013

@author: dominiquef
'''
def read_input(input_file):
    """Read input files for the lp_norm script"""
    import os
    import re
    from numpy import asarray
    
    os.system('pwd')
    fid = open(input_file,'r')
    
    line = fid.readline()
    l_input  = line.split('!')
    Home_dir = l_input[0].rstrip()
    
    line = fid.readline()
    l_input  = line.split('!')
    UBC_dir = l_input[0].rstrip()
    
    line = fid.readline()
    l_input = line.split('!') 
    iter_start = int(l_input[0])
    
    line = fid.readline()
    l_input = line.split('!') 
    chifact = float(l_input[0])
    
    line = fid.readline()
    l_input = line.split('!')
    pp = l_input[0].split(' ')
    p = asarray(pp)
    
    line = fid.readline()
    l_input = line.split('!')
    qq = l_input[0].split(' ')
    q = asarray(qq)
    
    line = fid.readline()
    l_input = line.split('!')
    ll = l_input[0].split(' ')
    l = asarray(ll)
    
    line = fid.readline()
    l_input = line.split('!') 
    cool_beta = float(l_input[0])
    
    line = fid.readline()
    l_input = line.split('!') 
    iter_max = float(l_input[0])
      
    return Home_dir, UBC_dir, iter_start, chifact, p, q, l, cool_beta, iter_max