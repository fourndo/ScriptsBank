# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 20:02:50 2018

@author: DominiqueFournier
"""
import re
import numpy as np
import matplotlib.pyplot as plt

work_dir = 'C:\\Users\\DominiqueFournier\\Desktop\\GoFundMe\\'

fileName = 'Boushi.dat'

with open(work_dir + fileName, encoding="utf8") as f:
    content = f.readlines()
    
amount = []
name = []
time = []
block = []
check1 = False
check2 = False
check3 = False
for line in content:
    
    if np.all(['supporter-amount' in line, '$' in line]):
       temp = re.findall('\d+',line)
       if len(temp) > 1:
           block += ["".join(map(str,temp))]
       else:
           block += temp
       check1 = True
       continue
       
           
    if check1:
        if "supporter-name" in line: 
            block += re.findall("\>(.*?)\<",line)
            check1 = False
            check2 = True
        else:
            check1 = False
            block = []
            continue
        
        continue
    
    if check2:        
        if "supporter-time" in line:
           
           if 'day' in line:
               block += [np.asarray(re.findall('\d+',line), dtype=int)*24]
           else:
               block += [np.asarray(re.findall('\d+',line), dtype=int)]
           
           check3 = True
           
        else:
            check2 = False
            block = []
            continue
        
    if check3:  
        amount += [block[0]]
        name += [block[1]]
        time += [block[2]]
        block = []
        check2 = False
        check3 = False
            
amount = np.asarray(amount).astype(float)
time = np.vstack(time).flatten()

# Re-order in time
tInd = np.argsort(time)[::-1]
clock = time[tInd]
cash = amount[tInd]
fund = np.cumsum(cash)

# Get unique time
gate, ind = np.unique(clock,return_index=True)
ind = ind[1:]-1

## Repeat for stanley


fileName = 'Stanley.dat'

with open(work_dir + fileName, encoding="utf8") as f:
    content = f.readlines()
    
amount = []
name = []
time = []
block = []
check1 = False
check2 = False
check3 = False
for line in content:
    
    if np.all(['supporter-amount' in line, '$' in line]):
       temp = re.findall('\d+',line)
       if len(temp) > 1:
           block += ["".join(map(str,temp))]
       else:
           block += temp
       check1 = True
       continue
       
           
    if check1:
        if "supporter-name" in line: 
            block += re.findall("\>(.*?)\<",line)
            check1 = False
            check2 = True
        else:
            check1 = False
            block = []
            continue
        
        continue
    
    if check2:        
        if "supporter-time" in line:
           
           if 'day' in line:
               block += [np.asarray(re.findall('\d+',line), dtype=int)*24]
           else:
               block += [np.asarray(re.findall('\d+',line), dtype=int)]
           
           check3 = True
           
        else:
            check2 = False
            block = []
            continue
        
    if check3:  
        amount += [block[0]]
        name += [block[1]]
        time += [block[2]]
        block = []
        check2 = False
        check3 = False
            
amount = np.asarray(amount).astype(float)
time = np.vstack(time).flatten()

# Re-order in time
tInd = np.argsort(time)[::-1]
clockStan = time[tInd]
cashStan = amount[tInd]
fundStan = np.cumsum(cashStan)

# Get unique time
gate, indStan = np.unique(clockStan,return_index=True)
indStan = indStan[1:]+1


plt.figure(figsize=(10,4))
axs = plt.subplot()
plt.grid(True)
# plot time scale
plt.bar(clock[ind], fund[ind], 10, color='r')
# plot time scale
plt.bar(clockStan[indStan], fundStan[indStan], 10, color='b', alpha=0.2)
axs.set_xlabel('Time (hours)')
axs.set_ylabel('Total raised (CAD)')
plt.gca().invert_xaxis()

#%%
plt.figure(figsize=(10,4))
axs = plt.subplot()
plt.hist(cashStan,200)
plt.hist(cash,200, color='r', alpha=0.75)
axs.set_xlim([0,500])
plt.grid(True)