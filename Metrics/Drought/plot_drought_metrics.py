# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 12:14:01 2023

@author: rg727
"""

import pandas as pd 
import seaborn as sns
import matplotlib.pyplot as plt


#Drought Occurrence

paleo=pd.read_csv("drought_occurrence.csv")


fig = plt.figure()
ax = plt.subplot(111)
ax=sns.lineplot(data=paleo, x="years", y="percentage",color="#BC6C25")
ax.set_xlim(1400,1986)
#ax.set_ylim(3,8)
ax.tick_params(axis='both', which='major', labelsize=14)
ax.set_ylabel('Percentage of Window',fontsize=14)
ax.set_xlabel("Year",fontsize=14)

fig.savefig('drought_occurrence.pdf')  


#Drought severity 

paleo=pd.read_csv("drought_severity.csv")


fig = plt.figure()
ax = plt.subplot(111)
ax=sns.lineplot(data=paleo, x="years", y="min_ssi",color="#BC6C25")
ax.set_xlim(1400,1986)
ax.set_ylabel('Min SSI',fontsize=14)
ax.set_xlabel("Year",fontsize=14)
ax.invert_yaxis()
ax.tick_params(axis='both', which='major', labelsize=14)

fig.savefig('drought_severity.pdf')  

#Drought Duration

paleo=pd.read_csv("drought_duration_one.csv")

        
fig = plt.figure()
ax = plt.subplot(111)
ax=sns.lineplot(data=paleo, x="years", y="months",color="#BC6C25")
ax.set_xlim(1400,1986)
ax.tick_params(axis='both', which='major', labelsize=14)
ax.set_ylabel('Months',fontsize=14)
ax.set_xlabel("Year",fontsize=14)

fig.savefig('drought_duration.pdf')  
