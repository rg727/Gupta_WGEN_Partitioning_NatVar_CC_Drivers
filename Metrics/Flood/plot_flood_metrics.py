# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 11:51:02 2023

@author: rg727
"""

import pandas as pd 
import seaborn as sns
import matplotlib.pyplot as plt

#Plot flood metrics 

#10 Year event w/ model reference 
paleo=pd.read_csv("Tuolumne_GEV_10yr_30yrMW.csv")


for i in range (28900):
    if paleo.iloc[i,1]>100:
        paleo.iloc[i,1]=float("NaN")
        

area_mi2_TLG = 1538


paleo.iloc[:,1]=paleo.iloc[:,1]*area_mi2_TLG*5280*5280*(1/4356)*4046.86*(1/1000)

#Historical 10-year flow derived from sac-sma baseline 
historical=	46.81393*area_mi2_TLG*5280*5280*(1/4356)*4046.86*(1/1000)


fig = plt.figure()
ax = plt.subplot(111)
ax=sns.lineplot(data=paleo, x="years", y="flow",color="#606C38")
ax.axhline(y=historical, color='black', linestyle='--')
ax.set_xlim(1400,2017)
ax.tick_params(axis='both', which='major', labelsize=14)
ax.set_xlabel("Year",fontsize=14)
ax.set_ylabel("3-Day Flow Volume ($m^3$)",fontsize=14)
fig.savefig('Tuolumne_10_yr_flood.pdf')


#100 Year Event w/ model reference


#100 Year event w/ model reference 
paleo=pd.read_csv("Tuolumne_GEV_100yr_30yrMW.csv")

for i in range (28900):
    if paleo.iloc[i,1]>330:
        paleo.iloc[i,1]=float("NaN")
        

paleo.iloc[:,1]=paleo.iloc[:,1]*area_mi2_TLG*5280*5280*(1/4356)*4046.86*(1/1000)

#Historical flow 
historical=	95.66451*area_mi2_TLG*5280*5280*(1/4356)*4046.86*(1/1000)


fig = plt.figure()
ax = plt.subplot(111)
ax=sns.lineplot(data=paleo, x="years", y="flow",color="#606C38")
ax.axhline(y=historical, color='black', linestyle='--')
ax.set_xlim(1400,2017)
ax.set_xlabel("Year")
ax.set_ylabel("3-Day Flow Volume ($m^3$)")
fig.savefig('Tuolumne_100_yr_flood.pdf')  


