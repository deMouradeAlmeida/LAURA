import LAURA_aux as aux
"""Cite De Moura et al. 2018,
Spectroscopic and asteroseismic analysis of the secondary clump red giant HD226808"""
"""Bruno Lustosa de Moura and Leandro de Almeida 2017-2018"""

print('#-###################################################')
print('Initiating the LAURA package')
print('#-###################################################')

#All packages
aux.TS_PS(ID,Qtin,Qtfin,cleaning,oversample,path = 'TEMP/')
aux.FIT_PS(ID)
aux.LS(ID,kernerforca,Pbeg,Pend,PeriodMod,path='TEMP/')
aux.OBS(ID)
aux.MOD(ID, largeSep,smallSep,size,Radial_Orders=5)
aux.ECHELLE(ID, deltav)
#Everything is going to be placed in the TEMP\ folder
print ('Done!')
