from mbiomas import fca
# mbiomas is our own toolbox published here: https://gitlab.uliege.be/mipi/published-software/mbiomas-core/-/tree/92e072ff8032a7048bebea9ffe1738c548b5327a/mbiomas
from datetime import datetime
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


print("Read Online Flow cytometer data!")

FCF, FCD, FCN = fca.read.multitel(fields=['FSC-A','SSC-A','FL1-A','FSC-H','SSC-H','FL1-H','Width','Time'],clean=True,delete=False)
print("Doing statistics")
FCStatistics = fca.statistics.calculate(FCF,FCD,scale='log10',fields=['FSC-A','SSC-A','FL1-A'])
FCStatistics=fca.statistics.entropy(FCF,FCD,FCStatistics, fields=['FL1-A'])
print('Gating and Selecting GFP- and GFP+ populations')
FCF,FCD,INFO=fca.gating.line(FCF,FCD,fields=['FL1-A'],thresholds=[4],name='GFP')
fca.write.FCFilescsv(FCF,name='FCF')


#####Plotting
batch=9.97
Glc_Ara=46.7
Ara_Glc=96.4
BioEnd=122
#####################################



fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [7, 0.8]}, sharex=True)

TS = fca.plot.Plot(FCF)
TS.time_scatter(ax1, FCF, FCD, y_label='FL1-A', y_log_lim=[1, 7], y_res=50)
ax1.set_ylabel('Fluorescence [-]', fontsize=20)

ax1.plot(FCStatistics.data["FL1-A"]["Time"], 10**FCStatistics.data["FL1-A"]["Median"], color="coral") #, label="Median FL1-A"

ax1.axvline(batch, color='r', label='Initiate Chemostat')
ax1.axvline(Ara_Glc, color='black') 
ax1.axvline(Glc_Ara, color='black') 

ax1.tick_params(labelbottom=False)

ax2.set_xlim(0, BioEnd) 

ax2.axvspan(0, Glc_Ara, color='#B452CD', alpha=0.5)  
ax2.axvspan(Glc_Ara, Ara_Glc, color='#FFD700', alpha=0.5)  
ax2.axvspan(Ara_Glc, BioEnd, color='#B452CD', alpha=0.5)  

ax2.set_yticks([])

ax2.set_xlabel('Time [h]', fontsize=20)

ax1.tick_params(axis='both', which='major', labelsize=20)
ax1.tick_params(axis='both', which='minor', labelsize=15)
ax2.tick_params(axis='both', which='major', labelsize=20)
ax2.tick_params(axis='both', which='minor', labelsize=15)
ax2.set_xticks([0, 47,96,122]) 

fig.set_size_inches(11.36, 5.62)
fig.subplots_adjust(hspace=0.01, right=0.56) 
plt.show()
fig.savefig('FL_TimeScatter-2.png', dpi=100)
