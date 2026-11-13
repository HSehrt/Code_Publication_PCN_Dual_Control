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
start=10.78
batch=16.28
Temp42=63.62
Temp30=135.78
BioEnd=160 #181.78
#####################################

#  3 Subplots
fig, (ax1, ax2, ax3) = plt.subplots(
    3, 1,
    gridspec_kw={'height_ratios': [7, 0.8, 0.8]},
    sharex=True
)
# --- 1. Subplot ---
TS = fca.plot.Plot(FCF)
TS.time_scatter(ax1, FCF, FCD, y_label='FL1-A', y_log_lim=[1, 7], y_res=50)
ax1.set_ylabel('Fluorescence', fontsize=20)
# Median plot
ax1.plot(
    FCStatistics.data["FL1-A"]["Time"],
    10**FCStatistics.data["FL1-A"]["Median"],
    color="coral"
)

ax1.axvline(batch, color='r')
ax1.axvline(Temp42, color='black')
ax1.axvline(Temp30, color='black')
ax1.tick_params(labelbottom=False)
# --- 2. Subplot ---
ax2.set_xlim(start, BioEnd)
ax1.set_xlim(start, BioEnd)
ax1.set_yticks([100, 10000, 1000000])

ax2.axvspan(0, batch, color="#5C5C5C", alpha=0.5)
ax2.axvspan(batch, Temp42, color="#C6C5C3", alpha=0.5)  
ax2.axvspan(Temp42, Temp30, color="#292929CF", alpha=0.5)  
ax2.axvspan(Temp30, BioEnd, color='#C6C5C3', alpha=0.5)
ax2.set_yticks([])
ax2.tick_params(labelbottom=False) 

# --- 3. Subplot ---
ax3.set_xlim(start, BioEnd)

ax3.axvspan(start, batch, color='#B452CD', alpha=0.5)
ax3.axvspan(batch, Temp42, color='#FFD700', alpha=0.5)
ax3.axvspan(Temp42, Temp30, color='#B452CD', alpha=0.5)
ax3.axvspan(Temp30, BioEnd, color='#FFD700', alpha=0.5)

ax3.set_yticks([])
ax3.set_xlabel('Time [h]', fontsize=20)

ax1.tick_params(axis='both', which='major', labelsize=20)
ax1.tick_params(axis='both', which='minor', labelsize=15)

ax2.tick_params(axis='both', which='major', labelsize=20)
ax2.tick_params(axis='both', which='minor', labelsize=15)

ax3.tick_params(axis='both', which='major', labelsize=20)
ax3.tick_params(axis='both', which='minor', labelsize=15)

ax3.set_xticks([16, 64, 136, 160])
# --- Layout ---
fig.set_size_inches(8, 4)
fig.subplots_adjust(hspace=0.01, right=0.56)
plt.show()
fig.savefig('FL_TimeScatter-3.png', dpi=100)