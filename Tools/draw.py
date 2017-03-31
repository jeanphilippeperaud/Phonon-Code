import numpy as np
import matplotlib
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt

fig, ax = plt.subplots()
patches = []
xmin = 0
xmax = 0
ymin = 0
ymax = 0
with open("T_detectors.txt") as T:
    Nd = int(T.readline())
    det = [[] for i in xrange(Nd)]
    for i in xrange(Nd):
	Line = T.readline()
	words = Line.split(' ')
        for word in words:
            det[i].append(float(word))
        polygon = Polygon(np.array([[det[i][0],det[i][1]],[det[i][2],det[i][3]],[det[i][4],det[i][5]],[det[i][6],det[i][7]]]))
        xmin = min(min([det[i][2*j] for j in xrange(4)]),xmin) 
        xmax = max(max([det[i][2*j] for j in xrange(4)]),xmax)
        ymin = min(min([det[i][2*j+1] for j in xrange(4)]),xmin) 
        ymax = max(max([det[i][2*j+1] for j in xrange(4)]),xmax)

        patches.append(polygon)
        
p = PatchCollection(patches)
colors = []
with open("results_T.txt") as Tr:
    for i in xrange(Nd):
	Line = Tr.readline()
	words = Line.split(' ')
        colors.append(float(words[3]))

print colors
p.set_array(np.array(colors))
ax.add_collection(p)
fig.colorbar(p, ax=ax)
axes = plt.gca()
axes.set_xlim([xmin,xmax])
axes.set_ylim([ymin,ymax])
plt.show()            

