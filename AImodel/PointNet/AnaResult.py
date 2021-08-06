#%%
import numpy as np
def GetCrystal(Cate_File0,Cate_File1,Cate_File2,Cate_File3):
    Cateresult0,Cateresult1,Cateresult2,Cateresult3 = np.load(Cate_File0),np.load(Cate_File1),np.load(Cate_File2),np.load(Cate_File3)
    CateSta = []
    for i in range(0,4):
        CateSta.append((np.sum(Cateresult0==i))/len(Cateresult0))
    for i in range(0,4):
        CateSta.append((np.sum(Cateresult1==i))/len(Cateresult1))
    for i in range(0,4):
        CateSta.append((np.sum(Cateresult2==i))/len(Cateresult2))
    for i in range(0,4):
        CateSta.append((np.sum(Cateresult3==i))/len(Cateresult3))
    CateSta=np.array(CateSta).reshape(4,-1)
    return CateSta
Crystal_Pointnet=GetCrystal(r'./Crystal/Pred_Crystal0.npy',r'./Crystal/Pred_Crystal1.npy',r'./Crystal/Pred_Crystal2.npy',r'./Crystal/Pred_Crystal3.npy')  
#%%
import matplotlib.pyplot as plt
f, (ax, ax2) = plt.subplots(2,1,sharex=True)
Crylabel = ('LiquidCluster','BCCCluster','FCCCluster','HCPCluster')
Colormap = ('c','red','g','b') 
for i in range(4):
    ax.bar(np.arange(4)+i*0.2,Crystal_Pointnet[:,i],label=Crylabel[i],width=0.2,color=Colormap[i])
    ax2.bar(np.arange(4)+i*0.2,Crystal_Pointnet[:,i],label=Crylabel[i],width=0.2,color=Colormap[i])
ax.set_ylim(.92,1.01)  # outliers only
ax.set_yticks(np.arange(3)*0.04+0.92)
ax2.set_ylim(0,.08)
ax2.set_yticks(np.arange(3)*0.04)
plt.xticks(np.arange(4)+0.3,('Liquid','BCC','FCC','HCP'))
ax.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax.xaxis.tick_top()
ax.tick_params(labeltop=False)
ax2.xaxis.tick_bottom()
d = .02  # how big to make the diagonal lines in axes coordinates
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal
kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
ax.legend(loc='upper left', fontsize='x-small')
plt.show(f)
f.savefig('HistogramofCry.jpg')
plt.close()
#%%
Pred_PointNet=np.load(r'./Predictset/Pred_PointNet.npy')
Pred_PointNetS = Pred_PointNet[len(Pred_PointNet)-18000:]
Pred_W4W6=np.load(r'./W4W6/Pred_W4W6.npy')
Pred_W4W6S = Pred_W4W6[len(Pred_W4W6)-18000:]
Pred_NewW4W6=np.load(r'./New_W4W6/Pred_W4W6.npy')
Pred_NewW4W6S = Pred_NewW4W6[len(Pred_NewW4W6)-18000:]
Pred_MLP=np.load(r'./MLP/Pred_MLP.npy')
Pred_MLPS = Pred_MLP[len(Pred_MLP)-18000:]
Pred_Xcate=np.load(r'./Xcate/Pred_Xcate.npy')
Pred_AveXcate=np.load(r'./AveXcate/Pred_AveXcate.npy')
#%%
from sklearn.metrics import confusion_matrix
from sklearn.metrics import ConfusionMatrixDisplay
Normalornot='true' #None 'true'
cnf_matrix = confusion_matrix(Pred_Xcate,Pred_AveXcate,normalize=Normalornot)
cm_display = ConfusionMatrixDisplay(cnf_matrix,display_labels=('Liquid','BCC','FCC','HCP')).plot()
plt.show()
