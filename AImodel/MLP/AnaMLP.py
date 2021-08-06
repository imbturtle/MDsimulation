import numpy as np
import glob
from sklearn.model_selection import train_test_split
from keras.utils import np_utils
from tensorflow.keras import models
from tensorflow.keras import layers
from tensorflow.keras.callbacks import TensorBoard
#%% Data preparation
def Get_dataset(Pre_file):    
    Pre_list=glob.glob(Pre_file)
    Pre_list.sort()
    Data = np.loadtxt(Pre_list[0],delimiter=' ')
    for i in range(1,len(Pre_list)):Data = np.concatenate((Data,np.loadtxt(Pre_list[i], delimiter=' ')))
    return Data
Train_X,Test_X,Train_y,Test_ya = train_test_split(Get_dataset(r'./Trainset/BOOP*')[:,(1,8,9)],Get_dataset(r'./Trainset/Crystal*'),test_size=0.0005,random_state=0,stratify=Get_dataset(r'./Trainset/Crystal*'))
Train_y= np_utils.to_categorical(Train_y, num_classes=4)
np.save("./Trainset/Train_X.npy",Train_X)
np.save("./Trainset/Train_y.npy",Train_y)
np.save("./Trainset/Test_X.npy",Test_X)
np.save("./Trainset/Test_ya.npy",Test_ya)
#%% Draw Train
Test_X = np.load('./Trainset/Test_X.npy')
Test_ya = np.load('./Trainset/Test_ya.npy')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
LiquidGrid = np.delete(Test_X,np.where(Test_ya != 0),axis=0)
BCCGrid = np.delete(Test_X,np.where(Test_ya != 1),axis=0)
FCCGrid = np.delete(Test_X,np.where(Test_ya != 2),axis=0)
HCPGrid = np.delete(Test_X,np.where(Test_ya != 3),axis=0)
fig = plt.figure(figsize=[10.,10.])
gridplt = fig.add_subplot(111, projection='3d')
gridplt.scatter(LiquidGrid[:,1],LiquidGrid[:,2],LiquidGrid[:,0],c='m',s=10,alpha=0.8,label='Liquid')
gridplt.scatter(BCCGrid[:,1],BCCGrid[:,2],BCCGrid[:,0],c='b',s=10,alpha=0.8,label='BCC')
gridplt.scatter(FCCGrid[:,1],FCCGrid[:,2],FCCGrid[:,0],c='g',s=10,alpha=0.8,label='FCC')
gridplt.scatter(HCPGrid[:,1],HCPGrid[:,2],HCPGrid[:,0],c='r',s=10,alpha=0.8,label='HCP')
gridplt.set_xlabel('W4')
gridplt.set_ylabel('W6')
gridplt.set_zlabel('Connect')
gridplt.legend(loc='upper right', fontsize='large')
gridplt.view_init(elev=5,azim=340)
#plt.show(gridplt)
plt.savefig('Traindata_W4W6scatter.jpg')
plt.close()
fig = plt.figure(figsize=[10.,10.])
gridplt = fig.add_subplot(111, projection='3d')
#gridplt.scatter(LiquidGrid[:,1],LiquidGrid[:,2],LiquidGrid[:,0],c='c',s=10,alpha=0.8,label='Liquid')
gridplt.scatter(BCCGrid[:,1],BCCGrid[:,2],BCCGrid[:,0],c='b',s=10,alpha=0.8,label='BCC')
gridplt.scatter(FCCGrid[:,1],FCCGrid[:,2],FCCGrid[:,0],c='g',s=10,alpha=0.8,label='FCC')
gridplt.scatter(HCPGrid[:,1],HCPGrid[:,2],HCPGrid[:,0],c='r',s=10,alpha=0.8,label='HCP')
gridplt.set_xlabel('W4')
gridplt.set_ylabel('W6')
gridplt.set_zlabel('Connect')
gridplt.legend(loc='upper right', fontsize='large')
gridplt.view_init(elev=90,azim=360)
#plt.show(gridplt)
plt.savefig('Traindata_W4W6scatter_1.jpg')
plt.close()
#%% Build model
def MLP(Train_X,Train_y,Train_epochs=50):
    model = models.Sequential() 
    model.add(layers.Dense(units=64,activation='relu',name='Den1',input_shape=(Train_X.shape[1],)))
    model.add(layers.Dense(units=64,activation='relu',name='Den2'))
    model.add(layers.Dense(units=64,activation='relu',name='Den3'))
    model.add(layers.Dropout(rate=0.5,name='Drop'))
    model.add(layers.Dense(units=4,activation='softmax',name='Den4'))
    model.compile(optimizer = 'adam',loss='categorical_crossentropy',metrics=['accuracy'])
    model.summary()
    tbCallBack = TensorBoard(log_dir=r'./Tensorboard',update_freq="epoch",write_graph=True,write_images=False,)
    model.fit(Train_X,Train_y,validation_split=0.2,epochs=Train_epochs,batch_size=128,callbacks=[tbCallBack])
    model.save(r'./classic_mlpmodel.h5')
    return model
model=MLP(Train_X,Train_y,5)
Test_y=np.argmax(model.predict(Test_X),axis=1)
np.save("./Trainset/Test_y.npy",Test_y)
#%% MD prediction
Pred_X=Get_dataset(r'./Predictset/BOOP*')[:,(1,8,9)]
Pred_MLP=np.argmax(model.predict(Pred_X),axis=1) 

Pred_Model(Pred_X)
Pred_W4W6=Get_dataset(r'./W4W6/TanaR*')
Pred_NewW4W6=Get_dataset(r'./New_W4W6/TanaR*')
Pred_Xcate=Get_dataset(r'./Xcate/XR*')
Pred_AveXcate=Get_dataset(r'./AveXcate/XR*')
np.save('./Predictset/Pred_X.npy',Pred_X)
np.save('./Predictset/Pred_MLP.npy',Pred_MLP)
np.save('./W4W6/Pred_W4W6.npy',Pred_W4W6)
np.save(r'./New_W4W6/Pred_W4W6.npy',Pred_NewW4W6)
np.save(r'./Xcate/Pred_Xcate.npy',Pred_Xcate)
np.save(r'./AveXcate/Pred_AveXcate.npy',Pred_AveXcate)


def Get_dataset(Train_file=r'./Predictset/BOOP*',Ans_file=r'./W4W6/Tana*'):  
    Train_list ,Ans_list= glob.glob(Train_file),glob.glob(Ans_file)
    Train_list.sort(),Ans_list.sort()
    Data,Answer = np.loadtxt(Train_list[0], delimiter=' '),np.loadtxt(Ans_list[0],delimiter=' ')
    for i in range(1,len(Train_list)):
        Data = np.concatenate((Data,np.loadtxt(Train_list[i], delimiter=' ')))
        Answer = np.concatenate((Answer,np.loadtxt(Ans_list[i],delimiter=' ')))
    Data=Data[:,(1,8,9)]
    return Data,Answer
Pred_X,Pred_W4W6a = Get_dataset()
Pred_MLP = np.argmax(model.predict(Pred_X),axis=1)
np.save("./Predictset/Pred_MLP.npy",Pred_MLP)
Pred_X,Pred_NewW4W6a=Get_dataset(Ans_file=r'./New_W4W6/Tana*')
#%%
#%%
def Get_dataset(Pre_file):    
    Pre_list=glob.glob(Pre_file)
    Pre_list.sort()
    Data = np.loadtxt(Pre_list[0],delimiter=' ')
    for i in range(1,len(Pre_list)):Data = np.concatenate((Data,np.loadtxt(Pre_list[i], delimiter=' ')))
    return Data
def Pred_Model(Pred_X,Pred_file):
    Pred_Y = np.argmax(model.predict(Pred_X),axis=1)
    np.save(Pred_file,Pred_Y)
    return 
Pred_X=Get_dataset(r'./Predictset/NB*.txt').reshape((-1,NearestP,3,1))
Pred_W4W6=Get_dataset(r'./W4W6/TanaR*')

np.save("./Predictset/Pred_X.npy",Pred_X)
np.save("./W4W6/Pred_W4W6.npy",Pred_W4W6)
Pred_PointNet=np.argmax(model.predict(Pred_X),axis=1)
np.save("./Predictset/Pred_PointNet.npy",Pred_PointNet)
#%% Draw connect<9
Pred_Pointnety = np.load('./Pointnet/Predict_y.npy')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
w4space = np.linspace(0.2,-0.2,20)
w6space = np.linspace(0.02,-0.02,20)
conspace = np.linspace(0,8,9)
W4space, W6space = np.meshgrid(w4space,w6space)
W4space, Conspace = np.meshgrid(W4space,conspace)
W6space, Conspace = np.meshgrid(W6space,conspace)
W4space, W6space, Conspace=W4space.reshape(-1,1),W6space.reshape(-1,1),Conspace.reshape(-1,1)
grid_X = np.concatenate((Conspace,W4space,W6space),axis=1)
grid_Y = np.array(model.predict_classes(grid_X))
LiquidGrid = np.delete(grid_X,np.where(grid_Y != 0),axis=0)
BCCGrid = np.delete(grid_X,np.where(grid_Y != 1),axis=0)
FCCGrid = np.delete(grid_X,np.where(grid_Y != 2),axis=0)
HCPGrid = np.delete(grid_X,np.where(grid_Y != 3),axis=0)
fig = plt.figure(figsize=[10.,10.])
gridplt = fig.add_subplot(111, projection='3d')
gridplt.scatter(LiquidGrid[:,1],LiquidGrid[:,2],LiquidGrid[:,0],c='c',s=10,alpha=0.8,label='Liquid')
gridplt.scatter(BCCGrid[:,1],BCCGrid[:,2],BCCGrid[:,0],c='r',s=10,alpha=0.8,label='BCC')
gridplt.scatter(FCCGrid[:,1],FCCGrid[:,2],FCCGrid[:,0],c='g',s=10,alpha=0.8,label='FCC')
gridplt.scatter(HCPGrid[:,1],HCPGrid[:,2],HCPGrid[:,0],c='b',s=10,alpha=0.8,label='HCP')
gridplt.set_xlabel('W4')
gridplt.set_ylabel('W6')
gridplt.set_zlabel('Connect')
gridplt.legend(loc='upper right', fontsize='large')
gridplt.view_init(elev=5,azim=230)
#plt.show(gridplt)
plt.savefig('MLP_GRID_0to8.jpg')
plt.close()
fig = plt.figure(figsize=[10.,10.])
gridplt = fig.add_subplot(111, projection='3d')
gridplt.scatter(LiquidGrid[:,1],LiquidGrid[:,2],LiquidGrid[:,0],c='c',s=10,alpha=0.8,label='Liquid')
gridplt.scatter(BCCGrid[:,1],BCCGrid[:,2],BCCGrid[:,0],c='r',s=10,alpha=0.8,label='BCC')
gridplt.scatter(FCCGrid[:,1],FCCGrid[:,2],FCCGrid[:,0],c='g',s=10,alpha=0.8,label='FCC')
gridplt.scatter(HCPGrid[:,1],HCPGrid[:,2],HCPGrid[:,0],c='b',s=10,alpha=0.8,label='HCP')
gridplt.set_xlabel('W4')
gridplt.set_ylabel('W6')
gridplt.set_zlabel('Connect')
gridplt.legend(loc='upper right', fontsize='large')
gridplt.view_init(elev=90,azim=270)
#plt.show(gridplt)
plt.savefig('MLP_GRID_0to8_1.jpg')
plt.close()
#%% Draw connect>9
w4space = np.linspace(0.2,-0.2,20)
w6space = np.linspace(0.02,-0.02,20)
conspace = np.linspace(9,14,6)
W4space, W6space = np.meshgrid(w4space,w6space)
W4space, Conspace = np.meshgrid(W4space,conspace)
W6space, Conspace = np.meshgrid(W6space,conspace)
W4space, W6space, Conspace=W4space.reshape(-1,1),W6space.reshape(-1,1),Conspace.reshape(-1,1)
grid_X = np.concatenate((Conspace,W4space,W6space),axis=1)
grid_Y = np.array(model.predict_classes(grid_X))
LiquidGrid = np.delete(grid_X,np.where(grid_Y != 0),axis=0)
BCCGrid = np.delete(grid_X,np.where(grid_Y != 1),axis=0)
FCCGrid = np.delete(grid_X,np.where(grid_Y != 2),axis=0)
HCPGrid = np.delete(grid_X,np.where(grid_Y != 3),axis=0)
fig = plt.figure(figsize=[10.,10.])
gridplt = fig.add_subplot(111, projection='3d')
gridplt.scatter(BCCGrid[:,1],BCCGrid[:,2],BCCGrid[:,0],c='r',s=10,alpha=0.8,label='BCC')
gridplt.scatter(FCCGrid[:,1],FCCGrid[:,2],FCCGrid[:,0],c='g',s=10,alpha=0.8,label='FCC')
gridplt.scatter(HCPGrid[:,1],HCPGrid[:,2],HCPGrid[:,0],c='b',s=10,alpha=0.8,label='HCP')
gridplt.scatter(LiquidGrid[:,1],LiquidGrid[:,2],LiquidGrid[:,0],c='c',s=10,alpha=0.8,label='Liquid')
gridplt.set_xlabel('W4')
gridplt.set_ylabel('W6')
gridplt.set_zlabel('Connect')
gridplt.legend(loc='upper right', fontsize='large')
gridplt.view_init(elev=5,azim=230)
#plt.show(gridplt)
plt.savefig('MLP_GRID_9to14.jpg')
plt.close()
fig = plt.figure(figsize=[10.,10.])
gridplt = fig.add_subplot(111, projection='3d')
gridplt.scatter(BCCGrid[:,1],BCCGrid[:,2],BCCGrid[:,0],c='r',s=10,alpha=0.8,label='BCC')
gridplt.scatter(FCCGrid[:,1],FCCGrid[:,2],FCCGrid[:,0],c='g',s=10,alpha=0.8,label='FCC')
gridplt.scatter(HCPGrid[:,1],HCPGrid[:,2],HCPGrid[:,0],c='b',s=10,alpha=0.8,label='HCP')
gridplt.scatter(LiquidGrid[:,1],LiquidGrid[:,2],LiquidGrid[:,0],c='c',s=10,alpha=0.8,label='Liquid')
gridplt.set_xlabel('W4')
gridplt.set_ylabel('W6')
gridplt.set_zlabel('Connect')
gridplt.legend(loc='upper right', fontsize='large')
gridplt.view_init(elev=90,azim=270)
#plt.show(gridplt)
plt.savefig('MLP_GRID_9to14_1.jpg')
plt.close()
#%% Confusion matrix AIvstanaka
from sklearn.metrics import confusion_matrix
from sklearn.metrics import ConfusionMatrixDisplay
Normalornot=None #None 'true'
cnf_matrix = confusion_matrix(Test_ya,Test_y,
                              normalize=Normalornot)
cm_display = ConfusionMatrixDisplay(cnf_matrix,display_labels=
                                    ('Liquid','BCC','FCC','HCP')).plot()
#plt.show()
plt.savefig('Conf_MLP_test.jpg')
plt.close()
cnf_matrix = confusion_matrix(Pred_Mlpa,Pred_Pointnety,
                              normalize=Normalornot)
cm_display = ConfusionMatrixDisplay(cnf_matrix,display_labels=
                                    ('Liquid','BCC','FCC','HCP')).plot()
#plt.show()
plt.savefig('Conf_MLP_to_Pointnet.jpg')
plt.close()
cnf_matrix = confusion_matrix(Pred_W4W6a,Pred_Mlpa,
                              normalize=Normalornot)
cm_display = ConfusionMatrixDisplay(cnf_matrix,display_labels=
                                    ('Liquid','BCC','FCC','HCP')).plot()
#plt.show()
plt.savefig('Conf_W4W6_to_MLP.jpg')
plt.close()
cnf_matrix = confusion_matrix(Pred_NewW4W6a,Pred_Mlpa,
                              normalize=Normalornot)
cm_display = ConfusionMatrixDisplay(cnf_matrix,display_labels=
                                    ('Liquid','BCC','FCC','HCP')).plot()
#plt.show()
plt.savefig('Conf_NewW4W6_to_MLP.jpg')
plt.close()
#%%
Crystal_X,Crystal_W4W6a=Get_dataset(r'./Crystal/BOOP*',r'./Crystal/Tana*')
Crystal_Mlpa=model.predict_classes(Crystal_X)
#%%
cnf_matrix = confusion_matrix(Crystal_W4W6a,Crystal_Mlpa,
                              normalize=Normalornot)
cm_display = ConfusionMatrixDisplay(cnf_matrix,display_labels=
                                    ('Liquid','BCC','FCC','HCP')).plot()
#plt.show()
plt.savefig('Conf_Crystal_W4W6_to_MLP.jpg')
plt.close()
#%%
Pred_X=Pred_X.reshape((10,-1,3))
Pred_y=Pred_y.reshape((10,-1))
Pred_ya=Pred_ya.reshape((10,-1))
Sta_Pred_y=[]
for i in range(0,4):
    for j in range(10):Sta_Pred_y.append((np.sum(Pred_y[j]==i)))
Sta_Pred_y=np.array(Sta_Pred_y).reshape((4,10))
Sta_Pred_ya=[]
for i in range(0,4):
    for j in range(10):Sta_Pred_ya.append((np.sum(Pred_ya[j]==i)))
Sta_Pred_ya=np.array(Sta_Pred_ya).reshape((4,10))
label=('Liquid','BCC','FCC','HCP')
for i in range(4):plt.bar(np.arange(10)+i*0.2,Sta_Pred_y[i],label=label[i],width=0.2)
plt.xticks(np.arange(10),('1.0','0.9','0.8','0.7','0.6','0.5','0.4','0.3','0.2','0.1'))
plt.xlabel('T*')
plt.legend()
plt.grid(1)
plt.title('AI_Class')
plt.savefig('Class_chart_MLP.jpg')
plt.close()
for i in range(4):plt.bar(np.arange(10)+i*0.2,Sta_Pred_ya[i],label=label[i],width=0.2)
plt.xticks(np.arange(10),('1.0','0.9','0.8','0.7','0.6','0.5','0.4','0.3','0.2','0.1'))
plt.xlabel('T*')
plt.legend()
plt.grid(1)
plt.title('W4W6_Class')
plt.savefig('Class_chart_newW4W6.jpg')
plt.close()

#%%
#Pointnetcate 
MLPcate = Pred_y[4]
Xcate = np.loadtxt('./AveXcate/XR'+'4'+'.txt',delimiter=' ')
W4W6cate = Pred_ya[4]
X_data = Pred_X[4]
Normalornot = 'true' #None 'true'
for i in (5,6,7,8,9):
#    MLPcate = np.concatenate((MLPcate,Pred_y[i]))
    Xcate = np.concatenate((Xcate,np.loadtxt('./AveXcate/XR'+(str)(i)+'.txt',delimiter=' ')))
#    W4W6cate = np.concatenate((W4W6cate,Pred_ya[i]))
#    X_data = np.concatenate((X_data,Pred_X[i]))
cnf_matrix = confusion_matrix(W4W6cate,MLPcate,
                              normalize=Normalornot)
cm_display = ConfusionMatrixDisplay(cnf_matrix,display_labels=
                                    ('Liquid','BCC','FCC','HCP')).plot()
#plt.show()
plt.savefig('Conf_Solid_NewW4W6_to_MLP1')
plt.close()
cnf_matrix = confusion_matrix(W4W6cate,Xcate,
                              normalize=Normalornot)
cm_display = ConfusionMatrixDisplay(cnf_matrix,display_labels=
                                    ('Liquid','BCC','FCC','HCP')).plot()
#plt.show()
plt.savefig('Conf_Solid_NewW4W6_to_Xcate1')
plt.close()
cnf_matrix = confusion_matrix(MLPcate,Xcate,
                              normalize=Normalornot)
cm_display = ConfusionMatrixDisplay(cnf_matrix,display_labels=
                                    ('Liquid','BCC','FCC','HCP')).plot()
#plt.show()
plt.savefig('Conf_Solid_MLP_to_Xcate1')
plt.close()
#%%
Liquiderror_Both= X_data[np.intersect1d(np.where(W4W6cate==0),np.where(MLPcate==0)),:]
Liquiderror_MLP = X_data[np.setdiff1d(np.where(W4W6cate==0),np.where(MLPcate==0)),:]
fig = plt.figure(figsize=[10.,10.])
gridplt = fig.add_subplot(111, projection='3d')
gridplt.scatter(Liquiderror_Both[:,1],Liquiderror_Both[:,2],Liquiderror_Both[:,0],c='b',s=10,alpha=0.8,label='Both_L')
gridplt.scatter(Liquiderror_MLP[:,1],Liquiderror_MLP[:,2],Liquiderror_MLP[:,0],c='r',s=10,alpha=0.8,label='MLP_Solid')
gridplt.set_xlabel('W4')
gridplt.set_ylabel('W6')
gridplt.set_zlabel('Connect')
gridplt.legend(loc='best', fontsize='large')
gridplt.view_init(elev=5,azim=0)
#plt.show(gridplt)
plt.savefig('ErrorScatter_W4W6toMLP.jpg')
plt.close()
#%%
#BCCerror_Both = X_data[np.intersect1d(np.where(MLPcate==1),np.where(Xcate==1),:]
BCCerror_X = X_data[np.setdiff1d(np.where(MLPcate==1),np.where(Xcate==1)),:]
BCCerror_Both = X_data[np.setdiff1d(np.where(MLPcate==0),np.where(Xcate==1)),:]
#BCCerror_X=np.where(Xcate==1)
fig = plt.figure(figsize=[10.,10.])
gridplt = fig.add_subplot(111, projection='3d')
gridplt.scatter(BCCerror_Both[:,1],BCCerror_Both[:,2],BCCerror_Both[:,0],c='b',s=10,alpha=0.8,label='MLP_Liquid')
gridplt.scatter(BCCerror_X[:,1],BCCerror_X[:,2],BCCerror_X[:,0],c='r',s=10,alpha=0.8,label='MLP_BCC')
gridplt.set_xlabel('W4')
gridplt.set_ylabel('W6')
gridplt.set_zlabel('Connect')
gridplt.legend(loc='best', fontsize='large')
gridplt.view_init(elev=90,azim=270)
#plt.show(gridplt)
plt.savefig('Scatter_MLPtoXN.jpg')
plt.close()
fig = plt.figure(figsize=[10.,10.])
gridplt = fig.add_subplot(111, projection='3d')
gridplt.scatter(BCCerror_Both[:,1],BCCerror_Both[:,2],BCCerror_Both[:,0],c='b',s=10,alpha=0.8,label='Both_Liquid')
gridplt.scatter(BCCerror_X[:,1],BCCerror_X[:,2],BCCerror_X[:,0],c='r',s=10,alpha=0.8,label='MLP_BCC')
gridplt.set_xlabel('W4')
gridplt.set_ylabel('W6')
gridplt.set_zlabel('Connect')
gridplt.legend(loc='best', fontsize='large')
gridplt.view_init(elev=5,azim=270)
#plt.show(gridplt)
plt.savefig('Scatter_MLPtoX_1N.jpg')
plt.close()