import numpy as np
import glob
from sklearn.model_selection import train_test_split
from keras.utils import np_utils
from tensorflow.keras import models
from tensorflow.keras import layers
from tensorflow.keras.callbacks import TensorBoard
import tensorflow as tf
config = tf.compat.v1.ConfigProto(gpu_options=tf.compat.v1.GPUOptions(allow_growth=True))
sess = tf.compat.v1.Session(config=config)
NearestP=16 #16,43,83
rnd = np.random.RandomState(223)
#%% Data preparation
def Get_dataset(Pre_file):    
    Pre_list=glob.glob(Pre_file)
    Pre_list.sort()
    Data = np.loadtxt(Pre_list[0],delimiter=' ')
    for i in range(1,len(Pre_list)):Data = np.concatenate((Data,np.loadtxt(Pre_list[i], delimiter=' ')))
    return Data
def Rotation(Data):
    for i in range(len(Data)):
        Rot_mat=RMatrix(rnd.uniform(0,360),rnd.uniform(0,180),rnd.uniform(0,360))
        Data[i]=np.transpose(np.dot(Rot_mat,np.transpose(Data[i])))
    return Data
def RMatrix(Alpha=0,Beta=0,Gammar=0):    
    Theta_0,Theta_1,Theta_2=np.pi/180*Alpha,np.pi/180*Beta,np.pi/180*Gammar
    Rot_mat=np.array(((np.cos(Theta_0)*np.cos(Theta_2)-np.cos(Theta_1)*np.sin(Theta_0)*np.sin(Theta_2),
                       np.sin(Theta_0)*np.cos(Theta_2)+np.cos(Theta_1)*np.cos(Theta_0)*np.sin(Theta_2),
                       np.sin(Theta_1)*np.sin(Theta_2)),
                      (-np.cos(Theta_0)*np.sin(Theta_2)-np.cos(Theta_1)*np.sin(Theta_0)*np.cos(Theta_2),
                       -np.sin(Theta_0)*np.sin(Theta_2)+np.cos(Theta_1)*np.cos(Theta_0)*np.cos(Theta_2),
                       np.sin(Theta_1)*np.cos(Theta_2)),
                      (np.sin(Theta_1)*np.sin(Theta_0),-np.sin(Theta_1)*np.cos(Theta_0),np.cos(Theta_1))))
    return Rot_mat
def Split_Data(Data,Answer):
    Data = Data.reshape((-1,NearestP,3,1))
    Train_X,Test_X,Train_y,Test_y = train_test_split(Data,Answer,test_size=0.05,random_state=0,stratify=Answer)
    Train_y = np_utils.to_categorical(Train_y, num_classes=4)
    return Train_X,Test_X,Train_y,Test_y
Train_X = Get_dataset(r'./Trainset/NB*').reshape((-1,NearestP,3))
Train_X = Rotation(Train_X)
Train_y = Get_dataset(r'./Trainset/Crystal*')
Train_X,Test_X,Train_y,Test_ya=Split_Data(Train_X,Train_y)
np.save(r'./Trainset/Train_X.npy',Train_X)
np.save(r'./Trainset/Train_y.npy',Train_y)
np.save(r'./Trainset/Test_X.npy',Test_X)
np.save(r'./Trainset/Test_ya.npy',Test_ya)
#%% Model build
def PointNet(Train_X,Train_y,Train_epochs=50):
    ''' Classification PointNet, input is BxNx3, output Bxn where n is num classes '''
    model = models.Sequential()
    model.add(layers.Convolution2D(filters=64,kernel_size=[1,3],strides=[1,1]
                ,padding='valid',name='Conv1',activation='relu'
                ,batch_input_shape=(None,Train_X.shape[1],Train_X.shape[2],1)))
    model.add(layers.BatchNormalization(name='Bn1'))
    model.add(layers.Convolution2D(filters=64,kernel_size=[1,1],strides=[1,1]
                ,padding='valid',activation='relu',name='Conv2'))
    model.add(layers.BatchNormalization(name='Bn2'))
    model.add(layers.Convolution2D(filters=64,kernel_size=[1,1],strides=[1,1]
                ,padding='valid',activation='relu',name='Conv3'))
    model.add(layers.BatchNormalization(name='Bn3'))
    model.add(layers.Convolution2D(filters=128,kernel_size=[1,1],strides=[1,1]
                ,padding='valid',activation='relu',name='Conv4'))
    model.add(layers.BatchNormalization(name='Bn4'))  
    model.add(layers.Convolution2D(filters=1024,kernel_size=[1,1],strides=[1,1]
                ,padding='valid',activation='relu',name='Conv5'))                                
    model.add(layers.BatchNormalization(name='Bn5'))
    model.add(layers.MaxPooling2D(pool_size=[NearestP,1],padding='valid',name='Pool'))
    model.add(layers.Flatten(name='Flat'))
    model.add(layers.Dense(units=512,activation='relu',name='Den1'))
    model.add(layers.BatchNormalization(name='Bn6'))
    model.add(layers.Dense(units=256,activation='relu',name='Den2'))
    model.add(layers.BatchNormalization(name='Bn7'))
    model.add(layers.Dropout(rate=0.5,name='Drop'))
    model.add(layers.Dense(units=4,activation='softmax',name='Den3'))
    model.compile(optimizer='adam',loss='categorical_crossentropy',metrics=['accuracy'])
    model.summary()
    tbCallBack = TensorBoard(log_dir=r'./Tensorboard',update_freq='epoch'
                             ,write_graph=True,write_images=False,)
    model.fit(Train_X,Train_y,validation_split=0.2,epochs=Train_epochs,batch_size=128,callbacks=[tbCallBack])
    model.save(r'./classic_modeln12.h5')
    return model
model=PointNet(Train_X,Train_y,100)
Test_y=np.argmax(model.predict(Test_X),axis=1)
np.save("./Trainset/Test_y.npy",Test_y)
#%% save result
Pred_X=Get_dataset(r'./Predictset/NB*.txt').reshape((-1,NearestP,3,1))
np.save(r'./Predictset/Pred_X.npy',Pred_X)
Pred_PointNet=np.argmax(model.predict(Pred_X),axis=1)
np.save(r'./Predictset/Pred_PointNet.npy',Pred_PointNet)

Pred_W4W6=Get_dataset(r'./W4W6/TanaR*.txt')
np.save(r'./W4W6/Pred_W4W6.npy',Pred_W4W6)
Pred_NewW4W6=Get_dataset(r'./New_W4W6/TanaR*.txt')
np.save(r'./New_W4W6/Pred_NewW4W6.npy',Pred_NewW4W6)

Pred_Xcate=Get_dataset(r'./Xcate/XR*.txt')
np.save(r'./Xcate/Pred_Xcate.npy',Pred_Xcate)
Pred_AveXcate=Get_dataset(r'./AveXcate/XR*.txt')
np.save(r'./AveXcate/Pred_AveXcate.npy',Pred_AveXcate)

Pred_CrystalX=Get_dataset(r'./Crystal/NBR*.txt').reshape((-1,NearestP,3,1))
Pred_CrystalY=np.argmax(model.predict(Pred_CrystalX),axis=1)
Pred_CrystalX=Get_dataset(r'./Crystal/NB0*.txt').reshape((-1,NearestP,3,1))
Pred_CrystalY=np.argmax(model.predict(Pred_CrystalX),axis=1)
np.save("./Crystal/Pred_Crystal0.npy",Pred_CrystalY)
Pred_CrystalX=Get_dataset(r'./Crystal/NB1*.txt').reshape((-1,NearestP,3,1))
Pred_CrystalY=np.argmax(model.predict(Pred_CrystalX),axis=1)
np.save("./Crystal/Pred_Crystal1.npy",Pred_CrystalY)
Pred_CrystalX=Get_dataset(r'./Crystal/NB2*.txt').reshape((-1,NearestP,3,1))
Pred_CrystalY=np.argmax(model.predict(Pred_CrystalX),axis=1)
np.save("./Crystal/Pred_Crystal2.npy",Pred_CrystalY)
Pred_CrystalX=Get_dataset(r'./Crystal/NB3*.txt').reshape((-1,NearestP,3,1))
Pred_CrystalY=np.argmax(model.predict(Pred_CrystalX),axis=1)
np.save("./Crystal/Pred_Crystal3.npy",Pred_CrystalY)
