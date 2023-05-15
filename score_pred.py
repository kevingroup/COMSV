from sklearn import svm
import pandas as pd
import numpy as np
import sys

df = pd.read_csv("all_sample_with_label.txt")
df2 = pd.read_csv(sys.argv[1],sep="\t")
df2["Size"] = (df2["Size"]-df["Size2"].min()) / (df["Size2"].max()-df["Size2"].min())
df2["Coverage"] = (df2["Coverage"]-df["Coverage"].min()) / (df["Coverage"].max()-df["Coverage"].min())
df2.Coverage[df2.Coverage<1/(df["Coverage"].max()-df["Coverage"].min())] = 1/(df["Coverage"].max()-df["Coverage"].min())
df2.Coverage[df2.Coverage>1] = 1
df2.Size[df2.Size<1/(df["Size2"].max()-df["Size2"].min())] = 1/(df["Size2"].max()-df["Size2"].min())
df2.Size[df2.Size>1] = 1
X1=np.array(df2)

df["Size2"] = (df["Size2"]-df["Size2"].min()) / (df["Size2"].max()-df["Size2"].min())
df["Support"] = df["Support"] / df["Coverage"] # change to support ratio
df["Coverage"] = (df["Coverage"]-df["Coverage"].min()) / (df["Coverage"].max()-df["Coverage"].min())

X=np.array(df.iloc[:,[7,9,11]])
Y=np.array(df.iloc[:,-1])

clf = svm.SVC(kernel='rbf', probability=True, random_state = np.random.RandomState(42), gamma='auto')

Y1 = clf.fit(X, Y).predict_proba(X1)

Y1[:,1] = Y1[:,1]+Y1[:,2]
pd.DataFrame(Y1[:,1:],columns=['Score1','Score2']).to_csv(sys.argv[2],index=None,sep='\t')
