import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import KFold
from sklearn.linear_model import LogisticRegression
from sklearn.svm import l1_min_c

kf = KFold(n_splits=5, shuffle=True, random_state=123)
logreg = LogisticRegressionCV(cv=kf, solver="liblinear", penalty="l1", n_jobs=-1)

df = pd.read_csv("data/abn_lrads_lag_prescr.csv")

df = df[2::3]

# df[2::3].isna().sum(axis=1)
# all(df.pid.value_counts() == 1)
df = df.dropna(thresh=165, axis=0)
df.isna().sum()
col_list = df.columns.to_list()
# col_list.index('any_nodule')
col_list.index('LRcat')
feature_list = col_list[36:77]
X = df[feature_list].select_dtypes(include="number")
y = df.case

logreg.coefs_paths_[1][0][0]

cs = l1_min_c(X, y, loss='log') * np.logspace(0, 7, 16)

clf = LogisticRegression(penalty='l1', solver='saga',
                         tol=1e-6, max_iter=int(1e6),
                         warm_start=True)
coefs_ = []
for c in cs:
    clf.set_params(C=c)
    clf.fit(X, y)
    coefs_.append(clf.coef_.ravel().copy())

coefs = np.array(coefs_)
plt.plot(np.log10(cs), coefs, marker='o')
ymin, ymax = plt.ylim()
plt.xlabel('log(C)')
plt.ylabel('Coefficients')
plt.title('Logistic Regression Path')
plt.axis('tight')
plt.show()
