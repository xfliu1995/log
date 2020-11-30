# 3.Machine Learning with Python

本章介绍如何用Python进行简单的机器学习，练习时可以在R中安装一下如下文所示的一些Python packages即可。

也可以使用我们提供的docker（下载链接如[附表](../appendix/appendix-iv.-teaching.md#teaching-docker)所示），里面已经安装好了很多Python packages：

```bash
docker load -i ~/Desktop/bioinfo_pca_machine.tar.gz
docker run --name=bioinfo_pca_machine -dt -h bioinfo_docker --restart unless-stopped -v ~/Downloads/data:/data gangxu/machine_learning:2.0
docker exec -it bioinfo_pca_machine bash
cd /home/test/machine_python/
```

## 1\) Practice Guide

### 1a\) 本章教程使用指南

读者将会发现，机器学习的核心模型已经被**scikit-learn**等工具包非常好的模块化了，调用起来非常简单，仅需要几行代码，但是一个完整的、有效的机器学习工程项目却包括很多步骤，可以包括**数据导入，数据可视化理解，前处理，特征选择，模型训练，参数调整，模型预测，模型评估，后处理**等多个步骤，一个在真实世界中有效的模型可能需要工作者对数据的深入理解，以选择各个步骤合适的方法。

通过本章教程，读者可以对机器学习的基本概念方法和具体流程有所了解，而且可以通过实践更好地掌握python相关工具包的使用，为后续的应用做好准备。

读者初次阅读和进行代码实践时，可以将重点放在对方法和概念的理解上，对于一些稍微复杂的代码，不需要理解代码里的每个细节。

### 1b\) 导入需要的Python工具包

这里我们会导入一些后续操作需要的python工具包，它们的相关文档如下，请有兴趣的读者重点学习和了解[scikit-learn](http://scikit-learn.org/)工具包。

* [numpy](https://docs.scipy.org/doc/numpy/): arrays
* [pandas](https://pandas.pydata.org/): data IO, DataFrame
* [scikit-learn](http://scikit-learn.org/): machine learning
* [statsmodels](https://www.statsmodels.org/): statistical functions
* [matplotlib](https://matplotlib.org/): plotting
* [seaborn](https://matplotlib.org/): high-level plotting based on _matplotlib_
* [jupyter](https://jupyter.org/): Python notebook

> tips: Anaconda并没有集成seaborn的最新版本，请使用pip更新seaborn：`pip install seaborn==0.9.0`，然后使用jupyter notebook新建文件，运行下面的代码

```python
# For data importing
import pandas as pd
# For machine learning
from sklearn.datasets import make_classification
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import KFold, train_test_split
from sklearn.metrics import accuracy_score, roc_auc_score, f1_score, recall_score, precision_score, \
    roc_curve, precision_recall_curve, average_precision_score, matthews_corrcoef, confusion_matrix
# For plotting
import seaborn as sns
sns.set()
sns.set_style('whitegrid')
```

### 1c\) 产生数据集

在处理真实世界的数据集之前，我们先产生一些模拟的数据集来学习机器学习的基本概念。 _scikit-learn_ 提供了很多方法\([sklearn.datasets](http://scikit-learn.org/stable/modules/classes.html#module-sklearn.datasets)\) 来方便地产生数据集。

我们可以产生一个标签为离散值的用于分类问题的数据集:

[sklearn.datasets.make\_classification](http://scikit-learn.org/stable/modules/generated/sklearn.datasets.make_classification.html#sklearn.datasets.make_classification) 可以从一个混合高斯分布中产生样本，并且可以控制样本数量，类别数量和特征数量。

我们会产生一个数据集，共有1000个样本，两种类别，四种特征。本章教程使用该数据作为演示。

* **产生数据:**

```python
random_state = np.random.RandomState(1289237)  #我们在本教程中固定numpy的随机种子，以使结果可重现
X, y = make_classification(n_samples=1000, n_classes=2, n_features=4,
                           n_informative=2, n_redundant=0, n_clusters_per_class=1,
                           class_sep=0.9, random_state=random_state)
X.shape, y.shape #查看特征和标签的shape
```

* **用matplotlib可视化样本数据的分布:**

```python
fig, ax = plt.subplots(figsize=(7, 7))
for label in np.unique(y):
    ax.scatter(X[y == label, 0], X[y == label, 1], s=10, label=str(label))
ax.legend(title='Class')
plot.show()
```

![](../.gitbook/assets/1.simple-machine-learning-basics_14_1.png)

### 1d\) 数据预处理

* **使用standard/z-score scaling 对数据做scaling:**

```text
X = StandardScaler().fit_transform(X)
```

这里我们展示了归一化前后数据变化

```python
#产生模拟数据，1000个数据点，均值为10，标准差为2
random_state = np.random.RandomState(1289237)
x = random_state.normal(10, 2, size=1000)
fig, ax = plt.subplots(1,2,figsize=(16, 6))
sns.distplot(x, ax=ax[0])
sns.distplot(x, ax=ax[1])
sns.distplot(np.ravel(x), ax=ax[0])
sns.distplot(np.ravel(StandardScaler().fit_transform(x.reshape((-1, 1)))), ax=ax[1])
ax[0].set_title('original data distribution',fontsize=20)
ax[1].set_title('scaled data distribution by standard scaling',fontsize=20)
plot.show()
```

![](../.gitbook/assets/1.simple-machine-learning-basics_37_1.png)

### 1e\) 划分数据得到训练集和测试集（training set & test set）

之前我们已经对数据进行了一些分析，并且做了一些基本的预处理，接下来我们需要对数据进行划分，得到训练集、验证集、测试集。

* 训练集：用于训练模型；
* 验证集：用于确定模型参数；
* 测试集：不参与模型训练过程，仅用于评价模型效果。

因为模型总是会在某种程度上过拟合训练数据，因此在训练数据上评估模型是有偏的，模型在训练集上的表现总会比测试集上好一些。

因为模型总是可以学到数据中隐藏的模式和分布，如果样本间彼此的差异比较大，过拟合问题就会得到一定程度的减轻。而如果数据的量比较大，模型在训练集和测试集上的表现差异就会减小。

**划分训练集和测试集** 这里我们使用[train\_test\_split](http://scikit-learn.org/stable/modules/generated/sklearn.model_selection.train_test_split.html) 方法来随机的将80%的样本设置为训练样本， 将其余20%设置为测试样本。

```python
random_state = np.random.RandomState(1289237)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=random_state)
print('number of training samples: {}, test samples: {}'.format(X_train.shape[0], X_test.shape[0]))
```

```text
number of training samples: 800, test samples: 200
```

**划分训练集和验证集** 验证集（validation set），是将训练集再随机划分为训练集和验证集，进行多折交叉验证（[cross validation](https://www.zhihu.com/question/39259296)），可以帮助我们评估不同的模型，调整模型的超参数等，此外交叉验证在数据集较小的时候也被用于直接评估模型的表现，我们在交叉验证部分还会详细讲解。

### 1f\) 使用机器学习模型进行分类

* **示例：调用逻辑回归模型并且训练模型**

> 逻辑斯谛回归是一个简单但是非常有效的模型，与它的名字不同，逻辑斯谛回归用于解决分类问题，在二分类问题中被广泛使用。对于二分类问题，我们需要对每一个样本预测它属于哪一类（0或者1）。
>
> 逻辑斯谛回归是一个线性分类模型，它会对输入的feature进行线性组合，然后将线性组合组合得到的值通过一个非线性的sigmoid函数映射为一个概率值\(范围为0~1\)。
>
> 模型训练过程中，模型内部的参数（线性模型的权重）会调整使得模型的损失函数（真实label和预测label的交叉熵）最小。
>
> $$p(y_i | \mathbf{x}_i) = \frac{1}{1 + \text{exp} \left( \sum_{j=1}^M x_{ij} w_{j} + b \right)}$$

使用sklearn封装好的模型进行模型的训练非常简单，以逻辑斯谛回归模型为例，只需要两行即可完成模型的训练，我们使用默认参数构建模型。

```python
model = LogisticRegression()
model
```

```text
LogisticRegression(C=1.0, class_weight=None, dual=False, fit_intercept=True,
          intercept_scaling=1, max_iter=100, multi_class='ovr', n_jobs=1,
          penalty='l2', random_state=None, solver='liblinear', tol=0.0001,
          verbose=0, warm_start=False)
```

### 1g）训练集上的交叉验证

我们首先在**训练集**上做**K折（k-folds）交叉验证**，在训练集上划分出一部分用于**训练**，另一部分用于**验证**，可以帮助我们挑选比较不同的模型，以及挑选模型中的超参数。

scikit-learn\_提供很多功能来[划分数据集](http://scikit-learn.org/stable/modules/classes.html#module-sklearn.model_selection).

这里我们使用[KFold](http://scikit-learn.org/stable/modules/generated/sklearn.model_selection.KFold.html) 来将_训练集_划分为10折，5和10是交叉验证中经常使用的折数。如果样本数量和计算资源允许，一般设置为10折。

下面的代码展示_KFold_是如何划分数据集的，图片中每一行为一个轮次，每一行中黑色的box为该轮次的测试集

```python
n_splits = 10

kfold = KFold(n_splits=n_splits, random_state=random_state)
is_train = np.zeros((n_splits, X_train.shape[0]), dtype=np.bool)
for i, (train_index, test_index) in enumerate(kfold.split(X_train, y_train)):
    is_train[i, train_index] = 1

fig, ax = plt.subplots(figsize=(15, 3))
ax.pcolormesh(is_train)
ax.set_yticks(np.arange(n_splits) + 0.5)
ax.set_yticklabels(np.arange(n_splits) + 1)
ax.set_ylabel('Round')
ax.set_xlabel('Sample')
plot.show()
```

![](../.gitbook/assets/1.simple-machine-learning-basics_31_1.png)

**利用10折交叉划分数据** 接下来我们在训练集上训练模型，对验证集进行预测，这样我们可以分析模型在10折交叉验证中每一轮时在训练集和验证集分别的表现。

```python
predictions = np.zeros((n_splits, X_train.shape[0]), dtype=np.int32)
predicted_scores = np.zeros((n_splits, X_train.shape[0]))

for i in range(n_splits):
    model.fit(X_train[is_train[i]], y_train[is_train[i]])
    predictions[i] = model.predict(X_train)
    predicted_scores[i] = model.predict_proba(X_train)[:, 1]
```

#### 10折交叉验证的评价指标

我们统计了模型10折交叉验证的指标

```python
scorers = {'accuracy': accuracy_score,
           'recall': recall_score,
           'precision': precision_score,
           'f1': f1_score,
           'mcc': matthews_corrcoef
}
cv_metrics = pd.DataFrame(np.zeros((n_splits*2, len(scorers) + 2)),
                          columns=list(scorers.keys()) + ['roc_auc', 'average_precision'])
cv_metrics.loc[:, 'dataset'] = np.empty(n_splits*2, dtype='U')
for i in range(n_splits):
    for metric in scorers.keys():
        cv_metrics.loc[i*2 + 0, metric] = scorers[metric](y_train[is_train[i]], predictions[i, is_train[i]])
        cv_metrics.loc[i*2 + 1, metric] = scorers[metric](y_train[~is_train[i]], predictions[i, ~is_train[i]])
    cv_metrics.loc[i*2 + 0, 'roc_auc'] = roc_auc_score(y_train[is_train[i]], predicted_scores[i, is_train[i]])
    cv_metrics.loc[i*2 + 1, 'roc_auc'] = roc_auc_score(y_train[~is_train[i]], predicted_scores[i, ~is_train[i]])
    cv_metrics.loc[i*2 + 0, 'average_precision'] = average_precision_score(y_train[is_train[i]], 
                                                                           predicted_scores[i, is_train[i]])
    cv_metrics.loc[i*2 + 1, 'average_precision'] = average_precision_score(y_train[~is_train[i]], 
                                                                           predicted_scores[i, ~is_train[i]])
    cv_metrics.loc[i*2 + 0, 'dataset'] = 'train'
    cv_metrics.loc[i*2 + 1, 'dataset'] = 'valid'

cv_metrics
```

|  | accuracy | recall | precision | f1 | mcc | roc\_auc | average\_precision | dataset |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |


| 0 | 0.883333 | 0.907609 | 0.869792 | 0.888298 | 0.767081 | 0.944641 | 0.936965 | train |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |


| 1 | 0.837500 | 0.894737 | 0.790698 | 0.839506 | 0.681520 | 0.897870 | 0.810110 | valid |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |


| 2 | 0.873611 | 0.901370 | 0.856771 | 0.878505 | 0.748032 | 0.936600 | 0.912893 | train |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |


| 3 | 0.925000 | 0.951220 | 0.906977 | 0.928571 | 0.850786 | 0.972483 | 0.973525 | valid |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |


| 4 | 0.869444 | 0.898352 | 0.851562 | 0.874332 | 0.739840 | 0.935370 | 0.912378 | train |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |


| 5 | 0.962500 | 0.952381 | 0.975610 | 0.963855 | 0.925196 | 0.981830 | 0.979835 | valid |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |


| 6 | 0.876389 | 0.907162 | 0.863636 | 0.884864 | 0.752664 | 0.934445 | 0.915088 | train |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |


| 7 | 0.900000 | 0.896552 | 0.838710 | 0.866667 | 0.787929 | 0.977688 | 0.960846 | valid |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |


| 8 | 0.879167 | 0.910811 | 0.861893 | 0.885677 | 0.759053 | 0.940394 | 0.919333 | train |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |


| 9 | 0.850000 | 0.861111 | 0.815789 | 0.837838 | 0.699376 | 0.935606 | 0.919198 | valid |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |


| 10 | 0.879167 | 0.909091 | 0.859375 | 0.883534 | 0.759494 | 0.935921 | 0.911489 | train |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |


| 11 | 0.887500 | 0.883721 | 0.904762 | 0.894118 | 0.774397 | 0.976744 | 0.981484 | valid |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |


| 12 | 0.873611 | 0.901370 | 0.856771 | 0.878505 | 0.748032 | 0.941486 | 0.922226 | train |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |


| 13 | 0.875000 | 0.926829 | 0.844444 | 0.883721 | 0.753015 | 0.922452 | 0.891856 | valid |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |


| 14 | 0.881944 | 0.909341 | 0.864230 | 0.886212 | 0.764789 | 0.943882 | 0.921491 | train |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |


| 15 | 0.837500 | 0.880952 | 0.822222 | 0.850575 | 0.674881 | 0.894737 | 0.889792 | valid |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |


| 16 | 0.877778 | 0.898305 | 0.859459 | 0.878453 | 0.756415 | 0.945255 | 0.918845 | train |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |


| 17 | 0.850000 | 0.865385 | 0.900000 | 0.882353 | 0.676665 | 0.888049 | 0.917479 | valid |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |


| 18 | 0.883333 | 0.906593 | 0.868421 | 0.887097 | 0.767282 | 0.940486 | 0.916185 | train |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |


| 19 | 0.825000 | 0.904762 | 0.791667 | 0.844444 | 0.654015 | 0.925439 | 0.933022 | valid |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |


#### 10折交叉验证ROC

```text
from scipy import interp
cv_metrics_mean = cv_metrics.groupby('dataset').mean()
fig, axes = plt.subplots(1, 2, figsize=(14, 7))
# ROC curve
ax = axes[0]
all_fprs = np.linspace(0, 1, 100)
roc_curves = np.zeros((n_splits, len(all_fprs), 2))
for i in range(n_splits):
    fpr, tpr, thresholds = roc_curve(y_train[~is_train[i]], predicted_scores[i, ~is_train[i]])
    roc_curves[i, :, 0] = all_fprs
    roc_curves[i, :, 1] = interp(all_fprs, fpr, tpr)
roc_curves = pd.DataFrame(roc_curves.reshape((-1, 2)), columns=['fpr', 'tpr'])
sns.lineplot(x='fpr', y='tpr', data=roc_curves, ci='sd', ax=ax,
             label='Test AUC = {:.4f}'.format(cv_metrics_mean.loc['valid', 'roc_auc']))
#ax.plot(fpr, tpr, label='ROAUC = {:.4f}'.format(roc_auc_score(y_test, y_score[:, 1])))
#ax.plot([0, 1], [0, 1], linestyle='dashed')
ax.set_xlabel('False positive rate')
ax.set_ylabel('True positive rate')
ax.plot([0, 1], [0, 1], linestyle='dashed', color='gray')
ax.set_title('ROC curve')
ax.legend()
plot.show()

# predision-recall curve
ax = axes[1]
all_precs = np.linspace(0, 1, 100)
pr_curves = np.zeros((n_splits, len(all_precs), 2))
for i in range(n_splits):
    fpr, tpr, thresholds = precision_recall_curve(y_train[~is_train[i]], predicted_scores[i, ~is_train[i]])
    pr_curves[i, :, 0] = all_precs
    pr_curves[i, :, 1] = interp(all_precs, fpr, tpr)
pr_curves = pd.DataFrame(pr_curves.reshape((-1, 2)), columns=['precision', 'recall'])
sns.lineplot(x='precision', y='recall', data=pr_curves, ci='sd', ax=ax,
             label='Test AP = {:.4f}'.format(cv_metrics_mean.loc['valid', 'average_precision']))

ax.set_xlabel('Precision')
ax.set_ylabel('Recall')
ax.plot([0, 1], [1, 0], linestyle='dashed', color='gray')
ax.set_title('Precision-recall curve')
ax.legend()
plot.show()
```

![](../.gitbook/assets/1.simple-machine-learning-basics_78_1.png)

### 1h\) 在整个训练集\(training set\)上进行模型训练

同样使用**Logistic Regression**模型。

**利用10折交叉验证优化参数**

```text
tree_para = {'C': [0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000]}
model = GridSearchCV(LogisticRegression(penalty='l2',solver='liblinear'), tree_para, cv=10)
```

```text
model.fit(X_train, y_train)
```

### 1i\) 在测试集\(test set\)上预测和评估整个训练集\(traning set\)得到的模型

#### 在测试集上预测样本类别

为了评估模型表现，我们需要对测试集样本进行预测，我们使用_predict_方法来预测样本类别，它会返回一个整数型array来表示不同的样本类别。

```python
y_pred = model.predict(X_test)
y_pred
```

```text
array([1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1,
       1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0,
       0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1,
       0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1,
       0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0,
       0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1,
       1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1,
       0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0,
       1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0,
       1, 0])
```

#### **构建预测结果的Confusion matrix：**

使用scikit-learn的confusion\_matrix方法即可得到模型预测结果的confusion matrix

```python
pd.DataFrame(confusion_matrix(y_test, y_pred), 
             columns=pd.Series(['Negative', 'Positive'], name='Predicted'),
             index=pd.Series(['Negative', 'Positive'], name='True'))
```

| Predicted | Negative | Positive |
| :--- | :--- | :--- |
| True |  |  |
| Negative | 92 | 13 |
| Positive | 11 | 84 |

```python
scorers = {'accuracy': accuracy_score,
           'recall': recall_score,
           'precision': precision_score,
           'f1': f1_score,
           'mcc': matthews_corrcoef
}
for metric in scorers.keys():
    print('{} = {}'.format(metric, scorers[metric](y_test, y_pred)))
```

```text
accuracy = 0.88
recall = 0.8842105263157894
precision = 0.865979381443299
f1 = 0.8749999999999999
mcc = 0.7597918897600666
```

#### 绘制模型评估性能图

**绘制ROC曲线和Precision-Recall曲线**

我们使用sklearn自带的_roc\_curve_和_precision\_recall\_curve_方法来计算绘图需要的指标，这两个方法需要的输入为测试集每个样本的真实标签和模型预测的每个样本的概率。

```python
fig, axes = plt.subplots(1, 2, figsize=(14, 7))
# ROC curve
y_score = model.predict_proba(X_test)
fpr, tpr, thresholds = roc_curve(y_test, y_score[:, 1])
ax = axes[0]
ax.plot(fpr, tpr, label='AUROC = {:.4f}'.format(roc_auc_score(y_test, y_score[:, 1])))
ax.plot([0, 1], [0, 1], linestyle='dashed')
ax.set_xlabel('False positive rate')
ax.set_ylabel('True positive rate')
ax.set_title('ROC curve')
ax.legend()
plot.show()
# predision-recall curve
precision, recall, thresholds = precision_recall_curve(y_test, y_score[:, 1])
ax = axes[1]
ax.plot(precision, recall, label='AP = {:.4f}'.format(average_precision_score(y_test, y_score[:, 1])))
ax.plot([0, 1], [1, 0], linestyle='dashed')
ax.set_xlabel('Precision')
ax.set_ylabel('Recall')
ax.set_title('Precision-recall curve')
ax.legend()
plot.show()
```

![](../.gitbook/assets/1.simple-machine-learning-basics_64_1_1_1.png)

可以看到AUROC和AP都接近于1，可以认为模型的分类效果很好。

## 2\) Homework

按照教程中的流程，利用不同分类器模型对给定数据集BreastCancer进行二分类，并且汇报每个模型的预测表现，包括accuracy,precision,recall, roc\_auc等指标，并绘制ROC曲线，sensitivity-specificity曲线。

> * 数据：R中mlbench包中的数据集BreastCancer（可从利用提供[文件](https://github.com/lulab/teaching_book/tree/master/files/PART_III/7.clinical_analyses/PCA_tSNE/)）。数据集包括9个特征，两种类别，良性（benign）和恶性（malignant）。以下是各特征的含义：

| Id | Cl.thickness | Cell.size | Cell.shape | Marg.adhesion | Epith.c.size | Bare.nuclei | Bl.cromatin | Normal.nucleoli | Mitoses | Class |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 1000025 | 5 | 1 | 1 | 1 | 2 | 1 | 3 | 1 | 1 | benign |
| 1002945 | 5 | 4 | 4 | 5 | 7 | 10 | 3 | 2 | 1 | benign |
| 1015425 | 3 | 1 | 1 | 1 | 2 | 2 | 3 | 1 | 1 | benign |
| 1016277 | 6 | 8 | 8 | 1 | 3 | 4 | 3 | 7 | 1 | benign |
| 1017023 | 4 | 1 | 1 | 3 | 2 | 1 | 3 | 1 | 1 | benign |
| 1017122 | 8 | 10 | 10 | 8 | 7 | 10 | 9 | 7 | 1 | malignant |

Cl.thickness: Clump Thickness Cell.size: Uniformity of Cell Size Cell.shape: Uniformity of Cell Shape Marg.adhesion: Marginal Adhesion Epith.c.size: Single Epithelial Cell Size Bare.nuclei: Bare Nuclei Bl.cromatin: Bland Chromatin Normal.nucleoli: Normal Nucleoli Mitoses: Mitoses 如需了解更多关于BreastCancer数据集信息，可参考[mlbench](https://cran.r-project.org/web/packages/mlbench/index.html)的文档。

> * 数据预处理：1）去除含有空缺值的样本 2）对数据进行归一化
> * 数据标签处理：正样本为benign，负样本为malignant
> * 数据集划分：训练集和测试集划分参考教程中的80%/20%划分方式，程序运行最开头加上random\_state = np.random.RandomState\(1289237\)保证划分一致
> * 分类器模型：LR,SVM,DT,RF
> * 编程工具：R/python
> * 作业要求：上传word/pdf文档附件，记录处理过程所用代码，并汇报每个模型在测试集的accuracy,precision,recall, roc\_auc等指标，并绘制ROC曲线，sensitivity-specificity曲线。

