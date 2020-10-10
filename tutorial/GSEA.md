# 3.3.GSEA

## 1\) Pipeline

![a8cf1bf98ac9bc82f7cae49abb771286.png](en-resource://database/3320:1)


## 2\) Getting software & data 

### \(1\) 程序下载与安装

**GSEA**(gene set enrichment analysis)即基因集富集分析，是常见富集分析之一。给定一个排序的基因表L和一个预先定义的基因集S (比如编码某个代谢通路的产物的基因, 基因组上物理位置相近的基因，或同一GO注释下的基因)，GSEA的目的是判断S里面的成员s在L里面是随机分布还是主要聚集在L的顶部或底部。

本文基于GSEA的图形化界面进行操作。
* 进入 [GSEA官方下载页面](http://software.broadinstitute.org/gsea/downloads.jsp)（需要登录方可进入，只需输入邮箱即可），可在安装列表中查找电脑对应版本安装 ；

### \(2\) 准备数据
#### \(2a\) Input format

| Format | Description | Notes |
| :--- | :--- | :--- |
| txt, res, gct, pcl | A _dataset_ containing _expression_ value for each feature in each sample. | 658个基因的表达值（所有基因的表达值） |
| cls | Associates each sample with a _phenotype label_. | 基因在不同样品下的表达值 |
| Chip | Lists each _probe_ and its matching HUGO _gene symbol_. | 基因 ID 和基因名字对应表 |
| gmx or gmt | gives the _gene set_ name and list of features in it. | 特定功能的基因集合 |

参见 [http://software.broadinstitute.org/gsea/doc/GSEAUserGuideTEXT.htm\#\_Preparing\_Data\_Files](http://software.broadinstitute.org/gsea/doc/GSEAUserGuideTEXT.htm#_Preparing_Data_Files)

#### \(2b\) Our input files

这里我们分析  [GSE19161](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19161) 数据集中与 [_EIF4G2_](https://www.ncbi.nlm.nih.gov/gene/1982) 基因的表达量显著相关的 gene set。其中输入文件包括所有基因的表达矩阵，EIF4G2在不同样本下的表达值，基因 ID 和基因名字对应表，以及特定功能的基因集合（可通过GSEA自带的数据库）。

* 我们已经准备好了文件，参考[文件获取方式]()，直接下载GSE19161.txt（表达矩阵），GSE19161.cls（EIF4G2在不同样本下的表达值），GSE19161.chip（基因 ID 和基因名字对应表）。

* 也可以通过qGSEA(R package)在 [GEO](https://www.ncbi.nlm.nih.gov/geo/browse/?view=series&display=20) 中远程下载感兴趣的数据集（本例中我们使用 [GSE19161](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19161)。

```bash
R
# 进入 R 操作界面
> library(qGSEA)
> dir.create('raw/')
> download.file(rGEO::gse_soft_ftp('GSE19161'), 'raw/GSE19161_family.soft.gz')
> download.file(rGEO::gse_matrix_ftp('GSE19161'), 'raw/GSE19161_series_matrix.txt.gz')

> dir.create('input/')
> qGSEA::make_gsea_input(
matrix_file = 'raw/GSE19161_series_matrix.txt.gz',
soft_file = 'raw/GSE19161_family.soft.gz',
output_dir = 'input/',
gene = 'EIF4G2'
)
> q()
Save workspace image? [y/n/c]: n # 按 n 再按 Enter
```

* 可得到文件目录结构如下：

```text
gsea
├── input
│   ├── GSE19161.chip
│   ├── GSE19161.cls
│   └── GSE19161.txt
└── raw
    ├── GSE19161_family.soft.gz
    └── GSE19161_series_matrix.txt.gz
```


## 3\) Run GSEA in a graphic interface

> **注意：由于 Docker 中无法启动图形化界面，本步骤需要在自己的电脑中进行。**

 这里我们分析 GSE19161 数据集中与 [_EIF4G2_](https://www.ncbi.nlm.nih.gov/gene/1982) 基因的表达量显著相关的 gene set。

1）提前准备好软件和数据（见 [getting software & data]()），[准备好输入文件]()。

2）打开 GSEA，点击**Load data**，选择文件夹下相应文件 GSE19161.txt，GSE19161.cls，GSE19161.chip。

![a72da496f267941d032374166051d526.gif](en-resource://database/3322:1)

3）点击**Run Gsea**。首先选择合适的参数，如下图所示。
* **Expression dataset**选择上传生成的**GSE19161**；
* **Gene sets database**选择**h.all.v6.1.symbols.gmt**([hallmark gene set collection](https://doi.org/10.1016/j.cels.2015.12.004))；
* **Phenoype labels**选择**GSE19161.cls#EIF4G2**；
* **Number of permutations**选择**1000**。置换检验的次数，数字越大结果越准确，但是太大会占用太多内存，软件默认检验1000次。软件分析时会得到一个基因富集的评分（ES），但是富集评分是否具有统计学意义，软件就会采用随机模拟的方法，根据指定参数随机打乱1000次，得到1000个富集评分，然后判断得到的ES是否在这1000个随机产生的得分中有统计学意义。测试使用时建议填一个很小的数如10，先让程序跑通。真正分析时再换为1000。
* **Collapse dataset to gene symbols**选择**True**。 如果表达数据集文件中NAME已经与gene sets database中名字一致，选择FALSE，反之选择TRUE。
* **Chip platform**选择**GSE19161.chip**；
* 点击**Basci fields**，在**Metric for ranking genes**选择**Pearson**（由于我们使用的是连续性表型（EIF4G2 基因的表达量）。**Metric for ranking genes**代表基因排序的度量。表型是**分组信息**，GSEA在计算分组间的差异值时支持5种统计方式，分别是signal2noise、t-Test、ratio_of_class、 diff_of_class(log2转换后的值计算倍数)和log2_ratio_of_class；表型是**连续数值信息**，使用pearson相关性、Cosine、Manhattan 或Euclidean指标之一计算两个配置文件之间的相关性。
* 其他为默认参数。

![57ea0979ac246c6da1e985cd379973da.gif](en-resource://database/3324:1)

4) 最后点击**Run**。

5)运行成功后，`Status` 会显示为 `Success`，点击其即可查看输出。在本例中，没有一个 gene set 通过了统计上的显著性检验\(FDR&lt;0.05\)。


![b8342e01af3762b7924e0d9c7d6e87dc.gif](en-resource://database/3326:1)


## 4\) output format

网页格式的输出信息，使用说明见 [http://software.broadinstitute.org/gsea/doc/GSEAUserGuideTEXT.htm\#\_Interpreting\_GSEA\_Results](http://software.broadinstitute.org/gsea/doc/GSEAUserGuideTEXT.htm#_Interpreting_GSEA_Results)。


## 5\) 更多教程

### \(1\) GSEA for a user-defined ranked
在上述过程中，我们利用GSEA计算EIF4G2与表达矩阵中其他基因的相关性，并进行排序生成rank list。在实际应用中，用户自定义一个rank list也是常有方法。主要使用**Run GSEAPreranked**模块。

#### \(1a\) 文件准备

我们只需要两个输入：

1. 所有gene的rank list，文件后缀为.rnk。内容如下：
```text
# Rank	CRC-2211756-N
ENSG00000177106.14|5693	1062.683366
ENSG00000005206.16|7130	469.7135917
ENSG00000169994.18|8800	390.5745216
ENSG00000185189.17|4924	368.5153284
ENSG00000167702.11|4447	359.6025582
ENSG00000203697.11|5192	295.0883717
ENSG00000148600.14|7633	285.2801585
ENSG00000132470.13|6761	283.3625902
ENSG00000170421.12|6881	265.3386578
ENSG00000089356.18|3334	265.1819968
ENSG00000164877.18|8843	258.0913928
ENSG00000074964.16|7817	258.0291112
ENSG00000130827.6|11710	233.2619224
ENSG00000073350.13|7909	224.5944785
ENSG00000162337.11|6164	221.4648265
ENSG00000266714.7|12761	218.084697
ENSG00000130702.15|14291	215.763963
ENSG00000225697.12|6373	213.5188218
ENSG00000166165.12|2669	209.8102726
ENSG00000143850.13|8251	205.5522214
ENSG00000136059.14|6476	187.6522776
ENSG00000119514.6|3137	177.6641996
ENSG00000162804.13|13592	173.2296691
ENSG00000008710.19|18532	169.6201757
ENSG00000100376.11|6512	162.4053748
ENSG00000095066.11|5950	155.1013834
ENSG00000169894.17|11374	153.0630688
ENSG00000142173.14|5386	150.1852479
ENSG00000139631.18|4390	146.7598616
ENSG00000142798.17|17684	143.2301143
ENSG00000108679.12|3584	141.7435469
ENSG00000179832.17|10903	141.7308791
ENSG00000010327.10|10588	138.0681808
ENSG00000155256.17|3899	137.8860788
ENSG00000145214.13|4970	132.9209602
ENSG00000051009.10|4736	131.0598402
ENSG00000176454.13|4799	130.7545967
ENSG00000161714.11|6740	120.79598
ENSG00000160191.17|5266	120.4307842
ENSG00000107404.19|3814	120.0137943
ENSG00000149418.10|4350	113.813452
ENSG00000103168.16|5577	113.3672329
ENSG00000186866.16|5862	112.4799634
ENSG00000099139.13|12705	111.7105489
ENSG00000241769.7|3604	110.8399847
ENSG00000127415.12|4000	110.1100513
ENSG00000125912.10|4175	105.8455017
ENSG00000096433.10|9870	105.1426635
ENSG00000076356.6|13781	99.74836321
ENSG00000105429.12|11372	98.36870785
ENSG00000107281.9|2737	96.31407399
ENSG00000019144.18|21722	96.1187369
ENSG00000014914.20|4562	95.99596141
ENSG00000174233.11|8410	95.28798157
ENSG00000185000.11|4131	94.92802401
ENSG00000172830.12|3641	94.56875219
ENSG00000006025.11|6572	93.58012722
ENSG00000171914.16|16209	92.63865604
ENSG00000086015.20|7877	92.34223213
ENSG00000111057.10|3218	91.21401706
ENSG00000213213.13|2964	90.85212777
ENSG00000250251.6|6015	90.83299637
ENSG00000088256.8|6242	89.9146758
ENSG00000205746.9|6804	86.88438031
ENSG00000235333.3|978	85.92648969
ENSG00000188735.12|8449	85.31206271
ENSG00000106397.11|5671	83.10072741
ENSG00000073605.18|5905	82.76440518
ENSG00000142949.16|10712	82.37841817
ENSG00000244257.5|6805	80.98201936
ENSG00000158158.11|5105	77.16146546
ENSG00000215908.10|5606	77.12748751
ENSG00000182087.12|4031	77.02419836
ENSG00000167280.16|6173	75.68284094
ENSG00000141736.13|10321	74.14570915
ENSG00000130158.13|12435	73.56693491
ENSG00000132361.16|6900	73.30071797
ENSG00000223705.9|4822	72.81680404
ENSG00000107331.16|11514	72.6896517
ENSG00000080947.14|6813	72.33305408
ENSG00000162341.16|10283	70.76091616
ENSG00000039068.18|6177	70.69215624
ENSG00000072310.16|8030	70.6023166
ENSG00000119042.16|10958	70.46374227
ENSG00000254681.6|10291	69.39335286
ENSG00000180900.18|6192	69.35600504
ENSG00000114841.17|13872	69.11274773
ENSG00000124574.14|8083	69.08881855
ENSG00000183458.13|6830	68.38163901
ENSG00000142910.15|5281	68.2872206
ENSG00000140474.12|3516	66.63447917
ENSG00000120756.12|5294	66.43898419
ENSG00000151208.16|10271	66.20665628
ENSG00000054793.13|8317	65.64611499
ENSG00000125775.14|2653	65.63825222
ENSG00000171105.13|10451	65.05852321
ENSG00000176248.8|5447	64.8204574
ENSG00000214021.15|10365	64.61527679
ENSG00000178809.11|5544	64.42569891
```
2. 用户定义的gene set文件(或者GSEA自带的gene set），后缀为.gmx。在参数**Gene sets database**中选择自己的GMX文件。

>*  [GMX文件](http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMX:_Gene_MatriX_file_format_.28.2A.gmx.29)。GMX文件格式是制表符分隔的文件格式，用于描述基因集。在GMX格式中，每一列代表一个基因集。

我们给出一个GMX文件例子，文件存为后缀.gmt。
```text
Tumor_vs_Normal_tissue_up_regulated_logFC_top20	Tumor_vs_Normal_tissue_down_regulated_logFC_top20
Normal	Tissue
ENSG00000278505.4|1919	G032605|628
ENSG00000169994.18|8800	ENSG00000088256.8|6242
ENSG00000148702.14|3657	ENSG00000277739.1|153
ENSG00000134827.7|1586	G061517|25253
ENSG00000148600.14|7633	ENSG00000235333.3|978
ENSG00000083782.7|1805	ENSG00000214021.15|10365
ENSG00000089356.18|3334	G002038|10030
G066936|4422	G040639|13337
ENSG00000251026.1|574	ENSG00000125851.9|7078
ENSG00000130827.6|11710	ENSG00000175785.12|3892
G028862|10826	ENSG00000123560.13|6078
ENSG00000185479.5|2282	ENSG00000142959.4|2047
ENSG00000137745.11|2834	ENSG00000183034.12|2778
ENSG00000266714.7|12761	ENSG00000108231.12|11108
ENSG00000123500.9|3816	ENSG00000054793.13|8317
ENSG00000166165.12|2669	ENSG00000111404.6|1549
ENSG00000143850.13|8251	ENSG00000180900.18|6192
ENSG00000136059.14|6476	ENSG00000136546.13|8196
ENSG00000119514.6|3137	ENSG00000107331.16|11514
ENSG00000257046.5|3774	ENSG00000215908.10|5606
```
#### \(1b\) 软件运行
1）提前准备好软件和数据，根据上述过程将RNK文件(rank list文件)，GMX文件(gene set文件)，分别存为test.rnk和test.gmt，存在电脑本地。

2）打开 GSEA，点击**Load data**，选择文件夹下相应文件test.rnk和test.gmx。
![0a18773a71d6f9b762fecd30baed5313.gif](en-resource://database/3434:0)


3）点击**Run GSEAPreranked**。首先选择合适的参数，如下图所示。

* **Gene set database**选择上传的GMX文件，test.gmt。
* **Number of permutations**选择**1000**。置换检验的次数，数字越大结果越准确，但是太大会占用太多内存，软件默认检验1000次。
* **Ranked List**选择上传的RNK文件，test.rnk。
* **Collapse/Remap to gene symbols**选择**No_Collapse**。因为这里不需要进行转换，如果GMX文件和RNK文件基因格式不同，咋需要上传.chip进行转换，如上述教程。
* **Basic fields**中**Max size**需大于GMX文件中gene set大小，**Min size**需小于GMX文件gene set大小。
* 其他为默认参数。

![65f6f1d5bf06f08567fd54a00675c4bd.gif](en-resource://database/3432:0)

4）运行成功后，`Status` 会显示为 `Success`，点击其即可查看输出。
![d2a972a2c786c303c3bc4fc440a3d401.gif](en-resource://database/3436:0)


### \(2\) 分析离散型表型

GSEA也可以分析与离散型表型显著相关的 gene set （与差异表达相似），与上文EIF4G2基因分析教程不同的地方在于 `.cls` 文件。

以 [GSE2251](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2251) 例，打开链接后，可以在 Samples 一栏找到样品信息（样本数过多时，网页会显示不全，这时可以解压并打开 `_series_matrix.txt.gz` 查看样本信息。）。
本例共有12个样本，其中第1, 2, 3, 7, 8, 9 号为"E2-8"，第4, 5, 6, 10, 11, 12 号为"vehicle"。针对这一表型的`.cls` 文件如下所示。

```text
12 2 1
# E2-8 vehicle
0 0 0 1 1 1 0 0 0 1 1 1
```
`12` 是样本总数，`2` 是表型类型数目，`0` 代表第一种表型，`1` 代表第二种表型。

**注意：由于表型文件(cls文件)为离散型，选择参数时Basci fields中Metric for ranking genes应选择signal2noise、t-Test、ratio_of_class、 diff_of_class、log2_ratio_of_class之一。**


### \(3\) 常见错误：筛选后无基因集

由于基因芯片仅分析一部分基因，有些 gene set 会因为包含的基因过少而被剔除 。在极端情况下，可能所有的 gene set 都会被剔除，这时就会发生如下错误：

![](../../.gitbook/assets/gsea-no-gene-set-error.png)
![6d47b4a59cf99e130a08075612d07107.png](en-resource://database/3328:1)

更多信息请参见[http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/1001](http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/1001)



## 5\) Homework

1. 分析 GSE2251 中与 _TP53_ 基因的表达量显著相关的 gene set，输入文件在 [这里](https://cloud.tsinghua.edu.cn/d/747db0edd36449289b6f/?p=%2FFiles%2FPART_II%2F3.2.gsea&mode=list) 下载。

> 注意：要求提交中间过程和最后结果，截图并解释。

