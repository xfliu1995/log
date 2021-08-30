# 1.1.Biological information Tools

本节教学介绍两个实用的生信分析工具，bedtools和Samtools。


## 1\) Files Needed

### 1a\) Method 1: Use docker

docker images的下载链接如[附表](https://lulab2.gitbook.io/teaching/appendix/appendix-iv.-teaching)所示，加载完我们提供的image后，文件都已经准备好了，可以这样查看：

```
docker load -i bioinfo_tools.tar.gz

docker run -dt --name bioinfo_tools --restart unless-stopped -v ~/Downloads/data:/data xfliu1995/bioinfo_tools:1.0

cd /home/test/bioinfo_tools
```

> 本教程docker使用方式：
> * 1) 运行容器:  docker exec -it motif bash
> * 2) 进行Linux系统的相关操作
> * 3) 退出容器：exit

### 1b\) Method 2: Directly Download Files Needed

* 如果不使用docker，也可以直接下载教程所需文件：Download Link

### 1c\) Method 3: Cluster 
如果使用P集群完成该节课程，不需要singularity，在集群中提供所需文件及软件绝对路径。

* 文件位置
```bash
/data/images/Biological_data
```

* 软件地址
```bash
/data/liuxiaofan/software/anaconda3/bin/samtools
/data/liuxiaofan/software/anaconda3/envs/r403/bin/bedtools
```

## 2\) bedtools

Bedtools是由犹他⼤学昆兰实验室开发的基因组算法⼯具集，它堪称是基因组分析⼯具中的瑞⼠军
⼑。Bedtools可以对基因组⼴泛使⽤的数据格式BAM,BED,GFF/GTF,VCF进⾏处理，进⾏取交集、并
集、补集、计数以及格式转变等操作。具体可见官方[教程](https://bedtools.readthedocs.io/en/latest/index.html)。

较典型而且常用的功能如下

| 函数 |功能  |
| --- | --- |
| bamToBed | 格式转换，bam转bed |
| bedToBam，bedToIgv | bed转其他格式 |
| intersectBed，windowBed |对基因组坐标的逻辑运算，交集  |
| closestBed | 对基因组坐标的逻辑运算，邻集 |
| complementBed | 对基因组坐标的逻辑运算，补集 |
| mergeBed | 对基因组坐标的逻辑运算，并集 |
| subtractBed | 对基因组坐标的逻辑运算，差集 |
| coverage，genomecov | 对基因组坐标的逻辑运算，计算覆盖度 |

本章中我们主要介绍bedtools常用的三个函数，intersect、merge、multicov。

### 2.1\) bedtools intersect

它的作用是比较两个或多个BED/BAM/VCF/GFF，然后找它们重叠的区域（就是指至少存在1个公共的碱基）
![8dd32ec1878510d60addbc83ccb65439.png](evernotecid://58A46CB2-04A5-4F32-9233-277C274EDE1B/appyinxiangcom/19737026/ENResource/p475)
![f54dd65f091fcb96d07bd02e1ef0184f.png](evernotecid://58A46CB2-04A5-4F32-9233-277C274EDE1B/appyinxiangcom/19737026/ENResource/p476)

* 没有参数：简单的返回A和B共有的区域
* -wa参数：把原始A文件的坐标输出来
* -wa、-wb参数：先输出A的坐标，然后输出与A有重叠的B1、B2的区域坐标
* -v参数：只输出A中有且B中没有的区域

#### (1) 默认的输出结果
```
bedtools intersect -a cpg.bed -b exons.bed | head -5

chr1	29320	29370	CpG:_116
chr1	135124	135563	CpG:_30
chr1	327790	328229	CpG:_29
chr1	327790	328229	CpG:_29
chr1	327790	328229	CpG:_29
```
结果并不是原始的CpG记录的区间，而是它和exon.bed的重叠区间，因此坐标只是原来CpG.bed的子集。

#### (2) 重叠双方的具体区域

```
bedtools intersect -a cpg.bed -b exons.bed -wa -wb | head -5

chr1	28735	29810	CpG:_116	chr1	29320	29370	NR_024540_exon_10_0_chr1_29321_r	0	-
chr1	135124	135563	CpG:_30	chr1	134772	139696	NR_039983_exon_0_0_chr1_134773_r	0	-
chr1	327790	328229	CpG:_29	chr1	324438	328581	NR_028322_exon_2_0_chr1_324439_f	0	+
chr1	327790	328229	CpG:_29	chr1	324438	328581	NR_028325_exon_2_0_chr1_324439_f	0	+
chr1	327790	328229	CpG:_29	chr1	327035	328581	NR_028327_exon_3_0_chr1_327036_f	0	+
```

重叠区域CpG.bed坐标和外显子坐标，例如第一行28735-29810是CpG.bed坐标，29320-29370是外显子坐标。

#### (3) 统计重叠的碱基数

使用参数-wo （write overlap），会在每一行结果末尾统计这个重叠区域有多长。

```
bedtools intersect -a cpg.bed -b exons.bed -wo | head -5

chr1	28735	29810	CpG:_116	chr1	29320	29370 NR_024540_exon_10_0_chr1_29321_r	0	-	50
chr1	135124	135563	CpG:_30	chr1	134772	139696	NR_039983_exon_0_0_chr1_134773_r	0	-	439
chr1	327790	328229	CpG:_29	chr1	324438	328581	NR_028322_exon_2_0_chr1_324439_f	0	+	439
chr1	327790	328229	CpG:_29	chr1	324438	328581	NR_028325_exon_2_0_chr1_324439_f	0	+	439
chr1	327790	328229	CpG:_29	chr1	327035	328581	NR_028327_exon_3_0_chr1_327036_f	0	+	439
```
例如第一行CpG.bed的28735-29810区域碱基数为50

#### (4) 统计发生重叠的features数量
使用参数-c ，就会统计对文件A中的每个feature，在文件B中有多少个features与之有交集。这个feature的含义很广，不限于基因，还可以是外显子、转录本等等与基因组有关的信息。
```
bedtools intersect -a cpg.bed -b exons.bed -c | head -5

chr1	28735	29810	CpG:_116	1
chr1	135124	135563	CpG:_30	1
chr1	327790	328229	CpG:_29	3
chr1	437151	438164	CpG:_84	0
chr1	449273	450544	CpG:_99	0
```

#### (5)统计没有重叠的区域
使用-v参数
```
bedtools intersect -a cpg.bed -b exons.bed -v | head -5

chr1	437151	438164	CpG:_84
chr1	449273	450544	CpG:_99
chr1	533219	534114	CpG:_94
chr1	544738	546649	CpG:_171
chr1	801975	802338	CpG:_24
```

#### (6)设定一个最小重叠比例
软件默认只要二者有1bp交叉，也算重叠区域。但这个我们可以自行调整。例如，设定至少满足文件A的百分之50要与文件B有重叠，才统计。使用参数-f （fraction）。
```
bedtools intersect -a cpg.bed -b exons.bed -wo -f 0.5 | head -5

chr1	135124	135563	CpG:_30	chr1	134772	139696	NR_039983_exon_0_0_chr1_134773_r	0	-	439
chr1	327790	328229	CpG:_29	chr1	324438	328581	NR_028322_exon_2_0_chr1_324439_f	0	+	439
chr1	327790	328229	CpG:_29	chr1	324438	328581	NR_028325_exon_2_0_chr1_324439_f	0	+	439
chr1	327790	328229	CpG:_29	chr1	327035	328581	NR_028327_exon_3_0_chr1_327036_f	0	+	439
chr1	788863	789211	CpG:_28	chr1	788770	794826	NR_047519_exon_5_0_chr1_788771_f	0	+	348
```

### 2.2\) bedtools merge
许多数据集的基因组feature坐标经常是连续的，就像下图的蓝色部分于是可以把这些连续的基因组小区间连接起来，拼成一个连续的大区间。注意bedtools merge需要拼接的输入文件（bed/gff/vcf，必须是排序（sort）过的。

![763551cc22f11487ec3bf137410bfca5.png](evernotecid://58A46CB2-04A5-4F32-9233-277C274EDE1B/appyinxiangcom/19737026/ENResource/p477)

以exon.bed为例，区间是有重叠的，因此它们可以进行merge。
```
head -n10 exons.bed

chr1    11873   12227   NR_046018_exon_0_0_chr1_11874_f 0       +
chr1    12612   12721   NR_046018_exon_1_0_chr1_12613_f 0       +
chr1    13220   14409   NR_046018_exon_2_0_chr1_13221_f 0       +
chr1    14361   14829   NR_024540_exon_0_0_chr1_14362_r 0       -
chr1    14969   15038   NR_024540_exon_1_0_chr1_14970_r 0       -
chr1    15795   15947   NR_024540_exon_2_0_chr1_15796_r 0       -
chr1    16606   16765   NR_024540_exon_3_0_chr1_16607_r 0       -
chr1    16857   17055   NR_024540_exon_4_0_chr1_16858_r 0       -
chr1    17232   17368   NR_024540_exon_5_0_chr1_17233_r 0       -
chr1    17605   17742   NR_024540_exon_6_0_chr1_17606_r 0       -
```

第3行和第4行的区间是有重叠的，因此它们可以进行merge。

#### (1)合并
```
bedtools merge -i exons.bed | head -10

chr1    11873   12227
chr1    12612   12721
chr1    13220   14829
chr1    14969   15038
chr1    15795   15947
chr1    16606   16765
chr1    16857   17055
chr1    17232   17368
chr1    17605   17742
chr1    17914   18061
```
可以看到，merge之后原来的第3行（13220 - 14409）和第4行（14361 - 14829）坐标合并成了13220 - 14829。

#### (2)merge后统计之前重叠区间的数量

我们希望能直观看出新的区间13220-14829 是不是合并后的；如果是，合并了多少重叠的区域。
```
# 仅仅统计次数就用count【可以指定第一列chr，但实际上还是按第2-3列的坐标进行识别】
bedtools merge -i exons.bed -c 1 -o count | head -n 10

chr1	11873	12227	1
chr1	12612	12721	1
chr1	13220	14829	2
chr1	14969	15038	1
chr1	15795	15947	1
chr1	16606	16765	1
chr1	16857	17055	1
chr1	17232	17368	1
chr1	17605	17742	1
chr1	17914	18061	1
```

第3行（13220 - 14409）和第4行（14361 - 14829）坐标合并成13220 - 14829，最后一列为区间个数，2。

#### (3)合并指定距离区间

使用-d参数（distance），认为相聚不超过distance的两个区间，依然表示同一个，那么也可以合并。以1000bp为例，如下。
```
bedtools merge -i exons.bed -d 1000 -c 1 -o count | head -10

chr1	11873	18366	12
chr1	24737	24891	1
chr1	29320	29370	1
chr1	34610	36081	6
chr1	69090	70008	1
chr1	134772	140566	3
chr1	323891	328581	10
chr1	367658	368597	3
chr1	621095	622034	3
chr1	661138	665731	3
```

#### (4)merge信息展示
我们可以对指定列的信息进行输出，使用-o collapse参数。-c参数是指对哪一列进行操作。如果存在合并的现象，那么其中所有的都会列出来，并以逗号分隔
```
bedtools merge -i exons.bed -d 90 -c 1,4 -o count,collapse | head -10

chr1	11873	12227	1	NR_046018_exon_0_0_chr1_11874_f
chr1	12612	12721	1	NR_046018_exon_1_0_chr1_12613_f
chr1	13220	14829	2	NR_046018_exon_2_0_chr1_13221_f,NR_024540_exon_0_0_chr1_14362_r
chr1	14969	15038	1	NR_024540_exon_1_0_chr1_14970_r
chr1	15795	15947	1	NR_024540_exon_2_0_chr1_15796_r
chr1	16606	16765	1	NR_024540_exon_3_0_chr1_16607_r
chr1	16857	17055	1	NR_024540_exon_4_0_chr1_16858_r
chr1	17232	17368	1	NR_024540_exon_5_0_chr1_17233_r
chr1	17605	17742	1	NR_024540_exon_6_0_chr1_17606_r
chr1	17914	18061	1	NR_024540_exon_7_0_chr1_17915_r
```
exon.bed中的第4列是exon的名称，这样在count结束后，就可以再根据第4列输出名称。

### 2.3\) bedtools multicov
bedtools multicov，报告来自多个位置排序和索引的 BAM 文件的对齐计数，这些 BAM 文件在 BED 文件中重叠间隔。具体来说，对于提供的每个 BED 间隔，它报告来自每个 BAM 文件的重叠对齐的单独计数。

```
multiBamCov [OPTIONS] -bams BAM1 BAM2 BAM3 ... BAMn -bed  <BED/GFF/VCF>

cat multicov.bed

chr1    15      20      a1      1       +
chr1    15      27      a2      2       +
chr1    15      20      a3      3       -
chr1    15      27      a4      4       -

#这里输入的bam文件必须是有index
samtools index one_block.bam
samtools index two_blocks.bam

bedtools multicov -bams one_block.bam  two_blocks.bam -bed multicov.bed

chr1    15      20      a1      1       +       1       1
chr1    15      27      a2      2       +       1       1
chr1    15      20      a3      3       -       1       1
chr1    15      27      a4      4       -       1       1
```
输出反映了文件中每个记录重叠对齐的独特报告。在上面的示例中，输出的每一行都反映了 a)文件中的原始行，然后是b)与每个输入文件的间隔重叠的对齐计数。在上面的例子中，输出由 7 列组成：前四列是来自文件的列，最后三列是来自 3 个输入 文件的重叠对齐的计数。


## 3\) samtools

samtools是一个用于操作sam和bam文件的工具合集。samtools最常用的功能是格式转换，排序，索引等等。
| 函数 |功能  |
| --- | --- |
| view | BAM-SAM/SAM-BAM 转换和提取部分比对 |
| sort | 比对排序 merge: 聚合多个排序比对 |
| index | 索引排序比对 |
| faidx | 建立FASTA索引，提取部分序列 |
| tview | 文本格式查看序列 |
| pileup | 产生基于位置的结果和 consensus/indel calling |
| rmdup | 过滤PCR重复 |

本章中我们介绍samtools常用函数view，sort，index，flagstat，

### 3.1 \) view
view命令的主要功能是：将sam文件转换成bam文件；然后对bam文件进行各种操作，比如数据提取(这些操作是对bam文件进行的，因而当输入为sam文件的时候，不能进行该操作)。

#### (1)将sam文件、bam文件互相转换
```
samtools view -h Shape01.bam > Shape01.sam
amtools view -b -S Shape01.sam > Shape01.bam

```

#### (2) 序列提取
```
#提取比对到参考序列上的比对结果
samtools view -bF 4 Shape01.bam > Shape01.F.bam
 
#提取paired reads中两条reads都比对到参考序列上的比对结果，只需要把两个4+8的值12作为过滤参数即可
samtools view -bF 12 Shape01.bam > Shape01.F12.bam
 
#提取没有比对到参考序列上的比对结果
samtools view -bf 4 Shape01.bam > Shape01.f.bam
 
#提取比对质量高的reads
samtools view -q 1 -O bam -o Shape01.highQual.bam Shape01.bam

```

### 3.2 \) sort
sort对bam文件进行排序。
```
#按默认的染色体位置进行排序,-@是线程数
samtools sort -@8 Shape01.bam Shape01.sort

-n参数则是根据read名进行排序,-t 根据TAG进行排序
samtools sort -n Shape01.bam Shape01.sort.read
```

### 3.3 \) index
必须对bam文件进行默认情况下的排序后，才能进行index。否则会报错。建立索引后将产生后缀为.bai的文件，用于快速的随机处理。很多情况下需要有bai文件的存在，特别是显示序列比对情况下。比如samtool的tview命令就需要；gbrowse2显示reads的比对图形的时候也需要。
```
samtools index Shape01.sort.bam
```

### 3.4 \) flagstat
给出BAM文件的比对结果。
```
samtools flagstat Shape01.bam

49301021 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
49301021 + 0 mapped (100.00% : N/A)
49301021 + 0 paired in sequencing
24650545 + 0 read1
24650476 + 0 read2
49300908 + 0 properly paired (100.00% : N/A)
49300908 + 0 with itself and mate mapped
113 + 0 singletons (0.00% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```

### 3.5 \) rmdup

过滤PCR重复,去除duplicate。

```
samtools rmdup Shape01.sorted.bam Shape01.rmdup.bam
```


## 4\) Homework

根据我们给出的exon.bed进行合并，对Shape01.bam进行排序并提取比对质量高的reads，然后计算BAM文件在的exon.bed上对齐计数，并输出前10行。
