### Structure-seq analysis of Ath (Shape-map)

以路径/Share2/home/lulab/liuxiaofan/shape_map/1.shapemapper/result_two/col_UV-_all_2 为例说明
#### 1.shapemapper
ShapeMapper automates the calculation of RNA structure probing reactivities from mutational profiling (MaP) experiments, in which chemical adducts on RNA are detected as internal mutations in cDNA through reverse transcription and read out by massively parallel sequencing.
The mutation rate (Mutr) at a given nucleotide is simply the mutation count (mismatches and unambiguously aligned deletions) divided by the read count at that location. 
Raw reactivities were generated for each nucleotide using thefollowing expression.
R = Mutrs-Mutru
S corresponds to a SHAPE modified sample, U to untreated.

由于内存问题，将不同染色体进行切分，详细文件见文件夹log

**chr1_.sh**
```
#!/bin/bash#BSUB -J col_UV-_1#BSUB -q Z-LU#BSUB -o output.%J#BSUB -e error.%J#BSUB -n 4#BSUB -R "span[hosts=1]"source ~/.bashrcset -e # exit on first error (if any)# Find the parent folder of this script,# resolving (possibly nested) symlinkscd /Share2/home/lulab/liuxiaofan/shape_map/1.shapemapper/result_two/col_UV-_all_2/BioII/lulab_b/shibinbin/projects/shapemap/shapemapper \--target /Share2/home/lulab/liuxiaofan/shape_map/3.reference_genome/Ath/my_tair10/Ensembl.tair10.fa.chr1.fa \--name "Ath_col0_UV-_all_chr1" \--min-depth 100 \--min-qual-to-count 20 \--overwrite \--modified --folder /Share2/home/lulab/liuxiaofan/shape_map/0.raw/modified_col0_UV-/two  \--untreated --folder /Share2/home/lulab/liuxiaofan/shape_map/0.raw/control_col0_UV-/two \--nproc 1 \--verboserm -rf /Share2/home/lulab/liuxiaofan/shape_map/1.shapemapper/result_two/col_UV-_all_2/shapemapper_temp/Ath_col0_UV-_all_chr1


```
##### input
1.Ensembl.tair10.fa.2.fa: target fasta file (all chr)
2.Shape02: modified fasta file
3.Shape01: untreated fasta file
##### output
1.Ath_all_100_chr_histograms.pdf
2.Ath_all_100_chr_profiles.pdf 
3.Ath_all_100_chr.shape
column1: nucleotide number 
column2: reactivity
4.example-results_chr.map
column1: nucleotide number 
column2: reactivity
column3: standard error
column4: nucleotide sequence
5.example-results_chr_profile.txt

|Nucleotide| Sequence| 	Modified_mutations|	Modified_read_depth	|Modified_effective_depth|Modified_rate|Untreated_mutations |Untreated_read_depth|	Untreated_effective_depth|	Untreated_rate|	Denatured_mutations|Denatured_read_depth|Denatured_effective_depth|	Denatured_rate|	Reactivity_profile|	Std_err|HQ_profile|HQ_stderr|Norm_profile|Norm_stderr
| --- | --- | --- | --- | --- | --- | --- | --- | ---  | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | 
1|	g|	0	|4788|	4335|	0.000000|	0|	5206|	4666|	0.000000|	0|	nan	|0.000000	|0.000000|	nan|	nan|	nan|	nan|
2|	g|0|	4837|	3405	|0.000000|	0|	5270	|3643|	0.000000|	0|	nan	|0.000000	|0.000000	|nan|	nan	|nan|	nan|

chr: chr name in target fasta file, like 1,2.

#### 2.Analyze only the exon region 
**readme_depth.sh :**
**input**
path0: code path
path1: output path
path2: exons.gtf path
path3: input path
**output**
1.final.modified : Each row is a gene, "c(0,1,0,...)". Each  point is Modified_mutations of each Nucleotide of this gene.
2.final.modified_depth : Each row is a gene, "c(0,1,0,...)". Each  point is Modified_effective_depth of each Nucleotide of this gene.
3.final.unmodified : Each row is a gene, "c(0,1,0,...)". Each  point is Untreated_mutations of each Nucleotide of this gene.
4.final.unmodified_depth : Each row is a gene, "c(0,1,0,...)". Each  point is Untreated_effective_depth of each Nucleotide of this gene.
```
#!/bin/bash#BSUB -J col_uv-  # job name#BSUB -q Z-LU      # queue#BSUB -o output.%J      # output file name in which %J is replaced by the job ID#BSUB -e error.%J       # error file name in which %J is replaced by the job ID#BSUB -n 1      # number of threads/processes in a job#BSUB -R "span[hosts=1]"        # do not run single jobs on multiple nodespath0=/Share2/home/lulab/liuxiaofan/shape_map/bin2;path1=/Share2/home/lulab/liuxiaofan/shape_map/1.shapemapper/result_two/col_UV-_all_2/shapemapper_out/trans_depth2;path3=/Share2/home/lulab/liuxiaofan/shape_map/1.shapemapper/result_two/col_UV-_all_2/shapemapper_out;mkdir -p $path1;path2=/Share2/home/lulab/liuxiaofan/shape_map/3.reference_genome/Ath/my_tair10/fa;input_gtf=/Share2/home/lulab/liuxiaofan/shape_map/3.reference_genome/Ath/my_tair10/Output/exons.gtf;source /BioII/lulab_b/containers/singularity/wrappers/bashrcrm $path1/foorm $path1/exons*awk -F '[\t;]' '{print $10"\t"$1"."$7"\t"$4"\t"$5}'             $input_gtf              |sed 's/ transcript_id //g;s/"//g'|sort -k1,1 -k2,2 
-k3,3n -k4,4n       >       $path1/foo;awk -F '\t' 'BEGIN{ID="";chr="";coor=""}{if(ID!=$1){print ID"\t"chr"\tc("coor")";ID=$1;chr=$2;coor=$3":"$4;}else{coor=coor","$3":"$4;}}END{print ID"\t"chr"\tc("coor")";}' $path1/foo             |sed '1d'                                                               
>       $path1/exons.coor;awk -F '\t' -v prefix="$path1/exons.coor"  '{print $0 >> prefix"."$2}' $path1/exons.coor;chmod 764 $path1/exons.coor;chmod 764 $path1/exons.coor.*;### transcriptomefor chr in 1 2 3 4 5 6 7;do echo $chr;if [ "$chr" == "6" ]; then Chr="Mt";elif [ "$chr" == "7" ]; then Chr="Pt";elseChr=$chr;fi;echo $chr"\t"$Chr;awk -F '\t' 'NR>1{print $1"\t"$2}'                $path3/"Ath_col0_UV-_all_chr"$chr"_"$chr"_profile.txt"         >       
$path1/Ath.$chr.Nucleotide.txt;awk -F '\t' 'NR>1{print $1"\t"$3}'                $path3/"Ath_col0_UV-_all_chr"$chr"_"$chr"_profile.txt"         >       
$path1/Ath.$chr.modified.txt;awk -F '\t' 'NR>1{print $1"\t"$5}'                $path3/"Ath_col0_UV-_all_chr"$chr"_"$chr"_profile.txt"        >       
$path1/Ath.$chr.modified_depth.txt;awk -F '\t' 'NR>1{print $1"\t"$7}'                $path3/"Ath_col0_UV-_all_chr"$chr"_"$chr"_profile.txt"        >       
$path1/Ath.$chr.unmodified.txt;awk -F '\t' 'NR>1{print $1"\t"$9}'                $path3/"Ath_col0_UV-_all_chr"$chr"_"$chr"_profile.txt"        >       
$path1/Ath.$chr.unmodified_depth.txt;geno=`awk -F '\t' -v hah=$chr 'NR==hah{print $2}' $path2/tair10.genome`;echo $geno;chmod 764 $path1/Ath.*;for state in modified unmodified;do echo $state;Rscript $path0/generate_HDF5.R          $path1/Ath.$chr.$state.txt      $geno   $path1/Ath.$chr.$state.HDF5;chmod 764 $path1/Ath*;Rscript $path0/generate_HDF5.R          $path1/Ath.$chr.$state"_depth".txt      $geno   $path1/Ath.$chr.$state"_depth".HDF5;chmod 764 $path1/Ath*;for str in + -;do echo $str;Rscript $path0/HDF5_TR_psite_track.R           $path1/Ath.$chr.$state.HDF5    $path1/exons.coor.$Chr.$str    $str    
$path1/exons.$state.$chr.$str;chmod 764 $path1/exons*;Rscript $path0/HDF5_TR_psite_track.R           $path1/Ath.$chr.$state"_depth".HDF5    $path1/exons.coor.$Chr.$str    $str    
$path1/exons.$state"_depth".$chr.$str;chmod 764 $path1/exons*;done;done;for state in Nucleotide;do echo $state;Rscript $path0/generate_HDF5.2.R          $path1/Ath.$chr.$state.txt      $geno   $path1/Ath.$chr.$state.HDF5;chmod 764 $path1/Ath*;for str in + -;do echo $str;Rscript $path0/HDF5_TR_psite_track.R           $path1/Ath.$chr.$state.HDF5    $path1/exons.coor.$Chr.$str    $str    
$path1/exons.$state.$chr.$str;chmod 764 $path1/exons*;done;done;done;### summarizefor chr in 1 2 3 4 5 6 7;do echo $chr;for str in + -;do echo $str;for state in modified unmodified;do echo $state;Rscript $path0/summarize.R              $path1/exons.$state.$chr.$str   $path1/exons.$state.$chr.$str.summary;chmod 764 $path1/exons*;Rscript $path0/summarize.R              $path1/exons.$state"_depth".$chr.$str   $path1/exons.$state"_depth".$chr.$str.summary;chmod 764 $path1/exons*;done;done;done;### combinerm -rf $path1/final.modified $path1/final.modified_depth $path1/final.unmodified $path1/final.unmodified_depth;for state in modified unmodified;do echo $state;for chr in 1 2 3 4 5 6 7;do echo $chr;for str in + -;do echo $str;cat $path1/exons.$state.$chr.$str >> $path1/final.$state;wc -l $path1/exons.$state.$chr.$str;chmod 764 $path1/final*;cat $path1/exons.$state"_depth".$chr.$str >> $path1/final.$state"_depth";wc -l $path1/exons.$state"_depth".$chr.$str;chmod 764 $path1/final*;done;done;done;### summarizefor chr in 1 2 3 4 5 6 7;do echo $chr;for str in + -;do echo $str;for state in Nucleotide;do echo $state;Rscript $path0/summarize2.R              $path1/exons.$state.$chr.$str   $path1/exons.$state.$chr.$str.summary;chmod 764 $path1/exons*;done;done;done;### combinefor state in Nucleotide;do echo $state;for chr in 1 2 3 4 5 6 7;do echo $chr;for str in + -;do echo $str;cat $path1/exons.$state.$chr.$str >> $path1/final.$state;wc -l $path1/exons.$state.$chr.$str;chmod 764 $path1/final*;done;done;done;


```

#### 3.hit level
The hit level metric quantifies the total background-subtracted signal per nucleotide of transcript.
![42dfca9c3692367939b123d0f1a2c617.png](en-resource://database/891:1)
where the subscripts S and B indicate the experimental sample and background control, respectively; events are either ligation-detected sequence stops or mutations, depending on readout method, and read depth corresponds to the median number of reads overlapping each nucleotide in the transcript. 

**readme_hit.sh**

```
#!/bin/shpath0=/Share2/home/lulab/liuxiaofan/shape_map/bin2;path1=/Share2/home/lulab/liuxiaofan/shape_map/1.shapemapper/result_two/col_UV-_all_2/shapemapper_out/trans_hit2;path2=/Share2/home/lulab/liuxiaofan/shape_map/1.shapemapper/result_two/col_UV-_all_2/shapemapper_out/trans_depth2;mkdir -p $path1;source /BioII/lulab_b/containers/singularity/wrappers/bashrc### hit level### calculate the average mutation events above background for each transcriptmodified=$path2/final.modified; #modified eventsunmodified=$path2/final.unmodified; # modified depthmodified_depth=$path2/final.modified_depth; #unmodified eventsunmodified_depth=$path2/final.unmodified_depth; # unmodified depthNucleotide=$path2/final.Nucleotide;### hit level > 5 can almost recover the RNA structure accuratelyrm -rf $path1/final.modified_unmodified.hit;python $path0/hit_level.py          --modified $modified \--unmodified $unmodified \--modified_depth $modified_depth \--unmodified_depth $unmodified_depth \--Nucleotide $Nucleotide \--savepath  $path1/final.modified_unmodified \--savepath_hit  $path1/final.modified_unmodified.hit;chmod 764 $path1/final.modified_unmodified.hit*;echo -e "cutoff\tmodified.median\tunmodified.median\thit"   >           $path1/cutoff.hit.group;awk -F '\t' '$3>0 && $2<=25 && $2>0{print "0\t"$0}'         $path1/final.modified_unmodified.hit  >>        $path1/cutoff.hit.group;awk -F '\t' '$3>0 && $2<=50 && $2 > 25{print "25\t"$0}'     $path1/final.modified_unmodified.hit  >>        $path1/cutoff.hit.group;awk -F '\t' '$3>0 && $2<=100 && $2 > 50{print "50\t"$0}' $path1/final.modified_unmodified.hit >>      $path1/cutoff.hit.group;awk -F '\t' '$3>0 && $2<=200 && $2 > 100{print "100\t"$0}' $path1/final.modified_unmodified.hit >>      $path1/cutoff.hit.group;awk -F '\t' '$3>0 && $2<=300 && $2 > 200{print "200\t"$0}' $path1/final.modified_unmodified.hit >>      $path1/cutoff.hit.group;awk -F '\t' '$3>0 && $2<=500 && $2 >300{print "300\t"$0}' $path1/final.modified_unmodified.hit >>      $path1/cutoff.hit.group;awk -F '\t' '$3>0 && $2<=750 && $2 >500{print "500\t"$0}' $path1/final.modified_unmodified.hit >>      $path1/cutoff.hit.group;awk -F '\t' '$3>0 && $2<=1000 && $2 >750{print "750\t"$0}' $path1/final.modified_unmodified.hit >>      $path1/cutoff.hit.group;awk -F '\t' '$3>0 && $2<=2000 && $2 >1000{print "1000\t"$0}' $path1/final.modified_unmodified.hit >>      $path1/cutoff.hit.group;awk -F '\t' '$3>0 && $2<=5000 && $2 >2000{print "2000\t"$0}' $path1/final.modified_unmodified.hit >>      $path1/cutoff.hit.group;awk -F '\t' '$3>0 && $2 > 5000{print "5000\t"$0}' $path1/final.modified_unmodified.hit >>      $path1/cutoff.hit.group;chmod 764 $path1/cutoff.hit.group*;

```

#### 4.数据归一化及基因改变区域
python代码在路径/Share2/home/lulab/liuxiaofan/shape_map/4.code

1.hit_level.py 画出hit level图
2.data_preprocessing.py 数据归一化
3.gini_summary.py 每个基因窗口下gini index的差值，最大值，是否存在结构改变区域4.gini_transcript.py 每个存在结构改变区域的具体信息
