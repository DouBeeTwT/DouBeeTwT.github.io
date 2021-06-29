---
layout: post
title: 'GATK4基础流程SNP筛选'
subtitle: ''
date: 2021-07-01
categories:
 - tutorial
 - Updating
tags: GATK
---

## 0 软件准备

| 软件名称 | 下载地址                                                     |
| -------: | :----------------------------------------------------------- |
|    fastp | https://github.com/OpenGene/fastp/releases                   |
|      bwa | https://github.com/lh3/bwa/releases                          |
| samtools | https://github.com/samtools/samtools/releases                |
|     gatk | https://github.com/broadinstitute/gatk/releases              |
|  annovar | https://annovar.openbioinformatics.org/en/latest/user-guide/download/ |

\* 软件版本可能随时间改变而改变，连接可能会失效

## 1 Fastq文件拼装

### 1.1 原始fastq质控

```bash
fastp -i [R1_In] -I [R2_In] \
      -o [R1_Out] -O [R2_Out] \
      --json [Name].json \
      --html [Name]].html
```



### 1.2 fastq组装

```bash
Thread_Num=[number]
bwa mem -t $Thread_Num -R  "[Head Infomation]" [Genome] [R1] [R2] | \
samtools view -b -S -@ $Thread_Num | \
samtools sort -@ $Thread_Num -o [Name].sorted.bam
```



### 1.3 BQSR

```bash
GATK MarkDuplicates -I [Name].sorted.bam -O [Name].marked.bam -M [Name].metrics

GATK FixMateInformation -I [Name].marked.bam -O [Name].fixed.bam -SO coordinate

GATK BaseRecalibrator \
     -R [Genome] -I [Name].fixed.bam \
     --known-sites [db_snp] \
     --known-sites [dn_indel] \
     -O [Name].recal.table
 
GATK ApplyBQSR -R [Genome] -I [Name].fixed.bam -bqsr [Name].recal.table -O [Name].bqsr.bam
```



### 1.4 SNP获取

- 肿瘤/正常配对数据（Mutect2）

  ```bash
  GATK Mutect2 \
       -R [Genome] -tumor [Tumor] -I [Tumor].bqsr.bam \
       --germline-resource [AF_ONLY_GENOMAD] \
       -O [Name].vcf
  
  GATK Mutect2 \
       -R [Genome] -normal [Normal] -tumor [Tumor] \
       -I [Normal].bqsr.bam -I [Tumor].bqsr.bam \
       --panel-of-normals [Name].vcf \
       --germline-resource [AF_ONLY_GENOMAD] \
       --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
       -O [Name].vcf
  ```

  

- 单一体细胞数据（HaplotypeCaller）

  ```bash
  GATK HaplotypeCaller -R [Genome] -I [Name].bqsr.bam --dbsnp [dbsnp] -O [Name].vcf
  ```

### 1.5 数据筛选

- 肿瘤、正常配对数据

  ```bash
  
  ```

  

- 单一体细胞数据

  ```bash
  ### ========== Filter of Snp ========== ###
  
  # Select
  GATK SelectVariants -R [Genome] -V [In].vcf -select-type SNP -O [Out].vcf
  
  # VQSR
  GATK VariantRecalibrator \
       -R [Genome] -V [Name].vcf -mode SNP \
       --resource:hapmap,known=false,training=true,truth=true,prior=15.0 [Hapmap] \
       --resource:omini,known=false,training=true,truth=false,prior=12.0 [Omini] \
       --resource:1000G,known=false,training=true,truth=false,prior=10.0 [1000G_snp] \
       --resource:dbsnp,known=true,training=false,truth=false,prior=6.0  [dbsnp] \
        -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \
        --tranches-file [Name].tranches -O [Name].recal
  
  GATK ApplyVQSR \
       -R [Genome] -V [Name].vcf -mode SNP \
       --truth-sensitivity-filter-level 99.5 \
       --tranches-file [Name].tranches \
       --recal-file [Name].recal \
       -O [Name].VQSR.vcf
  
  # Filter
  GATK VariantFiltration \
       -R [Genome] -V [Name].VQSR.vcf \
       --filter-expression "DP<30.0||QD<2.0||FS>60.0||MQ<40.0||MQRankSum< -12.5||ReadPosRankSum< -8.0" \
       --filter-name "My_Snp_Filter" \
       -O [Name].filtered.vcf
   
  GATK SelectVariants -R [Genome] --exclude-filtered -V [Name].filtered.vcf -O [Name].PASS.vcf
  
  ### ========== Filter of Ind ========== ###
   
  # Select
  GATK SelectVariants -R [Genome] -select-type INDEL -V [Name].vcf -O [Name].vcf
   
  # VQSR
  GATK VariantRecalibrator \
       -R [Genome] -V [Name].vcf -mode INDEL \
       --resource:mills,known=true,training=true,truth=true,prior=12.0 [db_indel] \
       -an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
       --tranches-file [Name].tranches -O [Name].recal
   
  GATK ApplyVQSR \
       -R [Genome] -V [Name].vcf -mode INDEL \
       --truth-sensitivity-filter-level 99.0 \
       --tranches-file [Name].tranches \
       --recal-file [Name].recal \
       -O [Name].VQSR.vcf \
   
  # Filter
  GATK VariantFiltration \
       -R [Genome] -V [Name].VQSR.vcf \
       --filter-expression "DP < 30.0 || QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"  \
       --filter-name "My_Indel_Filter" \
       -O [Name].filtered.vcf
   
  GATK SelectVariants -R [Genome] --exclude-filtered -V [Name].filtered.vcf -O [Name].PASS.vcf
  ```

  

### 1.6 基因注释

```bash
table_annovar.pl [Name].vcf [ANNOVA Home Pathway]/humandb/ \
            -buildver [version] \
            -out [Name].annovar \
            -protocol [refGene,cytoBand,exac03,avsnp147,dbnsfp30a] \
            -operation [gx,r,f,f,f] \
            -nastring . -vcfinput --otherinfo
```

