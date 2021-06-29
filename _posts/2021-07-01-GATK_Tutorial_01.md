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

## 1 Fastq文件拼装

### 1.1 原始fastq质控

```bash
fastp -i Input_R1.fastq.gz\
      -I Input_R2.fastq.gz\
      -o Output_R1.fastp.gz\
      -O Output_R2.fastp.gz\
      --jason
      --html
```

软件：fastp

下载地址：

常用指令：

- i/I
- o/O
- 

### 1.2 fastq组装

```bash
bwa mem | \
samtools view |\
samtools
```

软件：bwa

下载地址：

软件：samtools

下载地址：

### 1.3 BQSR

### 1.4 SNP获取

- 肿瘤/正常配对数据

  ```
  ```

  

- 单一体细胞数据

  ```
  
  
  ```

### 1.5 数据筛选

- 肿瘤、正常配对数据
- 单一体细胞数据

### 1.6 基因注释

