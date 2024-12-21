# S_cerevisiae_DLA_Response
Transcriptomic response of Saccharomyces cerevisiae to D-lactic acid
# Transcriptional Response of Saccharomyces cerevisiae to D-lactic Acid

## Резюме
Модель дрожжей S. cerevisiae является популярным объектом для разных исследований. Однако, некоторые транскрипционные ответы S. cerevisae на различные вещества не до конца изучены. Авторы данного исследования изучают транскрипционный ответ штамма BY4742 на широкий диапазон концентраций DLAи сравнивают его с ответом на L-молочную кислоту (LLA) (от 0,05 до 45 мМ) . Полученные ими данные не выявили естественные дрожжевые промоторы, количественно распознающие DLA, но дают первое описание транскриптомного ответа на DLA и дают лучшего понимания ответа LLA. Так же их анализы показали, что некоторые из генов активируются кислой формой DLA, что в свою очередь раскрывает роль pH. Целью данного микропроекта является изучение и сравнение транскрипционного ответа штамма BY4742 на концентрацию D-молочной кислоты (DLA) 0,05 мМ и L-молочную кислоту (LLA).
## Методы

Использованные программы: 

  * `mafft v7.490`;
  * `R v4.1.2`;
  * `ggplot v3.5.0`.
    
Скачиваем данные: 

`fasterq-dump --threads 2 -A --progress SRR24466389; fasterq-dump --threads 2 -A --progress SRR24466390; fasterq-dump --threads 2 -A --progress SRR24466391; fasterq-dump --threads 2 -A --progress SRR24466380; fasterq-dump --threads 2 -A --progress
SRR24466381; fasterq-dump --threads 2 -A --progress SRR24466382`

Выравниваем чтения на референс:

` wget https://ftp.ensembl.org/pub/release-108/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.108.gtf.gz
Wget
https://ftp.ensembl.org/pub/release-108/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.tople
vel.fa.gz`

Распаковываем архивы:

`gunzip Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
gunzip Saccharomyces_cerevisiae.R64-1-1.108.gtf.gz`

Cоздаём индекс и готовим файл с данными сплайсинга в hisat2:

`hisat2-build Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa yeast_index
hisat2_extract_splice_sites.py Saccharomyces_cerevisiae.R64-1-1.108.gtf > yeast_splice_sites.txt`

Выравниваем с помощью hisat2 и сортируем bam-файл с помощью samtools:

`for sample in `ls *_1.fastq`; do base=$(basename $sample "_1.fastq"); hisat2 -x yeast_index --known-splicesite-infile
yeast_splice_sites.txt -p 8 -1 ${base}_1.fastq -2 ${base}_2.fastq | samtools view --threads 2 -bS | samtools sort
--threads 2 -o $base.bam; done`

Графики строили с помощью R: 
```{r}

`install.packages("BiocManager")
BiocManager::install("EnhancedVolcano")
BiocManager::install("DESeq2")`

`count_table <- read.delim("allSamples.featureCounts.txt", skip=1, row.names="Geneid")
sample_table <- data.frame(condition=c("DL", "DL", "DL", "control", "control",
"control"))`

`library(DESeq2)
ddsFullCountTable <- DESeqDataSetFromMatrix(
 countData = count_table[,6:11], colData = sample_table, design = ~ condition)
dds <- DESeq(ddsFullCountTable)
res <- results(dds)`

`library(EnhancedVolcano)
EnhancedVolcano(res, lab = rownames(res),
 x = 'log2FoldChange', y = 'pvalue',
 pCutoff=0.05, pCutoffCol = 'padj', FCcutoff = 1,
 title="Large Title", subtitle="Subtitle",
 col = c("grey30", "grey30", "grey30", "red2"),
 xlab="", ylab = bquote(~-Log[10] ~ italic(p)),
 caption="", selectLab = "", legendPosition = 'none')`

`DEGs <- res[abs(res$log2FoldChange) > 1 & res$padj < 0.05 & complete.cases(res$padj), ]
DEGs <- DEGs[order(DEGs$log2FoldChange), ]
library(openxlsx)
write.xlsx(x = DEGs, file = "DEGs_yeast.xlsx", rowNames = TRUE)`





