## SpliceAI: An deep learning RNA splicing predictor

Ryan Izadshenas, Srisarada Ramesh, Emily Marschall-Niswonger

### Introduction

SpliceAI is a deep learning–based bioinformatics tool developed by Illumina to predict how genetic variants affect RNA splicing. Many disease-causing mutations do not alter protein-coding sequences directly but instead disrupt the splicing process that produces mature mRNA. SpliceAI was created to address this challenge by identifying whether a DNA variant will weaken, remove, or create splice sites.

Traditional splice prediction tools rely on short sequence motifs and often miss complex or “cryptic” splicing changes. In contrast, SpliceAI analyzes a very large sequence window—up to 10,000 nucleotides around a variant—and learns splicing patterns directly from large transcriptomic datasets. This enables it to capture long-range regulatory signals and make much more accurate predictions.

Overall, SpliceAI is widely used in genomic research, clinical variant interpretation, and large-scale annotation projects. It helps scientists prioritize variants that may affect gene expression, understand potential disease mechanisms, and decide which mutations warrant experimental validation.

### Brief Review of RNA Splicing

To understand SpliceAI, it is important to review how RNA splicing works. When a gene is transcribed, the initial RNA product (pre-mRNA) contains exons — the coding regions — and introns, which must be removed. Splicing is the process by which introns are cut out and exons are joined together to produce a mature mRNA capable of being translated into protein.

Splicing depends on specific sequence signals:

- Splice donor sites (5' splice sites) — located at the beginning of introns, typically starting with GU.
- Splice acceptor sites (3' splice sites) — located at the end of introns, typically ending with AG.
- Branch point and polypyrimidine tract — internal sequence features that help the spliceosome position itself correctly.

The spliceosome, a large RNA–protein complex, recognizes these signals and performs two precise cutting and joining reactions. Any disruption to these signals can alter how exons and introns are assembled, potentially creating abnormal mRNA transcripts.

These errors may lead to exon skipping, intron retention, activation of cryptic splice sites, or production of truncated or nonfunctional proteins. Because of this, accurate prediction of splicing changes is essential for understanding the impact of genetic variants.

### Deep Learning and Its Role in SpliceAI

Deep learning is a branch of machine learning that uses multi-layered neural networks to identify complex patterns in data. Unlike traditional algorithms that rely on hand-crafted rules, deep learning systems learn features automatically from large datasets. This makes them particularly powerful for biological problems where the underlying rules are complex or poorly understood.

A deep learning model is built from layers of interconnected “neurons,” each transforming input information in small ways. As data passes through multiple layers, the network learns increasingly abstract representations. For example, in image analysis, early layers learn edges, while deeper layers recognize shapes or objects. Training occurs by adjusting the strength of connections between neurons to minimize prediction errors.

SpliceAI applies this model to genomic sequences. Instead of learning visual patterns, it learns the subtle sequence characteristics that determine where splicing occurs. The network considers long-range sequence context—thousands of nucleotides upstream and downstream of a variant—allowing it to detect regulatory signals that older splice prediction tools often miss. During training, the model was shown large numbers of real human transcripts paired with their genomic sequences, enabling it to learn where true splice donor and acceptor sites exist. While the transcriptome and genome of one human has sufficient information to train SpliceAI, the algorithm was trained with over 100 human samples to reduce noise and offered very promising results when tested. SpliceAI was trained on chromosomes 1, 3, 5, 7, and 9 and tested with all other chromosomes.

By comparing the splicing predictions for the reference sequence and the mutated sequence, SpliceAI calculates whether a variant is likely to create or disrupt splice sites. This deep learning approach gives SpliceAI unmatched accuracy for detecting cryptic splice events and predicting the functional impact of many noncoding and synonymous variants.

### Workflow

The user workflow is very simple, with few inputs and parameters. The user needs a .vcf (variant call format) file containing a list of variants to investigate, as well as a .vcf file name that they want their output stored. .vcf file format contains metadata as a header and a row for each variant with relevant information. More can be learned about .vcf file format [here](https://samtools.github.io/hts-specs/VCFv4.2.pdf).![alt text](image-3.png)Users must also provide a genome .fa file and an annotation file. There are annotation files for the most common genomes and one just needs to specify which genome they are providing as input.

### Worked Example

Starting with an input `.vcf` file from Clinvar, a database of known genetic variants and their pathogenecity as found experimentally, use bcftools to trim the file (containing the entire database) to just the gene of interest: MAPT.

```bash
# Download ClinVar VCF (still large, ~2GB compressed)
wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz

# Extract by gene symbol using INFO annotations
bcftools view clinvar.vcf.gz -i 'GENEINFO="MAPT"' -o mapt_clinvar.vcf
```

Be sure to download the genome file that corresponds to the genome specified in the clinvar database. This gene required hg19. With the genome `.fa` file and the MAPT variants `.vcf`, SpliceAI is ready to be installed and ran:

```Bash
# installation:
pip install spliceai
pip install tensorflow
```

```Bash
# running spliceAI
spliceai -I mapt_clinvar.vcf -O output.vcf -R hg38.fa -A grch38
```

After about 30 minutes, `output.vcf` is obtained. For this analysis, we want to see what variants SpliceAI identifies as splice-altering that clinVar also identifies as pathogenic. This can be used as a validation of SpliceAI's efficacy in identifying pathogenic mutations and to provide an explanation for why a variant in clinvar is pathogenic. For this analysis, we will need pandas.

```
import pandas as pd
clinvar = pd.read_csv('output.vcf', sep="\t", comment='#')
```

![alt text](image.png)

The initial dataframe isn't very useful as all the relavent data is in one column, so more processing is required to separate the columns.

```python
#parse the INFO column for clinical significance and SpliceAI input
clinvar['CLNSIG'] = clinvar['INFO'].str.extract(r'CLNSIG=([^;]+)')
clinvar['spliceAI'] = clinvar['INFO'].str.extract(r'SpliceAI=([^;]+)')

#parse SpliceAI input and provide relevent columns (need 9 dummy columns due to formatting issues, can delete immediately)
clinvar[['ALLELE', 'SYMBOL', 'DS_AG', 'DS_AL','DS_DG','DS_DL', 'DP_AG', 'DP_AL', 'DP_DG', 'DP_DL',1,2,3,4,5,6,7,8,9]] = clinvar['spliceAI'].str.split('|', expand=True)
clinvar = clinvar.drop(columns = [1,2,3,4,5,6,7,8,9])

#typecast SpliceAI values to floats
clmns = ['DS_AG', 'DS_AL','DS_DG','DS_DL', 'DP_AG', 'DP_AL', 'DP_DG', 'DP_DL']
for name in clmns:
    clinvar[name] = pd.to_numeric(clinvar[name], errors='coerce')
```

![alt text](image-1.png)
Now the dataframe `clinvar` is ready for filtering. We will filter this dataset for all variants that clinvar has marked as 'Pathogenic' or 'Conflicting_classifications_of_pathogenicity' that also have some possible splicing impact:

```python
sig_results = clinvar[
    clinvar["CLNSIG"].isin(['Conflicting_classifications_of_pathogenicity', 'Pathogenic']) &
    (
        (clinvar['DS_AG'] > 0) |
        (clinvar['DS_AL'] > 0) |
        (clinvar['DS_DG'] > 0) |
        (clinvar['DS_DL'] > 0)
    )
].copy()
```

![alt text](image-2.png)
This leaves us with 76 variants that may be pathogenic because of their splice-altering effects. This can be a foundation for further review, or more analysis into the false negatives/positives produced by SpliceAI
