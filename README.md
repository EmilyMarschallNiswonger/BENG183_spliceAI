## SpliceAI: A deep learning RNA-splicing predictor

Ryan Izadshenas, Srisarada Ramesh, Emily Marschall-Niswonger

### Introduction

SpliceAI is a deep learning–based bioinformatics tool developed by Illumina to predict how genetic variants affect RNA splicing. Many disease-causing mutations do not alter protein-coding sequences directly but instead disrupt the splicing process that produces mature mRNA. SpliceAI was created to address this challenge by identifying whether a DNA variant will weaken, remove, or create splice sites. (Jaganathan et al., 2019)

Traditional splice prediction tools rely on short sequence motifs and often miss complex or “cryptic” splicing changes. In contrast, SpliceAI analyzes a very large sequence window—up to 10,000 nucleotides around a variant—and learns splicing patterns directly from large transcriptomic datasets. This enables it to capture long-range regulatory signals and make much more accurate predictions.

Overall, SpliceAI is widely used in genomic research, clinical variant interpretation, and large-scale annotation projects. It helps scientists prioritize variants that may affect gene expression, understand potential disease mechanisms, and decide which mutations warrant experimental validation.

### Brief Review of RNA Splicing

To understand SpliceAI, it is important to review how RNA splicing works. When a gene is transcribed, the initial RNA product (pre-mRNA) contains exons — the coding regions — and introns, which must be removed. Splicing is the process by which introns are cut out and exons are joined together to produce a mature mRNA capable of being translated into protein.

<img width="1367" height="500" alt="image" src="https://github.com/user-attachments/assets/dfc7fd1f-68f9-4417-8955-56f66703dff9" />

**Figure 1** Overview of RNA splicing. Exons are joined together while introns are removed to form mature mRNA. *Image credit: Khan Academy* 

Splicing depends on specific sequence signals:

- Splice donor sites (5' splice sites) — located at the beginning of introns, typically starting with GU.
- Splice acceptor sites (3' splice sites) — located at the end of introns, typically ending with AG.
- Branch point and polypyrimidine tract — internal sequence features that help the spliceosome position itself correctly.

The spliceosome, a large RNA–protein complex, recognizes these signals and performs two precise cutting and joining reactions. Any disruption to these signals can alter how exons and introns are assembled, potentially creating abnormal mRNA transcripts.

These errors may lead to exon skipping, intron retention, activation of cryptic splice sites, or production of truncated or nonfunctional proteins. Because of this, accurate prediction of splicing changes is essential for understanding the impact of genetic variants.

### Deep Learning and Its Role in SpliceAI

Deep learning is a branch of machine learning that uses multi-layered neural networks to identify complex patterns in data. Unlike traditional algorithms that rely on hand-crafted rules, deep learning systems learn features automatically from large datasets. This makes them particularly powerful for biological problems where the underlying rules are complex or poorly understood.

A deep learning model is built from layers of interconnected “neurons,” each transforming input information in small ways. As data passes through multiple layers, the network learns increasingly abstract representations. For example, in image analysis, early layers learn edges, while deeper layers recognize shapes or objects. Training occurs by adjusting the strength of connections between neurons to minimize prediction errors.

SpliceAI applies this model to genomic sequences. According to Jaganathan (2020), instead of learning visual patterns, it learns the subtle sequence characteristics that determine where splicing occurs. The network considers long-range sequence context—thousands of nucleotides upstream and downstream of a variant—allowing it to detect regulatory signals that older splice prediction tools often miss. During training, the model was shown large numbers of real human transcripts paired with their genomic sequences, enabling it to learn where true splice donor and acceptor sites exist. While the transcriptome and genome of one human has sufficient information to train SpliceAI, the algorithm was trained with over 100 human samples to reduce noise and offered very promising results when tested. SpliceAI was trained on chromosomes 1, 3, 5, 7, and 9 and tested with all other chromosomes.

By comparing the splicing predictions for the reference sequence and the mutated sequence, SpliceAI calculates whether a variant is likely to create or disrupt splice sites. This deep learning approach gives SpliceAI unmatched accuracy for detecting cryptic splice events and predicting the functional impact of many noncoding and synonymous variants.

### A Look Into Splice AI-10K

Now that we know how Splice AI is built upon a deep learning model, we can move forward to understanding the details. 

First and foremost, it is important to note that there are a couple of Splice AI versions that exist. These include Splice AI-80nt, Splice AI-400nt, Splice AI-2k, Splice AI-10K. Of all these models, Splice AI-10k is known to be the most accurate version (Jaganathan, 2024). This is due to its capability to analyze a very large input, an input that is up to 10,000 nucleotides long. The other models analyze smaller inputs. For example, Splice AI-80nt can analyze only up to 80 nucleotides, which is far less than the range of Splice AI-10k. Eventually, Splice AI-10k is known to have a 95% top-k accuracy, which is the highest of all versions. 

Therefore, Splice AI-10k is the prevalent model and the focus of our report.

<img width="500" height="500" alt="image" src="https://assets.illumina.com/content/dam/illumina-marketing/images/genomics-research/articles/splice-ai/figure-2.jpg"/>

**Figure 2.** Showcases the top-k accuracy values for different Splice AI versions as well as other similar splice-predicting softwares. *(Figure 2, n.d.)*

#### How does Splice AI-10k work?

The model takes in a pre-mRNA sequence as an input, traverses through the sequence, and generates a set of predicted output scores for each nucleotide. The output scores tell the likelihood that that specific nucleotide is an acceptor, donor, or neither. 

The model generates these scores for each nucleotide based on its positional context. In other words, the model splits the input sequence into two regions, known as the upstream and downstream regions. These are the regions that come before and after the nucleotide of interest. Since Splice AI-10k utilizes a deep learning algorithm, it senses a nucleotide’s positional context by generating a window of analysis that covers a certain number of bases upstream and downstream. Such a cohesive analysis procedure eventually enables the model to understand long intron and exon sequences relative to the nucleotide’s position that may otherwise remain invisible to the model. This drastically increases the chances of the model making accurate predictions regarding whether the nucleotide of interest is indeed a splice acceptor or splice donor site (Jaganathan, 2024). 

<img width="500" height="500" alt="image" src="https://assets.illumina.com/content/dam/illumina-marketing/images/genomics-research/articles/splice-ai/figure-1.jpg"/>

**Figure 3.** Gives an overview of how the deep learning model works and produces output scores. *(Figure 3, n.d.)*



<img width="500" height="600" alt="image" src="https://assets.illumina.com/content/dam/illumina-marketing/images/genomics-research/articles/splice-ai/figure-3.jpg"/>

**Figure 4.** Shows an example where Splice AI-10k accurately predicts all 26 acceptor sites and all 26 donor sites for the CFTR gene, unlike another software known as MaxEntScan, which utilizes a very narrow/local window during analysis. *(Figure 4, n.d.)*

#### How does Splice AI-10k compare reference transcripts with variant transcripts?

Although Splice AI-10k is trained on a set of reference transcripts, its purpose is to be efficient enough to compare reference sequences with variant sequences. Variant sequences are simply sequences that have been mutated in some form. SNV’s are the type of variant sequences that are often used for such comparisons. SNV stands for single nucleotide variant (Canson et.al., 2023). These variants contain a point mutation that changes a single base to another base. When inputted into Splice AI-10k, the model predicts the splice-altering nature of such a variant. 

For the reference vs. variant analysis, Splice-AI 10k produces 4 summary metrics known as the acceptor gain, acceptor loss, donor gain, and donor loss (Jaganathan, 2018). These values calculate any increment or decrements in the acceptor or donor scores predicted for the reference transcript. The acceptor gain (DS_AG) computes how much an acceptor score increased in the variant. In contrast, the acceptor loss (DS_AL) computes how much an acceptor score decreased in the variant. Similarly, the donor gain (DS_DG) computes the increase in a donor score while the donor loss (DS_DL) computes the decrease in a donor score. 

The final summary metric, also known as the delta score of the variant, is generated by taking the maximum of all gain and loss values (Jaganathan, 2020). This value quantitatively expresses the splice-altering change within the variant transcript and its magnitude.

### Workflow

The user workflow is very simple, with few inputs and parameters. The user needs a .vcf (variant call format) file containing a list of variants to investigate, as well as a .vcf file name that they want their output stored. .vcf file format contains metadata as a header and a row for each variant with relevant information. More can be learned about .vcf file format [here](https://samtools.github.io/hts-specs/VCFv4.2.pdf).
<img width="1414" height="532" alt="image" src="https://github.com/user-attachments/assets/d7940bd0-ddc6-43a6-94b0-1be97e5b3e93" />
Users must also provide a genome .fa file and an annotation file. There are annotation files for the most common genomes and one just needs to specify which genome they are providing as input.

### Worked Example

Starting with an input `.vcf` file from Clinvar (2018), a database of known genetic variants and their pathogenecity as found experimentally, use bcftools to trim the file (containing the entire database) to just the gene of interest: MAPT.

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

<img width="1798" height="726" alt="image" src="https://github.com/user-attachments/assets/1c880c0c-e60b-4ab2-968d-9ee4cf5cf825" />


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

<img width="1922" height="372" alt="image" src="https://github.com/user-attachments/assets/70d5d2ce-caf8-4d4a-8126-46c8c34eb622" />

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

<img width="2058" height="366" alt="image" src="https://github.com/user-attachments/assets/a1d516a0-c729-4d1b-9ff7-005e77f7ae93" />

This leaves us with 76 variants that may be pathogenic because of their splice-altering effects. This can be a foundation for further review, or more analysis into the false negatives/positives produced by SpliceAI.

### Limitations & Current Tools

Despite being a great tool for RNA-splicing analysis, Splice AI has its own limitations. First of all, Splice AI solely relies on a deep learning algorithm that may learn patterns/features that weren’t meant to be studied by the model (Jaganathan, 2020). Therefore, the model may possibly have learned about patterns/features present in the pre-mRNA sequence that are biologically irrelevant during training. This may eventually lead to the production of overestimated or underestimated acceptor and donor scores that don’t explain any biological phenomena, and hence, affecting the delta score of the variant. Although Splice AI has been proven to be an efficient model with a high top-k accuracy, this fact makes the model susceptible to potential bias. According to Jaganathan (2020), the susceptibility to bias requires the use of extra methods like RNA-seq analysis to validate the results generated by Splice AI. 

Furthermore, Splice AI doesn’t inherently produce any visuals. This requires the use of extra visualization software like Splice AI-Visual. Splice AI-Visual helps visualize the predicted acceptor and donor sites on the reference allele. In addition, this software can show any discrepancies between the acceptor/donor score and the delta score of the variant (de Sainte Agathe et.al., 2023). 

While employing a large window of nearly 10,000 nucleotides may increase the accuracy of the model, in some cases, this may lead to inaccurate results, given that Splice AI solely relies on a deep learning algorithm that can learn patterns that are biologically irrelevant. Eventually, when studying the positional context of a nucleotide, the model could possibly make predictions that depend on bases further away from the nucleotide of interest that have no intronic or exonic relationship whatsoever. This possibility for inaccuracy offers a huge disadvantage when inputting large sequences of pre-mRNA. 

Considering the limitations of Splice AI, new tools have been introduced to better predict splicing changes in variants. One such tool is Splam, which has been developed by John Hopkins University. Splam is known to use a more realistic window-size for analysis, which is around 800 nucleotides. Furthermore, this model has been tested across different species, contrary to Splice AI. This has enhanced Splam’s model and the biological relevance of its results. According to John Hopkins Biomedical Engineering (2024), Splam has been tested on sequences derived from three distant species and the results were accurate, showcasing that variants are correctly pinpointed for their splice-altering capabilities. 

### Conclusion

Overall, it can be concluded that Splice AI is a highly efficient tool with a simple user interface and workflow that produces the necessary results in seconds. Its deep learning model makes it a highly favored tool among its predecessors. It can predict splice acceptor and splice donor sites with outstanding accuracy, unlike its predecessors. However, like any other tool, Splice AI has its own setbacks. Fortunately, modern tools, like Splam, are using Splice AI as their inspiration and build upon its algorithm to make better predictions. Therefore, Splice AI has made a lasting impact in the field of RNA-splicing analysis. 

### References

Canson, D. M., Davidson, A. L., de la Hoya, M., Parsons, M. T., Glubb, D. M., Kondrashova, O., & Spurdle, A. B. (2023). SpliceAI-10k calculator for the prediction of pseudoexonization, intron retention, and exon deletion. Bioinformatics (Oxford, England), 39(4), btad179. https://doi.org/10.1093/bioinformatics/btad179

de Sainte Agathe, J. M., Filser, M., Isidor, B., Besnard, T., Gueguen, P., Perrin, A., Van Goethem, C., Verebi, C., Masingue, M., Rendu, J., Cossée, M., Bergougnoux, A., Frobert, L., Buratti, J., Lejeune, É., Le Guern, É., Pasquier, F., Clot, F., Kalatzis, V., Roux, A. F., … Baux, D. (2023). SpliceAI-visual: a free online tool to improve SpliceAI splicing variant interpretation. Human genomics, 17(1), 7. https://doi.org/10.1186/s40246-023-00451-1

Illumina. (2018). GitHub - Illumina/SpliceAI: A deep learning-based tool to identify splice variants. GitHub; Illumina. https://github.com/Illumina/SpliceAI

Illumina. (n.d.). Figure 2 [Digital]. Retrieved December 9, 2025, from https://assets.illumina.com/content/dam/illumina-marketing/images/genomics-research/articles/splice-ai/figure-2.jpg

Illumina. (n.d.-b). Figure 3 [Digital]. Retrieved December 9, 2025, from https://assets.illumina.com/content/dam/illumina-marketing/images/genomics-research/articles/splice-ai/figure-1.jpg

Illumina. (n.d.-c). Figure 4 [Digital]. Retrieved December 9, 2025, from https://assets.illumina.com/content/dam/illumina-marketing/images/genomics-research/articles/splice-ai/figure-3.jpg

Jaganathan, K. (2020, November 12). Predicting splicing from primary sequence with deep learning. Illumina.com; Illumina. https://www.illumina.com/science/genomics-research/articles/predict-splicing-primary-sequence-deep-learning.html

Jaganathan, K., Panagiotopoulou, S. K., McRae, J. F., Darbandi, S. F., Knowles, D., Li, Y. I., Kosmicki, J. A., Arbelaez, J., Cui, W., Schwartz, G. B., Chow, E. D., Kanterakis, A., Gao, H., Kia, A., Batzoglou, S., Sanders, S. J., & Farh, K. K.-H. (2019). Predicting splicing from primary sequence with deep learning. *Cell, 176*(3), 535–548.e24. https://doi.org/10.1016/j.cell.2018.12.015

Khan Academy. (n.d.). *RNA splicing* [Diagram]. Retrieved from https://cdn.kastatic.org/ka-perseus-images/3157cfd56297f8bc312c0b53bdf9dd8c09f07063.png

Landrum, M. J., Lee, J. M., Benson, M., Brown, G. R., Chao, C., Chitipiralla, S., Gu, B., Hart, J., Hoffman, D., Jang, W., Karapetyan, K., Katz, K., Liu, C., Maddipatla, Z., Malheiro, A., McDaniel, K., Ovetsky, M., Riley, G., Zhou, G., Holmes, J. B., … Maglott, D. R. (2018). ClinVar: improving access to variant interpretations and supporting evidence. Nucleic acids research, 46(D1), D1062–D1067. https://doi.org/10.1093/nar/gkx1153


