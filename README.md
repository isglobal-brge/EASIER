# EASIER

**E**w**AS**: quality control, meta-analys**I**s and **E**n**R**ichment

The EASIER package performs epigenetic wide-association study (EWAS) downstream analysis:

* **Quality control of EWAS results**
   - Folders: input and ouput
   - Configuration: array type, sample, ethnic, exclusion CpGs criteria
   - CpG filtering selection -> list of CpGs filtered and reason
   - QC with summaries -> summary SE, Beta, lambda, significatives…
   - QC with plots -> QQplot, Distribution plot, precision plot, …
   - CpG annotation and adjustment -> QCed EWAS results file
* **Meta-analysis of EWAS results (using GWAMA)**
   - Folders: input and output
   - Link to GWAMA
   - Format QCed EWAS results file
   - Run GWAMA -> EWAS meta-analysis results file
   - Meta-analysis with summaries 
   - Meta-analysis with plots -> Heterogeneity plot, distribution plots, QQ-plots, Volcano plots, Manhattan plots andForest plots,
* **Functional enrichment (pathway and molecular enrichments)**
   - GO and KEGG and MolSig
   - Pathways with Molecular Signatures Database (MSigDB)
   - ConsensusPath Data Base
   - Gene relative position
   - CpG island relative position
   - Specific for Blood : 
      - 15 ROADMAP chromatine states
      - eQTM enrichment with data from Helix : [Identification of blood autosomal cis-expression quantitative trait methylation (cis-eQTMs) in children](https://helixomics.isglobal.org/)
   - Specific for Placenta : 
      - ROADMAP chromatine states Fetal Placenta 15 and 18  
      - Partially methylated domains (PMDs)
      - Impreinted regions

