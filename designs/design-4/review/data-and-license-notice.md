# Data and attribution

Patient-derived evidence comes from [popvax/yash](https://github.com/popvax/yash) and earlier local research identified in `../exports/source-receipt.json`. Source availability does not establish unrestricted rights to redistribute patient data or reuse third-party software. This contribution makes no additional licence grant over source material.

Reference sequences are attributed to [GENCODE](https://www.gencodegenes.org/) releases 46 and 50, with recorded provider checksums and source URLs. See [GENCODE's data-access information](https://www.gencodegenes.org/pages/data_access.html). The patient-normal coding/HLA reference is a constructed, incomplete reference rather than a comprehensive measured normal proteome.

Methods and tool attribution include OpenVax MHCflurry, the Gfeller Lab's MixMHCpred, PRIME and MixMHC2pred, Google DeepMind AlphaGenome, LinearDesign, ViennaRNA, Biopython, OR-Tools and the codon-table provider. See [the source list](literature.tsv) and [`software-receipt.json.gz`](software-receipt.json.gz). Software, data, model and commercial-use terms differ; users must check the applicable provider terms. Tool names and source attribution do not imply endorsement.

This design directory does not bundle raw sequencing reads, credentials, complete predictor distributions or trained weights. The four clinical PDFs are shared repository files under `data/reports/`, outside this design directory. Historical paths and original identifiers are retained where needed to trace evidence. Recorded software fields and author-proposed research questions do not certify clinical safety, regulatory approval or manufacturing readiness.
