## PANDORA Microbiome Analysis

Title: Longitudinal Gut Microbiota Changes After Antiretroviral Therapy Initiation in People with Advanced HIV

Matthew Sedlock<sup>1</sup>\*, Brian P. Epling<sup>1</sup>\*, Kathryn McCauley<sup>2</sup>, Elisha Segrist<sup>3</sup>, Andrea Lisco<sup>1</sup>, Elizabeth Laidlaw<sup>1</sup>, Maura Manion<sup>4</sup>, Irini Sereti<sup>1</sup>

1: Division of Intramural Research, National Institute of Allergy and Infectious Diseases, National Institutes of Health, Bethesda, Maryland, USA.

2: Bioinformatics and Computational Biosciences Branch, National Institute of Allergy and Infectious Diseases, National Institutes of Health, Bethesda, Maryland, USA.

3: National Institute for General Medical Sciences, National Institutes of Health, Bethesda, Maryland, USA.

4: Division of Clinical Research, National Institute of Allergy and Infectious Diseases, National Institutes of Health, Bethesda, Maryland, USA.

*: contributed equally.

--

Code used for the above-titled manuscript can be found in the `code` folder. `figure_code_linker.csv` was used to link code to the specific output figure and table numbers/names.

Raw outputs from the Nephele pipeline are provided in `relevant_nephele_outputs`. This includes:

- `logfile.txt`: log file detailing options and log messages from the pipeline
- `OTU_table.txt`: the raw sample by ASV table
- `taxonomy_table.txt`: table of taxonomic lineages for each ASV
- `rooted_tree.nwk`: phylogenetic tree of distances between ASVs

The Nephele DADA2 pipeline was run on October 10, 2025.