# SpermMethylomeAnalysisScripts
Custom scripts for studying mammalian sperm methylome evolution

### Dependency ###
- [MethPipe](https://github.com/smithlabcode/methpipe)
- [Epiphyte](https://github.com/smithlabcode/epiphyte)
- [adssrc](https://github.com/andrewdavidsmith/adssrc) 
- [rmap](https://github.com/smithlabcode/rmap)
- [bedtools](https://github.com/arq5x/bedtools2/)
- [UCSC Genome Browser Utilities](http://hgdownload.soe.ucsc.edu/admin/exe/): `liftOver`, `mafSpeciesSubset`
- R and R libraries `ggplot2`, `grid`, `rphast`
- Python 2.7
- [Galaxy](https://github.com/galaxyproject/galaxy/): `maf_thread_for_species.py`, `maf_to_fasta_concat.py`


### Directory structure ###
Analysis scripts assume the following directory structure

```
.
+-- scripts
+-- OriginalMethylomes
|   +-- Human_hg19
|   |   +-- Human_Sperm.meth
|   |   +-- Human_Sperm.hmr
|   |   +-- Human_Sperm.methposterior
|   |   +-- [Human_Sperm.state_level_post]
|   +-- Chimp_panTro4
|   |   +-- ...
|   +-- Gorilla_gorGor3
|   |   +-- ...
|   +-- Rhesus_rheMac3
|   |   +-- ...
|   +-- Mouse_mm10
|   |   +-- ...
|   +-- Rat_rn5
|   |   +-- ...
|   +-- Dog_canFam3
|   |   +-- ...
+-- hg19Multiz100wayMaf
|   +-- chr<1-22/M/X/Y>.maf.gz
|   +-- cpgmaps_7orth
|   |   +-- species.lst
|   |   +-- [<assembly>_cpg_map_to_hg19.cpgmap]
|   +-- cpgmaps
|   |   +-- [<assembly>_cpg_map_to_hg19.cpgmap]
+-- DNATree
|   +-- [chr<1-22>.7sp.fasta]
|   +-- [chr<1-22>.7sp.UNREST_tree.nwk]
+-- hg19Multiz7orthMethylomes
|   +-- [<species>_Sperm_hg19.state_level_post]
|   +-- [<species>_Sperm_hg19.hmr]
|   +-- [<species>_Sperm_hg19.methposterior]
|   +-- [7orth_sites.methposterior]
|   +-- breakdown
|   |   +-- 7sp_binary.nwk
|   |   +-- [methposterior.table.part<0-32>]
|   +-- events
|   |   +-- [child_parent.txt]
|   |   +-- [7sp.<species>.HYPO]
|   |   +-- [7sp.<species>.HYPO.<event>]
|   |   +-- [7orth_complement.bed]
|   |   +-- filtered
|   |   |   +-- [child_parent.txt]
|   |   |   +-- [7sp.<species>.HYPO]
|   |   |   +-- [7sp.<species>.HYPO.<event>]
|   |   |   +-- [events_number_size]
|   |   |   +-- [events_size.pdf]
|   |   |   +-- dist_to_TSS
|   |   |   |   +-- [7sp.<species>.HYPO.<event>.dist2TSS]
|   |   |   +-- shuffle_events
|   |   |   |   +-- [child_parent.txt]
|   |   |   |   +-- [branch_<species>_<gain/loss>_dist_conserved.txt]
|   |   |   |   +-- [branch_<species>_<gain/loss>_conserved_shuffled.txt]
|   |   +-- Human_HMR_by_age_group
|   |   |   +-- [Human_HYPO_since_<species>]
|   |   |   +-- [Human_HYPO_since_<species>.cpgobex]
|   |   |   +-- [Human_HYPO_since_<species>.human_sperm_roi]
+-- mm10Multiz60wayMethylomes
|   +-- [<species>_Sperm_mm10.state_level_post]
|   +-- [<species>_Sperm_mm10.hmr]
|   +-- [<species>_Sperm_mm10.methposterior]
|   +-- [7orth_sites.methposterior]
+-- OriginalOrthVar
|   +-- OrthVar<Species_assembly>
|   |   +-- [<Species>_Sperm.hmr.<ORTH/VAR>_<promoter/repeat/other>]
+-- CladeSpecific
|   +-- [Human_Sperm_HMR.since_<ancestor>]
|   +-- [Human_Sperm_HMR.since_<ancestor>.nonpromoter]
|   +-- [Human_Sperm_HMR.since_<ancestor>.promoter]
|   +-- [Mouse_Sperm_HMR.since_<ancestor>]
|   +-- [Mouse_Sperm_HMR.since_<ancestor>.nonpromoter]
|   +-- [Mouse_Sperm_HMR.since_<ancestor>.promoter]
+-- LineageMutation
|   +-- Ensembl75_Compara15
|   |   +-- raw_data
|   |   |   +-- Compara.15_eutherian_mammals_EPO.chr*_*.emf.gz
|   |   +-- processed_data
|   |   |   +-- [Compara.15_eutherian_mammals_EPO.aligned.bed]
|   |   |   +-- [Compara.15_eutherian_mammals_EPO.mutation.bed]
|   |   |   +-- [Compara.15_eutherian_mammals_EPO.mutation.bed.node_<0-12>]
|   +-- random_sample_region_0001
|   |   +-- [sampling_results.txt]
|   +-- regions
|   |   +-- <lineage>.<hmr/nonhmr>.specific
|   +-- regions_mutation
|   |   +-- [region_size.txt]
|   |   +-- [node_num_species.txt]
|   |   +-- [<lineage>.<hmr/nonhmr>.specific.txt]
+-- GC_profile
|   +-- Human_Mouse
|   |   +-- mm10_profile
|   |   |   +-- [I-II_ref_HSBD.bed.max500bp.aligned.GCprofile.1bp]
|   |   |   +-- [II-III_ref_MSBD.bed.max500bp.aligned.GCprofile.1bp]
|   |   |   +-- [III-IV_ref_HEBD.bed.max500bp.aligned.GCprofile.1bp]
|   |   +-- hg19_profile
|   |   |   +-- [I-II_ref_HSBD.bed.max500bp.aligned.GCprofile.1bp]
|   |   |   +-- [II-III_ref_MSBD.bed.max500bp.aligned.GCprofile.1bp]
|   |   |   +-- [III-IV_ref_HEBD.bed.max500bp.aligned.GCprofile.1bp]
|   +-- Dog_Mouse
|   |   +-- mm10_profile
|   |   |   +-- [I-II_ref_HSBD.bed.max500bp.aligned.GCprofile.1bp]
|   |   |   +-- [II-III_ref_MSBD.bed.max500bp.aligned.GCprofile.1bp]
|   |   |   +-- [III-IV_ref_HEBD.bed.max500bp.aligned.GCprofile.1bp]
|   |   +-- canFam3_profile
|   |   |   +-- [I-II_ref_HSBD.bed.max500bp.aligned.GCprofile.1bp]
|   |   |   +-- [II-III_ref_MSBD.bed.max500bp.aligned.GCprofile.1bp]
|   |   |   +-- [III-IV_ref_HEBD.bed.max500bp.aligned.GCprofile.1bp]
+-- TFBS
|   +-- regions
|   |   +--[<Primate/Rodent>_lineage_specific.extension.filtered.at_promoter]
|   +-- TFBS_ENCODE_Human_hg19
|   |   +-- wgEncodRegTfbsClusteredV3.bed.gz
|   |   +-- <TF>.bed
|   +-- TFBS_ENCODE_Mouse_mm9
|   |   +-- <TF>-mouse.merged.bed
|   +-- TFBS_ENCODE_Mouse_hg19
|   |   +-- <TF>-mouse.merged.bed.lift2hg19
|   +-- human_lineage_counts
|   |   +-- [<region>.HumanTFcount]
|   +-- mouse_lineage_counts
|   |   +-- [<region>.MouseTFcount]

```
- The `scripts` directory contains internal scripts called by analysis scripts
- Files in ``[]`` are generated by analysis scripts
- Link: [hg19 Multiz100way alignment](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/multiz100way/maf/)
- Link: [Ensembl75 Compara 15 species alignment](ftp://ftp.ensembl.org/pub/release-75/emf/ensembl-compara/epo_15_eutherian/)
- Link: [wgEncodRegTfbsClusteredV3.bed.gz ](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegTfbsClustered/wgEncodeRegTfbsClusteredV3.bed.gz)
- Link: [Mouse TFBS metadata](https://www.encodeproject.org/metadata/type=experiment&replicates.library.biosample.donor.organism.scientific_name=Mus%20musculus&assay_term_name=ChIP-seq&target.investigated_as=transcription%20factor&files.file_type=bed%20narrowPeak/metadata.tsv
)

