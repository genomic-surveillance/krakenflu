-- This fixture set contains a nearly complete RSV taxonomy tree. 
-- The complete tree path down to e.g. 'Human respiratory syncytial virus B' is
-- Viruses; Riboviria; Orthornavirae; Negarnaviricota; Haploviricotina; Monjiviricetes; Mononegavirales; Pneumoviridae; Orthopneumovirus; Orthopneumovirus hominis; human respiratory syncytial virus; Human respiratory syncytial virus B
-- but we don't need all of the high-level nodes here, so the tree path constructed here goes:
-- Monjiviricetes; Mononegavirales; Pneumoviridae; Orthopneumovirus; Orthopneumovirus hominis; human respiratory syncytial virus; Human respiratory syncytial virus B
INSERT INTO taxonomy_names (id,tax_id,name,name_class,unique_name) VALUES
    (1,11250,'Human orthopneumovirus','scientific name',''),
    (2,208893,'Human respiratory syncytial virus A','scientific name',''),
    (3,208895,'Human respiratory syncytial virus B','scientific name',''),
    (4,3049954,'Orthopneumovirus hominis','scientific name',''),
    (5,3050248,'Orthopneumovirus bovis','scientific name',''),
    (6,11246,'Bovine orthopneumovirus','scientific name',''),
    (7,12814,'Respiratory syncytial virus','scientific name',''),
    (8,1868215,'Orthopneumovirus','scientific name',''),
    (9,11244,'Pneumoviridae','scientific name',''),
    (10,11157,'Mononegavirales','scientific name',''),
    (11,2497574,'Monjiviricetes','scientific name','');

INSERT INTO taxonomy_nodes (tax_id,parent_tax_id,rank,embl_code,division_id,inherited_div_flag,genetic_code_id,inherited_GC_flag,mitochondrial_genetic_code_id,inherited_MGC_flag,GenBank_hidden_flag,hidden_subtree_root_flag,comments) VALUES
	(11250,3049954,'no rank','',9,1,1,1,0,1,1,0,'code compliant; specified'),
    (208893,11250,'no rank','',9,1,1,1,0,1,1,0,'code compliant; specified'),
    (208895,11250,'no rank','',9,1,1,1,0,1,1,0,'code compliant; specified'),
    (3049954,1868215,'no rank','',9,1,1,1,0,1,1,0,'code compliant; specified'),
    (1868215,11244,'no rank','',9,1,1,1,0,1,1,0,'code compliant; specified'),
    (11244,11157,'no rank','',9,1,1,1,0,1,1,0,'code compliant; specified'),
    (3050248, 1868215,'no rank','',9,1,1,1,0,1,1,0,'code compliant; specified'),
    (11246, 3050248,'no rank','',9,1,1,1,0,1,1,0,'code compliant; specified'),
    (12814, 11246,'no rank','',9,1,1,1,0,1,1,0,'code compliant; specified'),
    (11157,2497574,'no rank','',9,1,1,1,0,1,1,0,'code compliant; specified');

-- RefSeq sequences without RSV A/B labels
-- loading with tax_ids to bring the DB into state after full linkage
-- DNA sequences not real but length set to pass the min RSV size filter (15000)
INSERT INTO sequences (id,tax_id,fasta_header,dna_sequence,percent_n,seq_length,segment_number,ncbi_acc,flu_name,flu_type,flu_a_h_subtype,flu_a_n_subtype,include,is_flu,category,original_tax_id) VALUES
	(1,11250,'Human orthopneumovirus Subgroup B','ATCGACTGACTGC',0,15100,NULL,'NC_001781.1',NULL,NULL,NULL,NULL,1,0,NULL,NULL),
	(2,11250,'Human orthopneumovirus Subgroup A','ATCGACTGACTGC',0,15100,NULL,'NC_038235.1',NULL,NULL,NULL,NULL,1,0,NULL,NULL),
	(3,12814,'Respiratory syncytial virus, complete genome','ATCGACTGACTGC',0,15100,NULL,'NC_001803.1',NULL,NULL,NULL,NULL,1,0,NULL,NULL);

-- add some RSV A/B labelled sequences to simulate upload of data e.g. from Nextstrain
-- these would not have tax_id set in real data, therefore not set here either
INSERT INTO sequences (id,tax_id,fasta_header,dna_sequence,percent_n,seq_length,segment_number,ncbi_acc,flu_name,flu_type,flu_a_h_subtype,flu_a_n_subtype,include,is_flu,category,original_tax_id) VALUES
	(4,NULL,'known RSV A','ATCGACTGACTGC',0,15100,NULL,NULL,NULL,NULL,NULL,NULL,1,0,"RSV A",NULL),
	(5,NULL,'known RSV A 2','ATCGACTGACTGC',0,15100,NULL,NULL,NULL,NULL,NULL,NULL,1,0,"RSV A",NULL),
	(6,NULL,'known RSV B','ATCGACTGACTGC',0,15100,NULL,NULL,NULL,NULL,NULL,NULL,1,0,"RSV B",NULL),
-- add a sequence that is too short and should be filtered out
	(7,NULL,'short RSV A','ATCGACTGACTGC',0,14000,NULL,NULL,NULL,NULL,NULL,NULL,1,0,"RSV A",NULL);
