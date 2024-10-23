/*  
This fixture simulates the state of the DB after running create_rsv_taxonomy and 
link_all_unlinked_sequences_to_taxonomy_nodes:
We have two sets of sequences, one simulating RefSeq and one set to simulate the result 
of linking pre-labelled RSV A/B (Nextstrain) sequences to the RSV taxonomy.  
The real-world RefSeq sequences link directly to taxonomy classifcation nodes such as 
'Human orthopneumovirus', whereas the Nextstrain sequences are linked (by kraken-flu) by 
creating new taxonomy nodes for each sequence record, which is also common practise for 
influenza isolate sequences on NCBI.  Both are valid ways of linking sequences to taxonomy 
and both are used in the NCBI taxonomy.  
*/

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

-- RefSeq RSV sequences, linked directly to the above taxonomy
INSERT INTO sequences (id,tax_id,fasta_header,dna_sequence,seq_length,segment_number,ncbi_acc,flu_name,flu_type,flu_a_h_subtype,flu_a_n_subtype,include,is_flu,category,original_tax_id) VALUES
	(1,11250,'Human orthopneumovirus Subgroup B','ATCGACTGACTGC',15100,NULL,'NC_001781.1',NULL,NULL,NULL,NULL,1,0,NULL,NULL),
	(2,11250,'Human orthopneumovirus Subgroup A','ATCGACTGACTGC',15100,NULL,'NC_038235.1',NULL,NULL,NULL,NULL,1,0,NULL,NULL),
	(3,12814,'Respiratory syncytial virus, complete genome','ATCGACTGACTGC',15100,NULL,'NC_001803.1',NULL,NULL,NULL,NULL,1,0,NULL,NULL);

-- add some RSV A/B labelled sequences to simulate upload of data e.g. from Nextstrain
-- create taxonomy nodes for each one as a child of hRSV A/B as would be the result of create_rsv_taxonomy
INSERT INTO taxonomy_names (id,tax_id,name,name_class,unique_name) VALUES
    (100,3000000,'known RSV A','scientific name',''),
    (101,3000001,'known RSV A 2','scientific name',''),
    (102,3000002,'known RSV B','scientific name',''),
    (103,3000003,'short RSV A','scientific name','');
INSERT INTO taxonomy_nodes (tax_id,parent_tax_id,rank,embl_code,division_id,inherited_div_flag,genetic_code_id,inherited_GC_flag,mitochondrial_genetic_code_id,inherited_MGC_flag,GenBank_hidden_flag,hidden_subtree_root_flag,comments) VALUES
	(3000000,208893,'no rank','',9,1,1,1,0,1,1,0,''),
    (3000001,208893,'no rank','',9,1,1,1,0,1,1,0,''),
    (3000002,208895,'no rank','',9,1,1,1,0,1,1,0,''),
    (3000003,208893,'no rank','',9,1,1,1,0,1,1,0,'');

INSERT INTO sequences (id,tax_id,fasta_header,dna_sequence,seq_length,segment_number,ncbi_acc,flu_name,flu_type,flu_a_h_subtype,flu_a_n_subtype,include,is_flu,category,original_tax_id) VALUES
	(4,3000000,'known RSV A','ATCGACTGACTGC',15100,NULL,NULL,NULL,NULL,NULL,NULL,1,0,"RSV A",NULL),
	(5,3000001,'known RSV A 2','ATCGACTGACTGC',15100,NULL,NULL,NULL,NULL,NULL,NULL,1,0,"RSV A",NULL),
	(6,3000002,'known RSV B','ATCGACTGACTGC',15100,NULL,NULL,NULL,NULL,NULL,NULL,1,0,"RSV B",NULL),
	(7,3000003,'short RSV A','ATCGACTGACTGC',14000,NULL,NULL,NULL,NULL,NULL,NULL,1,0,"RSV A",NULL);
