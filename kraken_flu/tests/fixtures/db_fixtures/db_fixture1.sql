INSERT INTO taxonomy_nodes( 
    tax_id, 
    parent_tax_id, 
    rank, 
    embl_code, 
    division_id, 
    inherited_div_flag, 
    genetic_code_id, 
    "inherited_GC_flag", 
    mitochondrial_genetic_code_id, 
    "inherited_MGC_flag", 
    "GenBank_hidden_flag", 
    hidden_subtree_root_flag, 
    comments) 
    VALUES 
    (1,1,'no rank',NULL,9,1,1,1,1,0,0,0,'this is the root node'),
    (2,1,'no rank',NULL,9,1,1,1,1,0,0,0,'child1 of root node'),
    (3,1,'no rank',NULL,9,1,1,1,1,0,0,0,'child2 of root node; this is a leaf node'),
    (4,2,'no rank',NULL,9,1,1,1,1,0,0,0,'child1 of child1; this is a leaf node')
