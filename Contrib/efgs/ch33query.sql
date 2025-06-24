select distinct
    md.chembl_id chid,
        cs.standard_inchi,
        md.pref_name molname,
        td.chembl_id tchid,
		td.pref_name tarname,
        pchembl_value pchembl,
        organism as organism,
		ac.activity_id,
        td.target_type
        from
        ch33.activities ac,
        ch33.molecule_dictionary md,
        ch33.assays ay,
        ch33.target_dictionary td,
        ch33.compound_structures cs,
        ch33.chembl_id_lookup cl
        where
        cl.entity_id = cs.molregno and
        cl.entity_type = 'COMPOUND' and
        ac.molregno = md.molregno and
        ac.molregno = cl.entity_id and
        ac.assay_id = ay.assay_id and
        pchembl_value > 5 and
        ac.pchembl_value is not null and td.tid = ay.tid and td.target_type in (
	'PROTEIN COMPLEX','SINGLE PROTEIN')
