select tanimoto_sml(avalon_fp('c1ccccc1'::mol),avalon_fp('c1ccccc1'::mol));
select tanimoto_sml(avalon_fp('c1ccccc1'::mol),avalon_fp('c1ccccc1'::mol,true));
select tanimoto_sml(avalon_fp('c1ccccc1'::mol),avalon_fp('c1ccccn1'::mol));
