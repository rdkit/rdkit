SET client_min_messages = warning;
\set ECHO none
drop extension if exists rdkit cascade;
drop table if exists mols,zinc_frags,zinc_leads,pubchem_pieces;
create extension rdkit;
\timing
\set ECHO all
RESET client_min_messages;
select id,mol_from_smiles(smiles::cstring) m into mols from raw_data;
create index molidx on mols using gist(m);
select id,mol_from_smiles(smiles::cstring) m into zinc_frags from zinc_frags_raw;
select id,mol_from_smiles(smiles::cstring) m into zinc_leads from zinc_leads_raw;
select mol_from_smiles(smiles::cstring) m into pubchem_pieces from pubchem_pieces_raw;

select count(*) from mols cross join zinc_frags qt where mols.m@>qt.m;
select count(*) from mols cross join zinc_leads qt where mols.m@>qt.m;
select count(*) from mols cross join pubchem_pieces qt where mols.m@>qt.m;

drop table if exists fps,zinc_frags_fps,zinc_leads_fps;
select id,morganbv_fp(m) mfp2 into fps from mols;
create index mfp2_idx on fps using gist(mfp2);
select id,morganbv_fp(m) mfp2 into zinc_leads_fps from zinc_leads;
select id,morganbv_fp(m) mfp2 into zinc_frags_fps from zinc_frags;

set rdkit.tanimoto_threshold=0.8;
select count(*) from fps cross join zinc_leads_fps qt where fps.mfp2%qt.mfp2;
select count(*) from fps cross join zinc_frags_fps qt where fps.mfp2%qt.mfp2;
set rdkit.tanimoto_threshold=0.6;
select count(*) from fps cross join zinc_leads_fps qt where fps.mfp2%qt.mfp2;
select count(*) from fps cross join zinc_frags_fps qt where fps.mfp2%qt.mfp2;
