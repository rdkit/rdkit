-- github #3687 : crash on bogus fmcs input
select fmcs(m::text) from pgmol where m@>'c1ncccc1';
select fmcs(m) from pgmol where m@>'c1ncccc1';

-- github #3687 : crash on bogus fmcs input
select fmcs('q');
select fmcs('p');
select fmcs('C1CC');
