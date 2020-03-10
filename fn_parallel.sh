#! /bin/bash
source activate kde

#contemporary
cat $1 | parallel "./fn_query.sh {} | Rscript fn_kde_vec.R | python fn_kde_cnt.py | ./fn_write_db.sh"
