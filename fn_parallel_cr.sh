#! /bin/bash
source activate png

#contemporary
cat $1 | parallel "./fn_query_cr.sh $2 {} | Rscript fn_kde_png.R $2 | python fn_kde_cnt.py | ./fn_write_db.sh"
