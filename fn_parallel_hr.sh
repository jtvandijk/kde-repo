#! /bin/bash
source activate kde

#consistent parish
if [ "$2" -lt "1901" ];then par="conpar1851"; else par="conpar1901"; fi

#historic
cat $1 | parallel "./fn_query_hr.sh $2 $par {} | Rscript fn_kde_png.R $2 | python fn_kde_cnt.py | ./fn_write_db.sh"
