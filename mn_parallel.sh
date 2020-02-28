#! /bin/bash

#historic
for i in {1851,1861,1881,1891,1901,1911}
do
  ls input/names$i* | parallel --sshloginfile nodeslist --workdir data/kde "./fn_parallel_hr.sh {} $i"
done

#contemporary
for i in {1997..2016..1}
do
  ls input/names$i* | parallel --sshloginfile nodeslist --workdir data/kde "./fn_parallel_cr.sh {} $i"
done
