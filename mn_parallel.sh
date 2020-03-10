#! /bin/bash
ls input/names* | parallel --sshloginfile nodeslist --workdir data/kde "./fn_parallel.sh {}"
