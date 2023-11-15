##参数，PDB 文件
pdb_file=$1
out_file=$2
protein=$(basename $pdb_file .pdb)

pdb_delhetatm /home/data/sda/wt/need_pdb/$pdb_file > /home/data/sda/wt/pqr/${protein}.pdb
pdb2pqr30 --ff=PARSE --whitespace --log-level=ERROR /home/data/sda/wt/pqr/${protein}.pdb /home/data/sda/wt/pqr/${protein}.pqr

FILE=/home/data/sda/wt/pqr/${protein}.pqr
if [ -f "$FILE" ]; then
    sed "s/need_replace/\/home\/data\/sda\/\/wt\/pqr\/$protein/g" /home/wt/PCP_model/scripts/temp_apbs_in > /home/data/sda/wt/tmp/${protein}_apbs_in_manual
    apbs --output-file=/home/data/sda/wt/tmp/${protein}_out --output-format=flat /home/data/sda/wt/tmp/${protein}_apbs_in_manual
    /home/wt/software/DMS/bin/dms /home/data/sda/wt/pqr/${protein}.pdb -o /home/data/sda/wt/tmp/${protein}_dms
    python /home/wt/PCP_model/scripts/cal_patch.py -d /home/data/sda/wt/pqr/${protein}.pqr.dx -s /home/data/sda/wt/tmp/${protein}_dms -p /home/data/sda/wt/pqr/${protein}.pdb -o ${out_file}_patches.csv
else
    echo "$FILE Not exits"
fi

#cp /home/data/sda/wt/pqr/${protein}.pdb ${out_file}_model.pdb