for rn in `cat ../../run2024_earlydata/run/fullgoodrunlist.list`; do
    N=`ls ../chi2/$rn | wc -l` #`ls /sphenix/tg/tg01/jets/jocl/chi2/${rn}/*chi2file* | wc -l`
    M=`ls ../chi2_bak/$rn | wc -l` #`cat chi2filesdat.txt | grep $rn | wc -l`
    if [ $N -ne $M ]; then
	echo $rn $N $M
    fi
done
