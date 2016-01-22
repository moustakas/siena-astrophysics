#set datafile = "dr11_ir4011_n10000.dat"
#set randfile = "dr11_randoms_ir4011_n10000.dat"
#set tag = "tenk"

set datafile = "dr11_ir4011_n50000.dat"
set randfile = "dr11_randoms_ir4011_n50000.dat"
set tag = "fiftyk"

echo time python calc_2pt_pair_counts.py --no-plots --outfilename "$tag"galsDD.dat $datafile $datafile --1d
echo time python calc_2pt_pair_counts.py --no-plots --outfilename "$tag"galsDR.dat $datafile $randfile --1d
echo time python calc_2pt_pair_counts.py --no-plots --outfilename "$tag"galsRR.dat $randfile $randfile --1d


