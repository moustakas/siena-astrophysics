set datafile = "ten_thousand_cmass.fits"
set randfile = "ten_thousand_rands.dat"
set tag = "tenk"

#set datafile = "twenty_thousand_cmass.fits"
#set randfile = "twenty_thousand_rands.dat"
#set tag = "twentyk"

#set datafile = "fifty_thousand_cmass.fits"
#set randfile = "fifty_thousand_rands.dat"
#set tag = "fiftyk"

time python calc_2pt_pair_counts.py --no-plots --outfilename "$tag"galsDD.dat $datafile $datafile --1d
time python calc_2pt_pair_counts.py --no-plots --outfilename "$tag"galsDR.dat $datafile $randfile --1d
time python calc_2pt_pair_counts.py --no-plots --outfilename "$tag"galsRR.dat $randfile $randfile --1d


