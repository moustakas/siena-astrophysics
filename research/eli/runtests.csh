#set datafile = "dr10_ir4011_n10000.dat"
#set randfile = "dr10_randoms_ir4011_n10000.dat"
#set tag = "dr10_manera_tenk"

#set datafile = "dr10_ir4011_n50000.dat"
#set randfile = "dr10_randoms_ir4011_n50000.dat"
#set tag = "dr10_manera_fiftyk"

set datafile = "dr10_ir4011_n200000.dat"
set randfile = "dr10_randoms_ir4011_n200000.dat"
set tag = "dr10_manera_manera_200k"

#set datafile = "dr11_ir4011_n10000.dat"
#set randfile = "dr11_randoms_ir4011_n10000.dat"
#set tag = "tenk"

#set datafile = "dr11_ir4011_n50000.dat"
#set randfile = "dr11_randoms_ir4011_n50000.dat"
#set tag = "fiftyk"

#set datafile = "dr11_ir4011_n200000.dat"
#set randfile = "dr11_randoms_ir4011_n200000.dat"
#set tag = "dr11_manera_manera_200k"

#set datafile = "dr11_ir4011_n100000.dat"
#set randfile = "dr11_randoms_ir4011_n100000.dat"
#set tag = "hundredk"

################################################################################
#set datafile = "dr10_cmass_selection_n100000.dat"
#set randfile = "dr10_randoms_ir4011_n100000.dat"
#set tag = "dr10_cmass_manera_hundredk"

#set datafile = "dr10_cmass_selection_n50000.dat"
#set randfile = "dr10_randoms_ir4011_n50000.dat"
#set tag = "dr10_cmass_manera_fiftyk"

set datafile = "dr10_cmass_selection_n200000.dat"
set randfile = "dr10_randoms_ir4011_n200000.dat"
set tag = "dr10_cmass_manera_200k"

################################################################################
#set datafile = "dr12_cmass_selection_n10000.dat"
#set randfile = "dr11_randoms_ir4011_n10000.dat"
#set tag = "dr12_tenk"

#set datafile = "dr12_cmass_selection_n100000.dat"
#set randfile = "dr11_randoms_ir4011_n100000.dat"
#set tag = "dr12_hundredk"

time python calc_2pt_pair_counts.py --no-plots --outfilename "$tag"galsDD.dat $datafile $datafile --1d >& dd"$tag".log &
time python calc_2pt_pair_counts.py --no-plots --outfilename "$tag"galsDR.dat $datafile $randfile --1d >& dr"$tag".log &
time python calc_2pt_pair_counts.py --no-plots --outfilename "$tag"galsRR.dat $randfile $randfile --1d >& rr"$tag".log &


