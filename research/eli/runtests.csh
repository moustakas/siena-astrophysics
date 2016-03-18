#set datafile = "dr10_ir4011_n10000.dat"
#set randfile = "dr10_randoms_ir4011_n10000.dat"
#set tag = "dr10_manera_10k"

#set datafile = "dr10_ir4011_n50000.dat"
#set randfile = "dr10_randoms_ir4011_n50000.dat"
#set tag = "dr10_manera_fiftyk"

#set datafile = "dr10_ir4011_n200000.dat"
#set randfile = "dr10_randoms_ir4011_n200000.dat"
#set tag = "dr10_manera_manera_ladofix_cartesian_voxelized_200k"

#set datafile = "dr10_ir4011_n10000.dat"
#set datafile = "test100.dat"
#set randfile = "dr10_randoms_ir4011_n10000.dat"
#set tag = "dr10_manera_manera_ladofix_cartesian_voxelized_10k"
#set tag = "dr10_manera_manera_ladofix_cartesian_10k"

################################################################################

#set datafile = "dr11_ir4011_n10000.dat"
#set randfile = "dr11_randoms_ir4011_n10000.dat"
#set tag = "10k"

#set datafile = "dr11_ir4011_n50000.dat"
#set randfile = "dr11_randoms_ir4011_n50000.dat"
#set tag = "fiftyk"

#set datafile = "dr11_ir4011_n100000.dat"
#set randfile = "dr11_randoms_ir4011_n100000.dat"
#set tag = "hundredk"

#set datafile = "dr11_ir4011_n200000.dat"
#set randfile = "dr11_randoms_ir4011_n200000.dat"
#set tag = "dr11_manera_manera_ladofix_cartesian_200k"

################################################################################
#set datafile = "dr10_cmass_selection_n100000.dat"
#set randfile = "dr10_randoms_ir4011_n100000.dat"
#set tag = "dr10_cmass_manera_hundredk"

#set datafile = "dr10_cmass_selection_n50000.dat"
#set randfile = "dr10_randoms_ir4011_n50000.dat"
#set tag = "dr10_cmass_manera_fiftyk"

#set datafile = "dr10_cmass_selection_n200000.dat"
#set randfile = "dr10_randoms_ir4011_n200000.dat"
#set tag = "dr10_cmass_manera_200k"

################################################################################
#set datafile = "dr12_cmass_selection_n10000.dat"
#set randfile = "dr11_randoms_ir4011_n10000.dat"
#set tag = "dr12_tenk"

#set datafile = "dr12_cmass_selection_n100000.dat"
#set randfile = "dr11_randoms_ir4011_n100000.dat"
#set tag = "dr12_hundredk"

# Eli's tests
#set num = "50k"
#set sampledir = "/home/elibeaudin/samples/"
set sampledir = "/home/elibeaudin/samples/"
set datafile = $sampledir"50k_weighted_north_cmass.dat"
#set datafile = $sampledir$num"_weighted_north_cmass.dat"
set randfile = $sampledir"200k_weighted_random.dat"
#set randfile = $sampledir$num"_weighted_random.dat"
set tag = "eli_voxelized_50k200k"

set cmd = "calc_2pt_pair_counts_factored_BELLIS.py"
#set cmd = "calc_2pt_pair_counts.py"

set python = `which python`
#set python = "~/anaconda/bin/python"

time $python $cmd --no-plots --outfilename "$tag"galsDD.dat $datafile $datafile --1d #>& dd"$tag".log &
time $python $cmd --no-plots --outfilename "$tag"galsDR.dat $datafile $randfile --1d #>& dr"$tag".log &
time $python $cmd --no-plots --outfilename "$tag"galsRR.dat $randfile $randfile --1d #>& rr"$tag".log &
