#!/bin/bash

# This scripts uses the pre-made gridpack generation directories for the Z' to
# ttbar analysis and actually produces the gridpacks one-by-one.
#
# With Python this would be a breeze, but why not take the opportunity to hone
# our skills with Bash scripting in this fairly simple case!
#
# Run instructions:
# First give the script execution rights with "chmod +x prod_A1_gridpacks.sh"
# Then execute: ./prod_V1_gridpacks.sh
#
# Dev started on 14 Aug 2020 by juska@cern.ch

# This is the command that we need to execute each scenario with:
# time ./gridpack_generation.sh ttbarReso__inclusive__DMsimp_s_spin1__mchi_50_mZp_3000_BenchmarkA1 cards/examples/ttbarReso__inclusive__DMsimp_s_spin1__mchi_50_mZp_3000_BenchmarkA1 local ALL slc7_amd64_gcc700 CMSSW_10_6_2

cmd="time ./gridpack_generation.sh"
cmd_end="local ALL slc7_amd64_gcc700 CMSSW_10_6_2"
path=cards/examples
log_apdx=output.txt


# The scenarios already processed can be just commented out

scenarios=(
ttbarReso__inclusive__DMsimp_s_spin1__mchi_1000_mZp_1000_BenchmarkA1
#ttbarReso__inclusive__DMsimp_s_spin1__mchi_1000_mZp_2000_BenchmarkA1
#ttbarReso__inclusive__DMsimp_s_spin1__mchi_1000_mZp_2500_BenchmarkA1
#ttbarReso__inclusive__DMsimp_s_spin1__mchi_1000_mZp_3000_BenchmarkA1
#ttbarReso__inclusive__DMsimp_s_spin1__mchi_1000_mZp_3500_BenchmarkA1
#ttbarReso__inclusive__DMsimp_s_spin1__mchi_1000_mZp_4000_BenchmarkA1
#ttbarReso__inclusive__DMsimp_s_spin1__mchi_1000_mZp_500_BenchmarkA1
#ttbarReso__inclusive__DMsimp_s_spin1__mchi_100_mZp_1000_BenchmarkA1
#ttbarReso__inclusive__DMsimp_s_spin1__mchi_100_mZp_2000_BenchmarkA1
#ttbarReso__inclusive__DMsimp_s_spin1__mchi_100_mZp_2500_BenchmarkA1
#ttbarReso__inclusive__DMsimp_s_spin1__mchi_100_mZp_3000_BenchmarkA1
#ttbarReso__inclusive__DMsimp_s_spin1__mchi_100_mZp_3500_BenchmarkA1
#ttbarReso__inclusive__DMsimp_s_spin1__mchi_100_mZp_4000_BenchmarkA1
#ttbarReso__inclusive__DMsimp_s_spin1__mchi_100_mZp_500_BenchmarkA1
#ttbarReso__inclusive__DMsimp_s_spin1__mchi_10_mZp_1000_BenchmarkA1
#ttbarReso__inclusive__DMsimp_s_spin1__mchi_10_mZp_2000_BenchmarkA1
#ttbarReso__inclusive__DMsimp_s_spin1__mchi_10_mZp_2500_BenchmarkA1
#ttbarReso__inclusive__DMsimp_s_spin1__mchi_10_mZp_3000_BenchmarkA1
#ttbarReso__inclusive__DMsimp_s_spin1__mchi_10_mZp_3500_BenchmarkA1
#ttbarReso__inclusive__DMsimp_s_spin1__mchi_10_mZp_4000_BenchmarkA1
#ttbarReso__inclusive__DMsimp_s_spin1__mchi_10_mZp_500_BenchmarkA1
#ttbarReso__inclusive__DMsimp_s_spin1__mchi_50_mZp_1000_BenchmarkA1
#ttbarReso__inclusive__DMsimp_s_spin1__mchi_50_mZp_2000_BenchmarkA1
#ttbarReso__inclusive__DMsimp_s_spin1__mchi_50_mZp_2500_BenchmarkA1
#ttbarReso__inclusive__DMsimp_s_spin1__mchi_50_mZp_3000_BenchmarkA1
#ttbarReso__inclusive__DMsimp_s_spin1__mchi_50_mZp_3500_BenchmarkA1
#ttbarReso__inclusive__DMsimp_s_spin1__mchi_50_mZp_4000_BenchmarkA1
#ttbarReso__inclusive__DMsimp_s_spin1__mchi_50_mZp_500_BenchmarkA1
)

echo "Starting to produce gridpacks..."
#echo ${scenarios[@]}

for i in "${scenarios[@]}"
do
	echo Producing gridpack "for" $i
    echo started at `date`
	$cmd $i $path/$i $cmd_end &> $i"_"$log_apdx
	echo finished at `date`
	echo Run output written to $i"_"$log_apdx
	echo
	#wait This should not be needed. Let's try without first.
done

echo "All scenarios should have been produced now!"
echo "Please check that the tarballs are roughly 63 Mb in size and that"
echo "The output.txt files do not show any errors!"
echo "(The 'NLO without showering' warning is OK, as we are doing only LO)"


