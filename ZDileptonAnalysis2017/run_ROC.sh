#!/bin/bash

args=("$@")

if [ $# -ne 3 ] ; then
  echo "Input 1:signal and 2:background samples and 3:region"
  exit
fi

sig=${args[0]}
bkg=${args[1]}
region=${args[2]}
vars=(
  "cosTheta1"
  "cosTheta2"
  "TOP_xl"
  "ANTITOP_xl"
  "MT2s"
  "sT_met"
  "rmin0"
  "rmin1"
)
for var in "${vars[@]}" ; do
  lines=( "sigName        $sig"
          "bkgName        $bkg"
          "varName        ${region}_${var}"
        )

  out=""
  for line in "${lines[@]}" ; do
    out="$out$line\n"
  done

  echo -e "$out" | column -t > ROC_pars.txt
  root -l -b -q 'ROC.c("ROC_pars.txt")'
done
