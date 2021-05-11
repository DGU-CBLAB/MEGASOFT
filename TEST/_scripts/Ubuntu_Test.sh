#
#!/bin/bash

$TestSettings
NumThreads=(1 2 3 4 5 6 7 8 9 10)
StrDate=$(date +'%F')
OS="Ubuntu"

$Directory
DirData="../_bin/"
DirOutput="../"$OS"/"$StrDate'/'
DirInput="../input_samples/"
DirMEGA="./MEGASOFT.o"

$DataSet
NumStudies=8
NumSNPs=10000

$MEGASOFT_Options
_method=mcmc
_samples=1000000
_seed=0
_sigma=0.05
_alpha=1
_beta=5
_thres=1.0


$Performance_Test
cp ../../MEGASOFT/MEGASOFT.o ./MEGASOFT.o
mkdir "../"$OS
mkdir $DirOutput
for nthr in ${NumThreads[@]}
do
    DirIn=$DirInput'inputMS_'$NumStudies'x'$NumSNPs'.txt'
    DirHan=$DirData"HanEskinPvalueTable.txt"
    DirOut=$DirOutput'posterior_'$NumStudies'_'$NumSNPs'_'$OS'_'$nthr".txt"
    DirLog=$DirOutput'Log_'$NumStudies'_'$NumSNPs'_'$OS'_'$nthr".log"
    
    Script=" "$DirMEGA
    Script=$Script" --input "$DirIn
    Script=$Script" --pvalue_table "$DirHan
    Script=$Script" --output "$DirOut
    Script=$Script" --log "$DirLog
    Script=$Script" --mvalue "
    Script=$Script" --mvalue_method "$_method
    Script=$Script" --mcmc_sample "$_samples
    Script=$Script" --seed "$_seed
    Script=$Script" --mvalue_prior_sigma "$_sigma
    Script=$Script" --mvalue_prior_alpha "$_alpha
    Script=$Script" --mvalue_prior_beta "$_beta
    Script=$Script" --mvalue_p_thres "$_thres
    Script=$Script" --thread "$nthr

    #echo $Script
    $Script

done