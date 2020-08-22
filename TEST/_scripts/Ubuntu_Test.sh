#
#!/bin/bash

$TestSettings
NumThreads=(1 3 5 10 15 20 25)
StrDate=$(date +'%F')
OS="Ubuntu"

$Directory
DirData="../_bin/"
DirOutput="../Ubuntu/"$StrDate'/'
DirInput="../input_samples/"
DirMEGA="./MEGASOFT.o"

$DataSet
NumStudies=
NumSNPs=

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

for nthr in ${NumThreads[@]}
do
    DirIn=$DirInput'inputMS_'$NumStudies'_'$NumSNPs'.txt'
    DirHan=$DirData"HanEskinPvalueTable.txt"
    DirOut=$DirOutput'posterior_'$NumStudies'_'$NumSNPs'_'$OS'_'$nthr".txt"
    DirLog=$DirOutput'Log_'$NumStudies'_'$NumSNPs'_'$OS'_'$nthr".log"
    
    Script="sudo "$DirMEGA
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