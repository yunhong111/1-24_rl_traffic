#!/bin/bash
# Trie script

keyfile=~/Documents/research/PCAP/mapped_caida_trace/BT_hao_Originalkey_prefix_key

ipfile_folder=~/Documents/research/PCAP/mapped_caida_trace/

ipfileNames=("mapped-caida-1" "mapped-caida-6" "mapped-caida-11" "mapped-caida-16"
                         "mapped-caida-21"
                         "mapped-caida-26"
                         "mapped-caida-31"
                         "mapped-caida-36"
                         "mapped-caida-41"
                         "mapped-caida-46"
                         "mapped-caida-51"
                         "mapped-caida-56")

outfileName=("outfile_No1" "outfile_No2" "outfile_No3" "outfile_No4")
fingerprint=("8" "10" "12" "14" "16");

for i in $(seq 0 3)
do
	ipfileName=${ipfileNames[0]}
	ipfilePathName=$ipfile_folder$ipfileName
	outfileNameSeq=${outfileName[i]}${fingerprint[0]}
	./cuckooFilterAgtMain ${keyfile} ${ipfilePathName} ${outfileNameSeq} ${fingerprint[0]}
done



