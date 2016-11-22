#!/bin/bash

folderList='D8 D16'

for elmt in $folderList:
do
	cd ${elmt}
	awk '{ if ($0 == "MT2" ) print }' ${elmt}.dcs.clipped.no_overlap.pileup  > ${elmt}.dcs.clipped.only_MT2.pileup
	python /Users/RRisques/Desktop/Duplex_Sequencing/Programs/plot_depth_by_position.py ${elmt}.dcs.clipped.no_overlap.region.pileup
	cd ..
done

	