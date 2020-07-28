#!/usr/bin/env bash
IMAF="./oncokb-annotator/data/example_maf.txt"
OMAF="./oncokb-annotator/data/example_maf.oncokb.txt"
IF="./oncokb-annotator/data/example_fusions.txt"
OF="./oncokb-annotator/data/example_fusions.oncokb.txt"
ICNA="./oncokb-annotator/data/example_cna.txt"
OCNA="./oncokb-annotator/data/example_cna.oncokb.txt"
IC="./oncokb-annotator/data/example_clinical.txt"
OC="./oncokb-annotator/data/example_clinical.oncokb.txt"
OCPDF="./oncokb-annotator/data/example_clinical.oncokb.pdf"
TOKEN="595cb11c-9394-4f21-91e4-a636b90634ae" #OncoKB API Token
README="./oncokb-annotator/data/example_README.txt"
python ./oncokb-annotator/MafAnnotator.py -i $IMAF -o $OMAF -c $IC -b $TOKEN
python ./oncokb-annotator/FusionAnnotator.py -i $IF -o $OF -c $IC -b $TOKEN
python ./oncokb-annotator/CnaAnnotator.py -i $ICNA -o $OCNA -c $IC -b $TOKEN
python ./oncokb-annotator/ClinicalDataAnnotator.py -i $IC -o $OC -a $OMAF,$OCNA,$OF
python ./oncokb-annotator/OncoKBPlots.py -i $OC -o $OCPDF -c ONCOTREE_CODE #-n 10
python ./oncokb-annotator/GenerateReadMe.py -o $README
