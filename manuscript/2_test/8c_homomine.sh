#!/bin/bash
qrygene=Zm00001eb001720
datadir=../1_data
perl ../../homomine \
  --qrygene $qrygene \
  --qrydir $datadir/B73 \
  --qrybase B73 \
  --tgtdir $datadir/A188 \
  --tgtbase A188 \
  --prefix B73vsA188HM
