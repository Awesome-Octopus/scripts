#!/bin/bash

sed -iE 's/\s/\n/g' memb_pnts.txt
sort -n memb_pnts.txt | head -n 1
sed -iE	's/\s/\n/g' prot_pnts.txt
sort -nr prot_pnts.txt | head -n 1
  
