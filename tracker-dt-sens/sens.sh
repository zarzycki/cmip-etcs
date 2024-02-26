#!/bin/bash

# Define arrays of parameter values
TRAJRANGES=(8.0 8.5 9.0)
MINLENGTHS=(12 15 18 21 24)
MAXGAPS=(3 4 5 6 9)
MINTIMES=("24h" "30h" "36h" "48h")

# Loop over each parameter
for TRAJRANGE in "${TRAJRANGES[@]}"; do
  for MINLENGTH in "${MINLENGTHS[@]}"; do
    for MAXGAP in "${MAXGAPS[@]}"; do
      for MINTIME in "${MINTIMES[@]}"; do
        echo "Testing with TRAJRANGE=$TRAJRANGE MINLENGTH=$MINLENGTH MAXGAP=$MAXGAP MINTIME=$MINTIME"
        ./testing.sh $TRAJRANGE $MINLENGTH $MAXGAP $MINTIME
      done
    done
  done
done

