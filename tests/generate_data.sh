#!/bin/bash

echo "a" > $1/vocab.txt
echo "b" >> $1/vocab.txt

for l in "a" "b"; do
  (for i in {1..20}; do echo -n "$l " ; done) > $1/${l}x20.txt
done

(for i in {1..20}; do echo -n "a b " ; done) > $1/abx20.txt