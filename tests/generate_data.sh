#!/bin/sh

echo "a\nb" > $1/vocab.txt

for l in "a" "b"; do
  (for i in {1..20}; do echo -n "$l " ; done) > $1/${l}x20.txt
done
