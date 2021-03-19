#!/bin/bash

cat rock_gen_spheres.config | ./rock_gen_spheres $1 

while [ $? -eq 2 ]; do
  echo "COSPROGEN cannot place spheres, retrying... "
  cat rock_gen_spheres.config | ./rock_gen_spheres $1
done
