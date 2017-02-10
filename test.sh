#!/bin/bash
# My first script

echo -e "\nTest N = 100"
./galsim 100 ./input_data/ellipse_N_00100.gal 200 1e-5 0

echo -e "\nTest N = 200"
./galsim 200 ./input_data/ellipse_N_00200.gal 200 1e-5 0

echo -e "\nTest N = 300"
./galsim 300 ./input_data/ellipse_N_00300.gal 200 1e-5 0

echo -e "\nTest N = 400"
./galsim 400 ./input_data/ellipse_N_00400.gal 200 1e-5 0

echo -e "\nTest N = 500"
./galsim 500 ./input_data/ellipse_N_00500.gal 200 1e-5 0

echo -e "\nTest N = 600"
./galsim 600 ./input_data/ellipse_N_00600.gal 200 1e-5 0

echo -e "\nTest N = 700"
./galsim 700 ./input_data/ellipse_N_00700.gal 200 1e-5 0

echo -e "\nTest N = 800"
./galsim 800 ./input_data/ellipse_N_00800.gal 200 1e-5 0

echo -e "\nTest N = 900"
./galsim 900 ./input_data/ellipse_N_00900.gal 200 1e-5 0

echo -e "\nTest N = 1000"
./galsim 1000 ./input_data/ellipse_N_01000.gal 200 1e-5 0
