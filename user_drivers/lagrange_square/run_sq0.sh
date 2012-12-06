#!/bin/bash

cd /home/mly/v327/src/ && make && make install && \
cd /home/mly/v327/user_drivers/lagrange_square/ && \
make square0 && \

./square0 --w_solver 0 --ns_solver 1 --p_solver 0 --f_solver 3 --visc 0 \
          --ang 30 --rey 100 --noel 4

#cd /home/mly/v327/src/navier_stokes/
#make && make install \
#&& cd /home/mly/v327/user_drivers/lagrange_square/ \
#&& make square3 && ./square3 --w_solver 0 --ns_solver 0 --visc Sim --ang 30 --rey 100 --noel 8 --diagw --doc_soln

#./square0.sh > square0.dat \
#&& ./square1.sh > square1.dat \
#&& ./square3.sh > square3.dat





