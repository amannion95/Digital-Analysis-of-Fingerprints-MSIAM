#!/bin/bash
cd build && cmake ../tests/ && make && mv geom_test1 geom_test2 geom_test3 gradient_descent main1 main3 rotation_and_translation_starter5 starter1 starter3 translation_starter5 ../demo && cd ../
echo && echo 
echo -e "\033[33m\033[1mExecutables are in \"../demo\" folder"
echo && echo
