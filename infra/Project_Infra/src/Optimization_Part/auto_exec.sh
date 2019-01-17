#!/bin/bash
cmake .
make
./program
fichier="topython_1D.txt"
fichier2="topython_2D.txt"
fichier3="topython_1D_opti.txt"
fichier4="topython_2D_opti.txt"
fichier5="topython_2D_opti_covariance.txt"
python3 plot_error_projet.py
if [ -f $fichier ]; then
   rm topython_1D.txt
else
   echo "$fichier n'est pas present"
fi

if [ -f $fichier2 ]; then
   echo "we do not delete topython_2D.txt bcs too long to compute."
   rm topython_2D.txt
else
   echo "$fichier2 n'est pas present"
fi

if [ -f $fichier3 ]; then
   rm topython_1D_opti.txt
else
   echo "$fichier3 n'est pas present"
fi
				
if [ -f $fichier4 ]; then
   rm topython_2D_opti.txt
else
   echo "$fichier4 n'est pas present"
fi

echo "Done"

if [ -f $fichier5 ]; then
   rm topython_2D_opti_covariance.txt
else
   echo "$fichier4 n'est pas present"
fi

echo "Done"
