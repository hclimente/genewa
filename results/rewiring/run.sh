for I in `seq 0 100`
do
    rewire_tab2.nf --tab2 ../preprocessing/hint.ht_complex.hgnc.pseudo.tab2
    mv rewired.tab2 nets/rewired_$I.tab2
done
