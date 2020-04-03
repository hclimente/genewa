while read -r method; do grep "$method " trace.txt | cut -f2 | while read -r line; do cp -f `ls work/$line*/genes` all_outputs/`echo $method | sed 's/ (gi\,//'`_$RANDOM.txt; done; done <methods.txt
