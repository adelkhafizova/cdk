#!/bin/sh
index=0
for word in $(cat int_signatures_active.txt); do 
    if [ $(($index % 3)) = 0 ]; then
        echo '\n'
        echo $word
        echo '\\begin{figure}[H]\n\\includegraphics[width=3cm]{'$word'}\n\\end{figure}' #| >> test_aggr.tex
    else
        echo $word #| >> 'test_aggr.tex' 
    fi; 
    index=$(($index + 1)); #echo $index; 
done
