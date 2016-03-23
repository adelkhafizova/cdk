#!/bin/sh
echo '\begin{center}\n'
echo '\begin{tabular}{ | c | c | c | c | c | c | c | c | }\n'
index=0
for word in $(cat cyp/statistics.csv); do 
    if [ $(($index % 8)) = 0 ]; then
        echo '\n\hline'
        echo $word
        echo '\\begin{figure}[H]\n\\includegraphics[width=3cm]{'$word'}\n\\end{figure} & ' #| >> test_aggr.tex
    else
        echo $word' & ' #| >> 'test_aggr.tex' 
    fi; 
    index=$(($index + 1)); #echo $index; 
done
echo '\end{tabular}\n'
echo '\end{center}\n'
