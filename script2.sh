threads=(8)

$(gcc -o lcs_par -O3 -fopenmp lcs_par.c)
#$(gcc -o lcs -fopenmp lcs.c)

#echo ' '

for i in 100 1000 20000 50000; do   
    
    echo '  entrada tamanho '$i
    echo 'tseq'  

    #for thread in ${threads[@]}; do
        
        $(cp A$i.in fileA.in)
        $(cp B$i.in fileB.in)

        for j in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do
            #t_s=$(./lcs)
            t_p=$(./lcs_par 6)
            echo $t_p
        done
           
    #done 

    echo ' '   



done



    