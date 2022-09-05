threads=(1 2 4 6 8)

$(gcc -o lcs_par -O3 -fopenmp lcs_par.c)
$(gcc -o lcs -fopenmp lcs.c)

echo ' '

for i in 10 100 1000 20000 50000; do   
    
    echo '  entrada tamanho '$i
    echo '     tseq                   tpar     %sequencial'  

    for thread in ${threads[@]}; do

        #$(export OMP_NUM_THREADS=$thread)          
        
        $(cp A$i.in fileA.in)
        $(cp B$i.in fileB.in)

        t_s=$(./lcs)
        t_p=$(./lcs_par $thread)
                
        echo '  '$t_s'  |  '$thread' threads: '$t_p
    
    done 

    echo ' '   



done



    