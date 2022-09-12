threads=(1 2 4 6)

echo ' '

for i in 100 1000 20000 50000; do   
    
    echo '  entrada tamanho '$i
    echo '  tpar     '  

    for thread in ${threads[@]}; do
                       
        $(cp A$i.in fileA.in)
        $(cp B$i.in fileB.in)
     
        t_p=$(./claudiogostoso $thread)
                
        echo '  '$t_s'  |  '$thread' threads: '$t_p
    
    done 

    echo ' '   

done



    