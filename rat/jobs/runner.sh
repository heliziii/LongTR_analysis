

for i in $(seq 0 9)
do
        echo $i
        ./stat.sh 0$i
done



for i in $(seq 10 30)
do
	echo $i
	./stat.sh $i
done
