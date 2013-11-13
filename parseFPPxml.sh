
for CONDITION in $(echo "S_L D_L S_R D_L" )
do
    echo $CONDITION
    for FILE in $(ls leastsq*)
    do
        grep $CONDITION $FILE | grep -o "b='[0-9.]\+" | grep -o "[0-9.]\+" >> ${CONDITION}_b.txt
        grep $CONDITION $FILE | grep -o "a='[0-9.]\+" | grep -o "[0-9.]\+" >> ${CONDITION}_a.txt
    done
done
