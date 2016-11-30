sleep 5h
for QIAOYIN in {102..141}
do
   cp haloden.f90 haloden_$QIAOYIN.f90
   sed -i -e 's/sim=5/sim='$QIAOYIN'/g' haloden_$QIAOYIN.f90
   ifort -O3 -mcmodel=medium haloden_$QIAOYIN.f90
   mv a.out a_$QIAOYIN.out
   rm haloden_$QIAOYIN.f90
#   rm a.out
   echo 'yes' $QIAOYIN
done
echo 'yes all'
qsub qsubhaloden
qsub qsubhaloden2
qsub qsubhaloden3
qsub qsubhaloden4
