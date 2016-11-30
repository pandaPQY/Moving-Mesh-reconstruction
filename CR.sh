pload
module load extras
export OMP_STACKSIZE=2000M
export OMP_NUM_THREADS=8
for QIAOYIN in {2..2}
do
  cp main.f90 main_$QIAOYIN.f90
  sed -i -e 's/_2/_'$QIAOYIN'/g' main.f90
  make clean
  make
  time  ./run5
  rm main.f90
  mv main_$QIAOYIN.f90 main.f90
  echo 'done' $QIAOYIN
done
echo 'done all'


