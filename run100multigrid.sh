module unload intel
module load intel/16.0.3
export OMP_STACKSIZE=2000M
export OMP_NUM_THREADS=8
for QIAOYIN in {3..141}
do
  cp main.f90 main_$QIAOYIN.f90
  sed -i -e 's/_2/_'$QIAOYIN'/g' main.f90
  cp Makefile ../Makefile_$QIAOYIN
  sed -i -e 's/run/run'$QIAOYIN'/g' Makefile 
  make clean
  make
  cp runmultigrid runmultigrid_$QIAOYIN
  sed -i -e 's/run/run'$QIAOYIN'/g' runmultigrid_$QIAOYIN
  qsub runmultigrid_$QIAOYIN
  qstat
  rm runmultigrid_$QIAOYIN
  rm Makefile
  mv ../Makefile_$QIAOYIN Makefile
  rm main.f90
  mv main_$QIAOYIN.f90 main.f90
  echo 'done' $QIAOYIN
done
echo 'done all'


