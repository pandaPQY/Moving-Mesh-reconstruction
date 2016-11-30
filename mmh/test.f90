program test
implicit none
real(4),allocatable::array(:,:,:,:)
real(4),allocatable::tmp(:,:,:)
allocate(array(3,512,512,512))
allocate(tmp(512,512,512))
open(31,file='./fields/0.000den0512.bin',access='stream',status='old')
write(*,*) 'reading'
read(31) tmp
close(31)
write(*,*)'mean:',sum(tmp)/512/512/512
array(1,:,:,:)=tmp
open(31,file='./fields/test.bin',access='stream',status='replace')
write(*,*) 'writing'
tmp=array(1,:,:,:)
write(31) tmp
close(31)
deallocate(array)
deallocate(tmp)
endprogram test
