program mergedata
implicit none

character(len=15) :: n,fname

integer :: i,ib,t

real(8),dimension(1:4001) :: pop1,pop2,pop3,pop,m1,m2,m3,mt

mt = 0d0
m1 = 0d0
m2 = 0d0
m3 = 0d0

do i = 1, 10
   write(n,'(i2)') i-1

   if (i <= 10) then
      write(n,'(i1)') i - 1
      fname = 'map-0'//trim(n)//'/temp.out'
   else
      write(n,'(i2)') i - 1
      fname = 'map-'//trim(n)//'/temp.out'
   end if

print *, fname

   open (11,file=fname)
   
   do ib = 1, 4001
      read(11,'(i10,4f20.9)') t, pop1(ib), pop2(ib), pop3(ib), pop(ib)
      
      mt(ib) = mt(ib) + pop(ib)

      m1(ib) = m1(ib) + pop1(ib)*pop(ib)
      m2(ib) = m2(ib) + pop2(ib)*pop(ib)
      m3(ib) = m3(ib) + pop3(ib)*pop(ib)
   end do

   close(11)
end do

open(12,file='final.log')
do ib = 1, 4001
   pop1(ib) = m1(ib)/mt(ib)
   pop2(ib) = m2(ib)/mt(ib)
   pop3(ib) = m3(ib)/mt(ib)
   
   pop(ib) = mt(ib)

   write(12,'(i10,4f18.9)') ib-1, pop1(ib), pop2(ib), pop3(ib), pop(ib)
end do
close(12)

end program mergedata
