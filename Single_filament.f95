!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!This Program gives output for L Vs T for a single filament system!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE numz
    IMPLICIT NONE
    INTEGER,PARAMETER:: DP=KIND(1.0D0)
    !REAL(DP),PARAMETER:: pi=3.14159265358979_DP
  END MODULE numz 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module parm

    use numz

    implicit none
    integer,parameter:: N=25000,Ttrans=0

    real(DP),parameter:: h=0.004_DP,GT=24._DP,GD=290._DP,ne=0.9_DP,q_1=0.003_DP,q_2=0.00003_DP,tstop=86400_DP
    

    
end module parm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program single_filament

    use parm                                                                            !Calling the parm module for parameyers
    
    implicit none
    integer:: texp,l1,NG,La1(N),N1,Lap1(N),j,i,code,m1,&                                !initialisation
    p1,idum,citr
    real(DP):: A(4),tam,tau,t,a0,ran2,w,k1,&                                            !Initialisation
    ti,tf,cms
    idum=49385
    La1=0                                                                               !Initalisation of the filamnet lattice (array)
    l1=0                                                                                !Initial length of the filament
    NG=0                                                                                !Number of GDP subunits in th pool trcaker

    open( unit = 10, file = "test.csv")                                                !Opening the file
    call cpu_time(ti)                                                                   !Calculating the CPU time

    do i =1,l1                                                                          !Filling the Number of of GTP momoners in the filament according to the given intial length
        La1(i)=1
    end do
    citr=0                                                                              !Tracker for counting the number of iterations

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    do                                                                                  !Startin tghe main loop
        CALL SYSTEM_CLOCK(COUNT=idum)                                                   !Calling the seed vlue for random  number genrator based on the cpu time
        N1=0                                                                            !Tracker for counting the number of GTP monomers in the lattice (array) of the filamnet 
        if(l1>0) then                                                                   !Tracking the Number of GTP only when filament length is >0
            do i = 1,l1
                if (La1(i)==1)then
                    N1 = N1 + 1
                    Lap1(N1)=i
                end if
            end do
        end if
        A(2) = h*DBLE(N1)                                                               !Hydrolysis rate
        A(4) = ne*DBLE(NG)                                                              !neucleotide exchange rate   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        if (l1 > 0 .and. La1(l1)==1 )then                                               !Conditioning the rates with the length values
            A(3)=GT
            A(1)=q_1*(N-NG-l1)
        else if(l1 > 0 .and. La1(l1)==-1)then
            A(3)=GD
            A(1)=q_2*(N-NG-l1)
        else if (l1 <= 0) then
            A(2)=0
            A(3)=0
            A(1)=q_1*(N-NG-l1)
        endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        a0=sum(A)
        cms=0.0_DP


        tam=ran2(idum)
        tau=(-1.0_DP)*(1/a0)*DLOG(tam)                                                  !Using gillsepie to get the case value for selection of the recations
        t=t+tau                                                                         !Based on the algo
        w=ran2(idum)

        do j=1,4
            code=j
            cms=cms+A(j)
            if(cms >= w*a0)exit
        end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        select case(code)                                                               !Case function
            case(1)
                l1=l1+1                                                                 !Growth 
                La1(l1)=1                                                               !Inserting a 1 in the lattice symbolising the presence of  a GTP momner
            case(2)
                k1=(N1*ran2(idum)+1)                                                    !calling a random number, based upon the total number of GTP monomer in the lattice
                m1=nint(k1)
                p1=Lap1(m1)                                                             !Hydrolysis
                La1(p1)=-1                                                              !Inserting a -1 in the lattice symbolising the presence of  a GDP momner
            case(3)
                La1(l1)=0                                                               !Decay
                l1=l1-1                                                                 !adding a zero in the lattice
                if(A(3)==GD)then
                    NG=NG+1                                                             !Decay leading to an incrase of GDP monomer in the pool if the decayed subunit was a GDP
                end if
            case(4)
                NG=NG-1                                                                 !Neuclotide exchange
            case default
                print*,"Error in code !!! ","Code=",code

        end select

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !if(t>500.0) write(10,*) t,texp,l1
        write(10,*) t,texp,l1                                                          !Writing the data
        citr=citr+1

        if(t >= tstop) exit                                            !condition for stopping the reaction

    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
close(10)


call cpu_time(tf)
print*,"Total time=",tf-ti
print*,"Total iterations=",citr

end program single_filament

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FUNCTION AND SUBROUTINES USED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!/////// Uniform Random number generators////////////////////////////////////

FUNCTION ran2(idum)
    USE numz
    IMPLICIT NONE
    REAL(DP):: ran2
    !INTEGER,INTENT(inout),OPTIONAL::idum
    INTEGER,INTENT(inout)::idum
    INTEGER,PARAMETER::IM1=2147483563,IM2=2147483399,IMM1=IM1-1
    INTEGER,PARAMETER::IA1=40014,IA2=40692,IQ1=53668
    INTEGER,PARAMETER::IQ2=52774,IR1=12211,IR2=3791   
    INTEGER,PARAMETER::NTAB=32,NDIV=1+IMM1/NTAB
    REAL(DP),PARAMETER::AM=1.0_DP/IM1,EPS=1.2e-7,RNMX=1.0_DP-EPS
    INTEGER::idum2,j,k,iv(NTAB),iy
    SAVE iv,iy,idum2
    DATA idum2/123456789/, iv/NTAB*0/, iy/0/
    IF (idum<0) THEN
       idum=MAX(-idum,1)
       idum2=idum
        DO j=NTAB+8,1,-1
           k=idum/IQ1
           idum=IA1*(idum-k*IQ1)-k*IR1
           IF (idum<0) idum=idum+IM1
           IF (j.LE.NTAB) iv(j)=idum
        ENDDO
        iy=iv(1)
     ENDIF
     k=idum/IQ1
     idum=IA1*(idum-k*IQ1)-k*IR1
     IF (idum<0) idum=idum+IM1
     k=idum2/IQ2
     idum2=IA2*(idum2-k*IQ2)-k*IR2
     IF (idum2<0) idum2=idum2+IM2
     j=1+iy/NDIV
     iy=iv(j)-idum2
     iv(j)=idum
     IF(iy.LT.1)iy=iy+IMM1
     ran2=MIN(AM*iy,RNMX)
     RETURN
   END FUNCTION ran2
   
