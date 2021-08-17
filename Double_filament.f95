!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!This Program gives output for L Vs T for a Double filament system!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

program Double_filament

    use parm                                                                                             !Calling the parm module for parameyers
    
    implicit none
    integer:: l1,l2,NG,La1(N),La2(N),N1,N2,Lap1(N),Lap2(N),j,i,&
    code,m1,m2,p1,p2,idum,citr
    real(DP):: A(7),tam,tau,t,a0,ran2,w,k1,k2,&                                                           !initialisation
    ti,tf,cms
    idum=49385                                                                                            !initialisation
    La1=0                                                                                                 !Initalisation of the filamnet 1 lattice (array)
    l2=0                                                                                                   !Initial length of the filament 1
    La2=0                                                                                                 !Initalisation of the filamnet 2 lattice (array)
    l1=0                                                                                                   !Initial length of the filament 2
    NG=0                                                                                                   !Number of GDP subunits in th pool trcaker

    open( unit = 10, file = "N=25000_2.csv")                                                                   !Opening the file
    call cpu_time(ti)                                                                                      !Calculating the CPU time

    do i =1,l1                                                                                             !Filling the Number of of GTP momoners in the filament 1 according to the given intial length
        La1(i)=1
    end do

    do i =1,l2                                                                                             !Filling the Number of of GTP momoners in the filament 2 according to the given intial length
        La2(i)=1
    end do
    citr=0                                                                                                 !Tracker for counting the number of iterations

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

    do                                                                                                      !Startin tghe main loop
        CALL SYSTEM_CLOCK(COUNT=idum)                                                                       !Calling the seed vlue for random  number genrator based on the cpu time
        N1=0                                                                                                !Tracker for counting the number of GTP monomers in the lattice (array) of the filamnet 1
        N2=0                                                                                                !Tracker for counting the number of GTP monomers in the lattice (array) of the filamnet 2
        if(l1>0) then                                                                                       !Tracking the Number of GTP only when filament 1 length is >0
            do i = 1,l1
                if (La1(i)==1)then
                    N1 = N1 + 1
                    Lap1(N1)=i
                end if
            end do
        end if
        if(l2>0) then                                                                                       !Tracking the Number of GTP only when filament 1 length is >0
            do i = 1,l2
                if (La2(i)==1)then
                    N2 = N2 + 1
                    Lap2(N2)=i
                end if
            end do
        end if
        A(2) = h*DBLE(N1)                                                                                   !Hydrolysis rate of filamnet 1
        A(5) = h*DBLE(N2)                                                                                   !Hydrolysis rate of filamnet 2
        A(7) = ne*DBLE(NG)                                                                                  !neucleotide exchange rate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (l1 > 0 .and. La1(l1)==1 )then                                                                   !Conditioning the rates with the length values
            A(3)=GT
            A(1)=q_1*(N-NG-l1-l2)
        else if(l1 > 0 .and. La1(l1)==-1)then
            A(3)=GD
            A(1)=q_2*(N-NG-l1-l2)
        else if (l1 <= 0) then
            A(2)=0
            A(3)=0
            A(1)=q_1*(N-NG-l1-l2)
        endif

        if (l2 > 0 .and. La2(l2)==1 )then
            A(6)=GT
            A(4)=q_1*(N-NG-l1-l2)
        else if(l2 > 0 .and. La2(l2)==-1)then
            A(6)=GD
            A(4)=q_2*(N-NG-l1-l2)
        else if (l2 <= 0) then
            A(5)=0
            A(6)=0
            A(4)=q_1*(N-NG-l1-l2)
        endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        a0=sum(A)
        cms=0.0_DP


        tam=ran2(idum)
        tau=(-1.0_DP)*(1/a0)*DLOG(tam)                                                                   !Using gillsepie to get the case value for selection of the recations
        t=t+tau                                                                                          !Based on the algo
        w=ran2(idum)

        do j=1,7
            code=j
            cms=cms+A(j)
            if(cms >= w*a0)exit
        end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        select case(code)                                                                               !Case function
            case(1)
                l1=l1+1                                                                                 !Growth of filamnet 1
                La1(l1)=1                                                                               !Inserting a 1 in the lattice of filament 1 symbolising the presence of  a GTP momner
            case(2)
                k1=(N1*ran2(idum)+1)                                                                    !calling a random number, based upon the total number of GTP monomer in the lattice of filamnet 1
                m1=nint(k1)                                                                             !Hydrolysis of filamnet 1
                p1=Lap1(m1)
                La1(p1)=-1                                                                              !Inserting a -1 in the lattice of filament 1 symbolising the presence of  a GDP momner
            case(3)
                La1(l1)=0                                                                               !Decay of filamnet 1
                l1=l1-1                                                                                 !adding a zero in the lattice of filamnet 1
                if(A(3)==GD)then
                    NG=NG+1                                                                             !Decay leading to an incrase of GDP monomer in the pool if the decayed subunit was a GDP
                end if
            case(4)                                                                                     !Growth of filamnet 2
                l2=l2+1                                                                                 !Inserting a 1 in the lattice of filament 1 symbolising the presence of  a GTP momner
                La2(l2)=1
            case(5)
                k2=(N2*ran2(idum)+1)                                                                    !calling a random number, based upon the total number of GTP monomer in the lattice of filamnet 2
                m2=nint(k2)                                                                             !Hydrolysis of filamnet 2
                p2=Lap2(m2)
                La2(p2)=-1                                                                              !Inserting a -1 in the lattice of filament 2 symbolising the presence of  a GDP momner
            case(6)
                La2(l2)=0                                                                               !Decay of filamnet 2
                l2=l2-1                                                                                 !adding a zero in the lattice of filamnet 2
                if(A(6)==GD)then
                    NG=NG+1                                                                             !Decay leading to an incrase of GDP monomer in the pool if the decayed subunit was a GDP
                end if    
            case(7)
                NG=NG-1                                                                                  !Neuclotide exchange
            case default
                print*,"Error in code !!! ","Code=",code

        end select

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !if(t>500.0) write(10,*) t,texp,l1
        write(10,*) t,l1,l2                                                                              !Writing the data
        citr=citr+1

        if(t >= tstop ) exit                                                            !condition for stopping the reaction

    end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
close(10)


call cpu_time(tf)
print*,"Total time_NT=50000_2=",tf-ti
print*,"Total iterations_NT=50000_2=",citr

end program Double_filament

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
   
