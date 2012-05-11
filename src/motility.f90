module motility

use global
use fields

implicit none

contains

!--------------------------------------------------------------------------------
! The jump site is selected from all neighbouring sites on the basis of jump
! direction probs for both this cell and for the neighbours.
! If go = .false. the cell doesn't move
! (Try adding a drift probability in the case of traffic)
!--------------------------------------------------------------------------------
subroutine jumper(kcell,indx1,kslot1,go,kpar)
integer :: kpar,kcell,indx1(2),kslot1
logical :: go
type (cell_type), pointer :: cell
integer :: fullslots1,fullslots2,site1(3),site2(3),kslot2,stype
integer :: irel,dir1,lastdir1,indx2(2),k,z
integer :: savesite2(3,MAXRELDIR), saveslots2(MAXRELDIR)
integer :: savesite2a(3,MAXRELDIR+1), saveslots2a(MAXRELDIR+1)
real(DP) :: psum, p(MAXRELDIR+1), R, pR, psumm 

cell => cellist(kcell)
site1 = cell%site
if (site1(1) < 1) then
    write(logmsg,*) 'jumper: bad site1: ',site1
    call logger(logmsg)
    stop
endif
fullslots1 = 0
do k = 1,2
    if (indx1(k) > 0) then
        fullslots1 = fullslots1 + k
    endif
enddo
if (fullslots1 /= BOTH) then
    R = par_uni(kpar)
    if (R <= dirprob(0)) then    ! case of no jump
	    go = .false.
        return
    endif
endif
! Now we must jump (if possible)

stype = struct_type(int(cell%ctype))     ! COG_TYPE_TAG or NONCOG_TYPE_TAG

z = site1(3)
lastdir1 = cell%lastdir
psum = 0
p = 0
savesite2a = 0
saveslots2a = 0
do irel = 1,nreldir
	dir1 = reldir(lastdir1,irel)
	site2 = site1 + jumpvec(:,dir1)
	if (inside_xyz(site2)) then
	    indx2 = occupancy(site2(1),site2(2),site2(3))%indx
		if (indx2(1) >= 0) then     ! not OUTSIDE_TAG or DC
            fullslots2 = 0
            do k = 1,2
                if (indx2(k) > 0) then
                    fullslots2 = fullslots2 + k
                endif
            enddo
            if (fullslots2 == BOTH) then
                cycle
            elseif (fullslots2 /= 0) then
!                p(irel) = dirprob(irel)*GAMMA
                p(dir1) = dirprob(irel)*GAMMA
            else
!                p(irel) = dirprob(irel)
                p(dir1) = dirprob(irel)
            endif
!            psum = psum + p(irel)
!            saveslots2(irel) = fullslots2
            saveslots2a(dir1) = fullslots2
		endif
	endif
!	savesite2(:,irel) = site2
	savesite2a(:,dir1) = site2
enddo
!if (kcell == 1) write(*,'(10f7.3)') p
psum = sum(p)
if (psum == 0) then
	go = .false.
	return
else
    go = .true.
endif

! Now choose a direction on the basis of these probs p()
R = 0
R = par_uni(kpar)
pR = psum*R
psumm = 0
!do irel = 1,nreldir
do dir1 = 1,njumpdirs
   	psumm = psumm + p(dir1)
   	if (pR <= psumm) then
   		exit
   	endif
enddo
!if (irel > nreldir) irel = nreldir
!dir1 = reldir(lastdir1,irel)
if (dir1 > njumpdirs) then
    dir1 = 0
    do k = 1,njumpdirs
        if (p(k) > 0) then
            dir1 = k
            exit
        endif
    enddo
    if (dir1 == 0) then
        write(logmsg,*) 'jumper: bad dir1: ',dir1
	    call logger(logmsg)
        write(logmsg,*) 'R, psum, psumm: ',R,psum,psumm
		call logger(logmsg)
        write(logmsg,*) p
	    call logger(logmsg)
        stop
    endif
endif
!site2 = savesite2(:,irel)
site2 = savesite2a(:,dir1)
!fullslots2 = saveslots2(irel)
fullslots2 = saveslots2a(dir1)
! new code
if (diagonal_jumps) then
	dir1 = fix_lastdir(dir1,kpar)
elseif (dir1 == 0) then
	dir1 = random_int(1,6,kpar)
endif

if (fullslots2 == 0) then       ! randomly select a slot
    R = par_uni(kpar)
    if (R < 0.5) then
        kslot2 = SLOT_NUM1
    else
        kslot2 = SLOT_NUM2
    endif
elseif (fullslots2 == SLOT_NUM1) then
    kslot2 = SLOT_NUM2
elseif (fullslots2 == SLOT_NUM2) then
    kslot2 = SLOT_NUM1
else
    write(logmsg,*) 'ERROR in jumper: jump to crowded site'
    call logger(logmsg)
    stop
endif
cell%site = site2
cell%lastdir = dir1
occupancy(site2(1),site2(2),site2(3))%indx(kslot2) = kcell
occupancy(site1(1),site1(2),site1(3))%indx(kslot1) = 0
end subroutine

!-----------------------------------------------------------------------------------------
! A cell may have chemotactic susceptibility to: S1P, CCL21, OXY, CXCL13.
! How to handle multiple chemokines?
! In the limit of all CHEMO_K = 0, this should be the same as jumper().
!-----------------------------------------------------------------------------------------
subroutine chemo_jumper(kcell,indx1,kslot1,go,kpar)
integer :: kpar,kcell,indx1(2),kslot1
logical :: go
type (cell_type), pointer :: cell
integer :: fullslots1,fullslots2,site1(3),site2(3),kslot2,stype
integer :: irel,dir1,lastdir1,indx2(2),k,kr,rv(3),id, ichemo
integer :: savesite2a(3,MAXRELDIR+1), saveslots2a(MAXRELDIR+1)
real(DP) :: p(MAXRELDIR+1),psum, R, pR, psumm, stay_prob,  psave(MAXRELDIR+1)
real :: tnow, v(3), vsum(3), f
logical :: ischemo, cognate

dbug = .false.
tnow = istep*DELTA_T
cell => cellist(kcell)
cognate = associated(cell%cptr)
	
id = cell%id
site1 = cell%site
if (site1(1) < 1) then
    write(logmsg,*) 'chemo_jumper: bad site1: ',site1
    call logger(logmsg)
    stop
endif
fullslots1 = 0
do k = 1,2
    if (indx1(k) > 0) then
        fullslots1 = fullslots1 + k
    endif
enddo
stay_prob = dirprob(0)

ischemo = .false.
do k = 1,MAX_RECEPTOR
    if (receptor(k)%used .and. (cell%receptor_level(k) > 0)) then
        ischemo = .true.
        exit
    endif
enddo

vsum = 0
if (ischemo) then
!	if (cognate) then
!		write(*,*) 'receptor used: ',receptor(:)%used
!	endif
    do kr = 1,MAX_RECEPTOR
        if (receptor(kr)%used) then
			ichemo = receptor(kr)%chemokine
            f = receptor(kr)%sign*cell%receptor_level(kr)*receptor(kr)%strength
			v = chemo(ichemo)%grad(:,site1(1),site1(2),site1(3))
!	    	vsum = vsum + (f/norm(v))*v
	    	vsum = vsum + f*v
	    	if (cognate) then
	    		if (f > 0) write(*,'(3i4,2x,a,f6.1,6f8.3)') site1,receptor(kr)%name,f,v,vsum
	    		if (kr == 5) write(*,*)
	    	endif
	    endif
	enddo
	! For exit chemotaxis:
	! Need to create estimate of v() that corresponds to the direction of vsum,
	! nearest discrete location on the 3D lattice (for chemo_p(x,y,z))
	! This is an approximation to increase speed - it enables a table lookup.
	! Note that we use only the direction of v (magnitude is insignificant at this stage,
	! since it is accounted for in f)
!	v = chemo_N*vsum/norm(vsum)	! this is a quick approximation, needs checking!!!!
    ! Note: need to ensure that chemo_N is set.  How big?
    if (norm(vsum) > 0) then
		f = min(1.0,norm(vsum))     ! Note: f is in (0-1)
		v = vsum/norm(vsum)
		rv = chemo_N*v
	else
		f = 0
		ischemo = .false.
	endif
    stay_prob = dirprob(0)
    stay_prob = (1-f)*stay_prob
else
    stay_prob = dirprob(0)
endif

if (fullslots1 /= BOTH) then
    R = par_uni(kpar)
    if (R <= stay_prob) then    ! case of no jump
	    go = .false.
        return
    endif
endif
! Now we must jump (if possible)

stype = struct_type(int(cell%ctype))     ! COG_TYPE_TAG or NONCOG_TYPE_TAG

! Compute jump probabilities in the absence of chemotaxis
site1 = cell%site
lastdir1 = cell%lastdir
p = 0
savesite2a = 0
saveslots2a = 0
do irel = 1,nreldir
	dir1 = reldir(lastdir1,irel)
	site2 = site1 + jumpvec(:,dir1)
	if (inside_xyz(site2)) then
	    indx2 = occupancy(site2(1),site2(2),site2(3))%indx
		if (indx2(1) >= 0) then     ! not OUTSIDE_TAG or DC
            fullslots2 = 0
            do k = 1,2
                if (indx2(k) > 0) then
                    fullslots2 = fullslots2 + k
                endif
            enddo
            if (fullslots2 == BOTH) then
                cycle
            elseif (fullslots2 /= 0) then
                p(dir1) = dirprob(irel)*GAMMA
            else
                p(dir1) = dirprob(irel)
            endif
            saveslots2a(dir1) = fullslots2
		endif
	endif
	savesite2a(:,dir1) = site2
enddo
if (sum(p) == 0) then
    go = .false.
    return
endif

if (ischemo) then
	psave = p
!	call chemo_probs(p,v,f)
	call chemo_probs_pre(p,rv,f)     ! this is the precomputed version
endif
psum = sum(p)

if (psum == 0) then
	go = .false.
	return
else
    go = .true.
endif

! Now choose a direction on the basis of these probs p()
R = par_uni(kpar)
pR = psum*R
psumm = 0
do dir1 = 1,njumpdirs
   	psumm = psumm + p(dir1)
   	if (pR <= psumm) then
   		exit
   	endif
enddo
if (dir1 > njumpdirs) then
    dir1 = 0
    do k = 1,njumpdirs
        if (p(k) > 0) then
            dir1 = k
            exit
        endif
    enddo
endif
site2 = savesite2a(:,dir1)

fullslots2 = saveslots2a(dir1)
! new code
if (diagonal_jumps) then
	dir1 = fix_lastdir(dir1,kpar)
elseif (dir1 == 0) then
	dir1 = random_int(1,6,kpar)
endif

if (fullslots2 == 0) then       ! randomly select a slot
    R = par_uni(kpar)
    if (R <= 0.5) then
        kslot2 = SLOT_NUM1
    else
        kslot2 = SLOT_NUM2
    endif
elseif (fullslots2 == SLOT_NUM1) then
    kslot2 = SLOT_NUM2
elseif (fullslots2 == SLOT_NUM2) then
    kslot2 = SLOT_NUM1
else
    write(logmsg,*) 'ERROR in jumper: jump to crowded site'
	call logger(logmsg)
    stop
endif
cell%site = site2
cell%lastdir = dir1
occupancy(site2(1),site2(2),site2(3))%indx(kslot2) = kcell
occupancy(site1(1),site1(2),site1(3))%indx(kslot1) = 0

end subroutine

!-----------------------------------------------------------------------------------------
! In this version the parallel section has been made into a subroutine.
! Try varying the sweep order, i.e. 0,1 and 1,0
!-----------------------------------------------------------------------------------------
subroutine mover(ok)
logical :: ok
integer :: kpar=0, sweep, nsweeps, sweep1, sweep2, dsweep

if (Mnodes > 1) then
!DEC$ IF .NOT. DEFINED (_OPENMP)
!    stop
!DEC$ ENDIF
endif

if (Mnodes == 1) then
    nsweeps = 1
else
    nsweeps = 2
endif

if (mod(istep,2) == 0) then
	sweep1 = 0
	sweep2 = nsweeps-1
	dsweep = 1
else
	sweep2 = 0
	sweep1 = nsweeps-1
	dsweep = -1
endif
!do sweep = 0,nsweeps-1
do sweep = sweep1,sweep2,dsweep
	if (Mnodes > 1) then
	!$omp parallel do
		do kpar = 0,Mnodes-1
			call par_mover(sweep,kpar)
		enddo
	!$omp end parallel do
	else
		call par_mover(sweep,kpar)
	endif
enddo
ok = .true.
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine par_mover(sweep,kpar)
integer :: sweep
integer :: kpar
integer :: site1(3), kcell, indx(2), slot, z, slice, site(3), stage, region
logical :: go
type(cell_type), pointer :: cell
integer :: z_lo,z_hi

slice = sweep + 2*kpar
do kcell = 1,nlist
	if (Mnodes == 1) then
	    z_lo = 1
	    z_hi = NZ
	else
	    z_lo = zoffset(slice) + 1
	    z_hi = zoffset(slice+1)
	endif
    cell => cellist(kcell)
    if (cell%ID == 0) cycle             ! skip gaps in the list
    if (associated(cell%cptr)) then
!		call get_stage(cell%cptr,stage,region)
		stage = get_stage(cell%cptr)
		region = get_region(cell%cptr)
		if (region /= FOLLICLE) cycle
	endif
    if (cell%step == istep) cycle
    site1 = cell%site
    z = site1(3)
    if (zdomain(z) /= slice) cycle      ! not in the slice for this processor
    indx = occupancy(site1(1),site1(2),site1(3))%indx
    if (indx(1) < 0) then	! outside or DC
        write(logmsg,*) 'Error: par_mover: OUTSIDE_TAG or DC: ',kcell,site1,indx
		call logger(logmsg)
        stop
    endif
    if (kcell == indx(1)) then
        slot = 1
    elseif (kcell == indx(2)) then
        slot = 2
    else
        write(logmsg,'(a,6i8)') 'Error: par_mover: bad indx: ',kcell,site1,indx
		call logger(logmsg)
        stop
    endif
	if (use_chemotaxis) then
		call chemo_jumper(kcell,indx,slot,go,kpar)
	else
		call jumper(kcell,indx,slot,go,kpar)
	endif
    cell%step = istep
enddo
end subroutine


!-----------------------------------------------------------------------------------------
! In this version the parallel section is in this subroutine, requiring the use of
! a PRIVATE clause.
! Note that if cell was not declared as a pointer, the FIRSTPRIVATE(cell) clause is
! needed to avoid an error when cell is assigned.
!-----------------------------------------------------------------------------------------
subroutine mover2
integer :: kpar
integer :: site1(3), kcell, indx(2), slot, x, xlocal
logical :: go
type(cell_type),pointer :: cell
integer :: x_lo,x_hi,sweep
integer :: i, nslice, sum, sump, nsweeps, cnt
integer, save :: xlim(0:8)
integer :: xtotal(0:8)
integer, allocatable :: xcount(:)

if (istep == 1) then    ! must be executed when the blob size changes
    allocate(xcount(NX))
    xcount = 0
    do kcell = 1,nlist
        cell => cellist(kcell)
        if (cell%ID == 0) cycle             ! skip gaps in the list
        i = cellist(kcell)%site(1)
        xcount(i) = xcount(i) + 1
    enddo
    nslice = NBcells/(2*Mnodes)
    i = 0
    sum = 0
    sump = 0
    i = 1
    do x = 1,NX
        sum = sum + xcount(x)
        if (sum > i*nslice) then
            xlim(i) = x
            xtotal(i) = sum - sump
            sump = sum
            i = i+1
            if (i == 2*Mnodes) exit
        endif
    enddo
    xlim(2*Mnodes) = NX
    xtotal(2*Mnodes) = NBcells - sum
    deallocate(xcount)
endif

if (Mnodes > 1) then
!DEC$ IF .NOT. DEFINED (_OPENMP)
    stop
!DEC$ ENDIF
endif

if (Mnodes == 1) then
    nsweeps = 1
else
    nsweeps = 2
endif

!write(*,*) 'start sweep loop'
do sweep = 0,nsweeps-1
!write(*,*) 'sweep: ',sweep

!$omp parallel PRIVATE(kcell,kpar,cell,x_lo,x_hi,site1,xlocal,indx,slot,go,cnt)
kpar = omp_get_thread_num()
!write(*,*) 'kpar: ',kpar
do kcell = 1,nlist
	if (Mnodes == 1) then
	    x_lo = 1
	    x_hi = NX
	else
	    x_lo = xlim(sweep+2*kpar) + 1
	    x_hi = xlim(sweep+2*kpar+1)
	endif
    cell => cellist(kcell)
    if (cell%ID == 0) cycle             ! skip gaps in the list
    if (cell%step == istep) cycle
    site1 = cell%site
    xlocal = site1(1)
    if (xlocal < x_lo .or. xlocal > x_hi) cycle      ! not in the slice for this processor
    indx = occupancy(site1(1),site1(2),site1(3))%indx
    if (indx(1) < 0) then	! outside or DC
        write(*,*) 'stage1: OUTSIDE_TAG: ',kcell,site1,indx
        stop
    endif
    if (kcell == indx(1)) then
        slot = 1
    elseif (kcell == indx(2)) then
        slot = 2
    else
        write(*,'(a,7i6)') 'ERROR: stage1: bad indx: ',kcell,site1,indx
        stop
    endif
    call jumper(kcell,indx,slot,go,kpar)
    cell%step = istep
enddo
!$omp end parallel
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine compute_Cm(tagsite,ntagged,nvar0,nvar,dt,Cm)
integer :: tagsite(3,ntagged,0:nvar)
integer :: ntagged,nvar0,nvar
real :: dt,Cm
integer :: k,j,d(3),d2,d2sum
real, allocatable :: r2mean(:)

allocate(r2mean(nvar))

do k = 1,nvar
    d2sum = 0
    do j = 1,ntagged
        d = tagsite(:,j,k) - tagsite(:,j,0)
        d2 = d(1)*d(1) + d(2)*d(2) + d(3)*d(3)
        d2sum = d2sum + d2
    enddo
    r2mean(k) = d2sum/ntagged
enddo
write(*,'(10f6.1)') r2mean
call bestfit(r2mean,nvar0,nvar,dt,Cm)
Cm = Cm*DELTA_X*DELTA_X
!write(*,*) 'Cm: ',Cm

deallocate(r2mean)

end subroutine

!--------------------------------------------------------------------------------------
! Use least squares to find the best straight line fit to r2mean() from n1 to n2.
! In fact we might need to reduce the range of points.
! We need to use a range within which the slope (dr2/dt) doesn't vary too much from
! an "average" value.  Note that the larger the persistence parameter rho the longer
! it takes for the r2 plot to become linear.
!--------------------------------------------------------------------------------------
subroutine bestfit(r2mean,n1,n2,dt,Cm)
real :: r2mean(:),dt,Cm
integer :: n1,n2
integer :: n,i
real :: sumt,sumt2,sumy,sumyt,t,y,a,b

n = n2 - n1 + 1		! number of data points

sumt = 0
sumt2 = 0
sumy = 0
sumyt = 0
do i = n1,n2
	t = i*dt
	y = r2mean(i)
	sumt = sumt + t
	sumt2 = sumt2 + t*t
	sumy = sumy + y
	sumyt = sumyt + y*t
enddo
a = (n*sumyt - sumy*sumt)/(n*sumt2 - sumt*sumt)		! slope of line
b = (sumy - a*sumt)/n								! intercept
Cm = a/6
write(*,*) 'a,b: ',a,b

end subroutine

!---------------------------------------------------------------------
! Need to precompute array reldir(:,:)
! FOR NRELDIR = 6 (diagonal_jumps = .false.)
! The directions 1 - 6 are defined according to the neumann array,
! i.e. 1 is -x, 2 is +x, 3 is -y, 4 is +y, 5 is -z, 6 is +z
! The relative directions are numbered as follows:
! irel = 1 gives the same direction as lastdir
! irel = 6 gives the opposite direction to lastdir
! irel = 2,3,4,5 cover the remaining directions - order not important
! at the moment since all directions normal to lastdir are equally likely
! FOR NRELDIR = 17 (diagonal_jumps = .true.)
! The reldir(:,:) array uses the same set of 6 previous jump directions,
! restricted to the axes.  Diagonal previous jumps are accommodated by
! making a random selection of an axis direction from either 2 or 3
! possibilities (depending on whether it was a D2 or D3 jump).
! The possible jumps, now including diagonal moves, are stored in jumpvec(:)
! The jumpvec(:) entries are ordered with z varying fastest, then y, then x,
! each ranging -1, 0, +1.
! For dirprob(:) calculation:
! There are 3 groups of jumps, with sets of probs
! Group 1: Q(1) = P(1), D1 jump in the same as last direction
! Group 2: Q(2) = P(2)+P(3)+P(4)+P(5) (D2 jumps) + P(6)+P(7)+P(8)+P(9) (D3 jumps)
! Group 3: Q(3) = P(10)+P(11)+P(12)+P(13) (D1 jumps) + P(14)+P(15)+P(16)+P(17) (D2 jumps)
! i.e. Q(1) = prob of no direction change, Q(2) = prob of direction change < 90 deg,
! Q(3) = prob of direction change = 90 deg.
! Setting D1 jumps (von Neumann) to 1, the D2 jumps have length L2 = sqrt(2)
! and the D3 jumps have length L3 = sqrt(3)
! Therefore jump probs must be scaled appropriately to give symmetry within a set
! In other words:
!		P(i)/P(j) = L2/L3 for j=6,7,8,9 j=2,3,4,5
! and	P(i)/P(j) = 1/L2  for i=14,15,16,17 j=10,11,12,13
! Then choosing P(2) to represent D2 jump prob in Group 2, and P(10) to
! represent D1 jump prob in Group 3:
!		Q(2) = 4(1 + L2/L3).P(2)
!		Q(3) = 4(1 + 1/L2).P(10)
! It is always the case that Q(1) + Q(2) + Q(3) = BETA
!
! When RHO = 1, Q(1) = BETA = p1, Q(2) + Q(3) = 0
!
! When RHO = 0, P(10) = P(1), and P(2) = P(1)/L2
! therefore:	Q(3) = 4(1 + 1/L2).P(1)
! and			Q(2) = 4(1 + L2/L3).P(1)/L2
! giving:		P(1).[1 + 4(1 + L2/L3)/L2 + 4(1 + 1/L2)] = BETA
!				P(1) = BETA/[1 + 4(1 + L2/L3)/L2 + 4(1 + 1/L2)] = p0
! and:			Q(3)/Q(2) = (1 + L2)/(1 + L2/L3) = q320
!
! Now for 0 < RHO < 1
! first set Q(1) = P(1) to vary linearly with RHO from
! the value at 0 (p0) to the value at 1 (p1)
!	P(1) = p0 + RHO*(p1-p0)
! Now we know that Q(2) + Q(3) = BETA - Q(1)
! and the third parameter ALPHA is used to determine how Q(3)/Q(2)
! varies with RHO, by making it change linearly from q23 at RHO=0
! to ALPHA*q23 at RHO = 1
! Now we have all the info to solve for Q(1), Q(2), Q(3) and all P(:)
!
! Revised Moore model.
! Now allow reverse directions, but disallow D3 diagonal jumps, giving
! 18 possible directions.
! In the preferred direction the surface is a prolate spheroid, with
! a = 1, b = 1-RHO.
! In the reverse direction the surface is either a sphere (a = b) or
! an oblate spheroid (a = b^2).
!---------------------------------------------------------------------
subroutine make_reldir
integer :: lastdir,irel,k,ix,iy,iz,i,site(3)
integer :: reldir18(6,18),reldir26(6,26)
real(DP) :: qsum,E, theta, p(8), psum, dd

if (MODEL == NEUMANN_MODEL) then
	diagonal_jumps = .false.
	njumpdirs = 6
	nreldir = 6
	reldir = 0
	do lastdir = 1,6
		reldir(lastdir,1) = lastdir
		if (mod(lastdir,2) == 0) then
			reldir(lastdir,6) = lastdir - 1
		else
			reldir(lastdir,6) = lastdir + 1
		endif
		irel = 2
		do k = 1,6
			if (k /= reldir(lastdir,1) .and. k /= reldir(lastdir,6)) then
				reldir(lastdir,irel) = k
				irel = irel + 1
			endif
		enddo
	enddo
	jumpvec(:,1:6) = neumann(:,1:6)

else
	if (MODEL == MOORE18_MODEL) then
    	nreldir = 18
	elseif (MODEL == MOORE26_MODEL) then
    	nreldir = 26
    endif
    njumpdirs = 27
	diagonal_jumps = .true.
	k = 0
	do ix = -1,1
		do iy = -1,1
			do iz = -1,1
				k = k+1
				jumpvec(:,k) = (/ix,iy,iz/)
			enddo
		enddo
	enddo

! Data for revised Moore model.  D3 jumps are excluded, but reverse directions are allowed.
!	reldir(1,:) = (/  5,  2, 4, 6, 8, 11,13,15,17, 10,12,16,18, 22,24,26,20, 23 /)	! -x
!	reldir(2,:) = (/ 23, 20,22,24,26, 11,13,15,17, 10,12,16,18,  2, 4, 6, 8,  5 /)	! +x
!	reldir(3,:) = (/ 11,  2,10,12,20,  5,13,15,23,  4, 6,22,24,  8,16,18,26, 17 /)	! -y
!	reldir(4,:) = (/ 17,  8,16,18,26,  5,13,15,23,  4, 6,22,24,  2,10,12,20, 11/)	! +y
!	reldir(5,:) = (/ 13,  4,10,16,22,  5,11,17,23,  2, 8,20,26,  6,12,18,24, 15 /)	! -z
!	reldir(6,:) = (/ 15,  6,12,18,24,  5,11,17,23,  2, 8,20,26,  4,10,16,22, 13 /)	! +z

! Data for revised Moore18 model.  D3 jumps are excluded, but reverse directions are allowed.
	reldir18(1,:) = (/  5,  2, 4, 6, 8, 11,13,15,17, 10,12,16,18, 22,24,26,20, 23 /)	! -x
	reldir18(2,:) = (/ 23, 20,22,24,26, 11,13,15,17, 10,12,16,18,  2, 4, 6, 8,  5 /)	! +x
	reldir18(3,:) = (/ 11,  2,10,12,20,  5,13,15,23,  4, 6,22,24,  8,16,18,26, 17 /)	! -y
	reldir18(4,:) = (/ 17,  8,16,18,26,  5,13,15,23,  4, 6,22,24,  2,10,12,20, 11/)	! +y
	reldir18(5,:) = (/ 13,  4,10,16,22,  5,11,17,23,  2, 8,20,26,  6,12,18,24, 15 /)	! -z
	reldir18(6,:) = (/ 15,  6,12,18,24,  5,11,17,23,  2, 8,20,26,  4,10,16,22, 13 /)	! +z

! Data for revised Moore26 model.  D3 jumps are excluded, but reverse directions are allowed.
	reldir26(1,:) = (/  5,  2, 4, 6, 8,  1, 3, 7, 9, 11,13,15,17, 10,12,16,18, 19,21,27,25, 20,22,24,26, 23 /)	! -x
	reldir26(2,:) = (/ 23, 20,22,24,26, 19,21,27,25, 11,13,15,17, 10,12,16,18,  1, 3, 7, 9,  2, 4, 6, 8,  5 /)	! +x
	reldir26(3,:) = (/ 11,  2,10,12,20,  1, 3,19,21,  5,13,15,23,  4, 6,22,24,  7, 9,25,27,  8,16,18,26, 17 /)	! -y
	reldir26(4,:) = (/ 17,  8,16,18,26,  7, 9,25,27,  5,13,15,23,  4, 6,22,24,  1, 3,19,21,  2,10,12,20, 11 /)	! +y
	reldir26(5,:) = (/ 13,  4,10,16,22,  1, 7,19,25,  5,11,17,23,  2, 8,20,26,  3, 9,21,27,  6,12,18,24, 15 /)	! -z
	reldir26(6,:) = (/ 15,  6,12,18,24,  3, 9,21,27,  5,11,17,23,  2, 8,20,26,  1, 7,19,25,  4,10,16,22, 13 /)	! +z

    if (MODEL == MOORE18_MODEL) then
        reldir(1:6,1:18) = reldir18
    elseif (MODEL == MOORE26_MODEL) then
        reldir = reldir26
    endif

endif
call compute_dirprobs
!write(logmsg,'(a,7f6.3)') dirprob(0:nreldir)
!call logger(logmsg)
qsum = 0
do i = 0,nreldir
    qsum = qsum + dirprob(i)
enddo
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine compute_dirprobs
real(DP) :: L2,L3,b,c,e,d(MAXRELDIR),R(MAXRELDIR)

if (MODEL == NEUMANN_MODEL) then
	diagonal_jumps = .false.
	nreldir = 6
! New Neumann model (with reverses)
	b = 1-RHO
	b = b*b
	e = BETA/(1 + 4*b + b**2)
	dirprob(0) = 1 - BETA	! this is the prob of no jump
	dirprob(1) = e
	dirprob(2:5) = e*b
	dirprob(6) = e*b**2
else
    ! New prolate spheroid approach
	diagonal_jumps = .true.
    if (MODEL == MOORE18_MODEL) then
	    nreldir = 18
	    L2 = sqrt(2.0d0)
	    b = 1 - RHO
	    b = b*b		! revised version
	    c = b*b		! oblate spheroid on reverse side
	    d(1) = 1
	    R(1) = 1
	    d(2:5) = L2
	    R(2:5) = L2*b/sqrt(b*b+1)
	    d(6:9) = 1
	    R(6:9) = b
	    d(10:13) = L2
	    R(10:13) = b
	    d(14:17) = L2
	    R(14:17) = L2*b*c/sqrt(b*b+c*c)
	    d(18) = 1
	    R(18) = c
	    e = BETA/(1 + 4*b/sqrt(b*b+1) + 4*b + 4*b/L2 + 4*b*c/sqrt(b*b+c*c) + c)
	    dirprob(0) = 1 - BETA
	    dirprob(1:nreldir) = e*R(1:nreldir)/d(1:nreldir)
    elseif (MODEL == MOORE26_MODEL) then
    	nreldir = 26
	    L2 = sqrt(2.0d0)
	    L3 = sqrt(3.0d0)
	    b = 1 - RHO
	    b = b*b		! revised version
	    c = b*b		! oblate spheroid on reverse side
	    d(1) = 1
	    R(1) = 1
	    d(2:5) = L2
	    R(2:5) = L2*b/sqrt(b*b+1)
	    d(6:9) = L3
	    R(6:9) = L2*b/sqrt(b*b+1)
	    d(10:13) = 1
	    R(10:13) = b
	    d(14:17) = L2
	    R(14:17) = b
	    d(18:21) = L3
	    R(18:21) = L2*b*c/sqrt(b*b+c*c)
	    d(22:25) = L2
	    R(22:25) = L2*b*c/sqrt(b*b+c*c)
	    d(26) = 1
	    R(26) = c
	    e = BETA/(1 + 4*(1+L2/L3)*b/sqrt(b*b+1) + 4*b + 4*b/L2 + 4*(1+L2/L3)*b*c/sqrt(b*b+c*c) + c)
	    dirprob(0) = 1 - BETA
	    dirprob(1:nreldir) = e*R(1:nreldir)/d(1:nreldir)
    endif
endif
end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine make_probvectors
integer :: icase, i, k, dir, jump(3)
real :: vec(3), veclen, scale
real :: xang = 0.3, yang = -0.15, zang = -0.5  ! radians
real :: scale_M18 = 10  ! for M18
real :: scale_N = 6     ! for N
real :: arrow_head = 0.07

write(*,*) 'make_probvectors'
call make_reldir
if (MODEL == NEUMANN_MODEL) then
    scale = scale_N
elseif (MODEL == MOORE18_MODEL) then
    scale = scale_M18
endif
open(nfvec,file='d:\immunesystem\lattice\motilitypaper\jump_vectors\cmgui\prob.exnode',status='replace')
beta = 1.0
k = 0
do icase = 0,4
    rho = (icase-1)*0.2
    call compute_dirprobs
    write(nfvec,'(a,i1)') 'Group name : probability',icase
    write(nfvec,'(a)') ' #Fields=1'
    write(nfvec,'(a)') ' 1) vector, field, rectangular cartesian, #Components=3'
    write(nfvec,'(a)') ' x.  Value index= 1, #Derivatives= 0'
    write(nfvec,'(a)') ' y.  Value index= 2, #Derivatives= 0'
    write(nfvec,'(a)') ' z.  Value index= 3, #Derivatives= 0'
    if (icase > 0) then
        do i = 1,nreldir
            k = k+1
            dir = reldir(2,i)
            jump = jumpvec(:,dir)
            vec = jump*dirprob(i)
            call rotate(vec,xang,yang,zang)
            vec = scale*vec
            veclen = norm(vec)
            if (veclen > arrow_head) then
                vec = vec*(veclen-arrow_head)/veclen
            else
                vec = 0
            endif
            write(nfvec,'(a,i3)') 'Node: ',k
            write(nfvec,'(3f8.4)') vec
        enddo
    else
        k = k+1
        vec = (/0.1,0.0,0.0/)
        call rotate(vec,xang,yang,zang)
        vec = scale*vec
        veclen = norm(vec)
        if (veclen > arrow_head) then
            vec = vec*(veclen-arrow_head)/veclen
        else
            vec = 0
        endif
         write(nfvec,'(a,i3)') 'Node: ',k
        write(nfvec,'(3f8.4)') vec
    endif
enddo
close(nfvec)
end subroutine

!---------------------------------------------------------------------
! Rotate the vector about the X-axis by xang, then about the Y-axis by yang
!---------------------------------------------------------------------
subroutine rotate(vec,xang,yang,zang)
real :: vec(3),xang,yang,zang
real :: new(3), cosa, sina

cosa = cos(xang)
sina = sin(xang)
new(1) = vec(1)
new(2) = cosa*vec(2) - sina*vec(3)
new(3) = sina*vec(2) + cosa*vec(3)
vec = new
cosa = cos(zang)
sina = sin(zang)
new(3) = vec(3)
new(1) = cosa*vec(1) - sina*vec(2)
new(2) = sina*vec(1) + cosa*vec(2)
vec = new
cosa = cos(yang)
sina = sin(yang)
new(2) = vec(2)
new(1) = cosa*vec(1) - sina*vec(3)
new(3) = sina*vec(1) + cosa*vec(3)
vec = new

end subroutine

!---------------------------------------------------------------------
! Choose one of the 6 axis directions to allocate to the direction of
! previous step given by jump, which takes values from 0 - 27.
! This could be already on an axis, or could be a D2 or D3 diagonal,
! in which case the choice is random.
!---------------------------------------------------------------------
integer function fix_lastdir(jump,kpar)
integer :: jump, kpar
integer :: k,nax
!                            0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27
integer :: naxes(0:27) = (/  0, 3, 2, 3, 2, 1, 2, 3, 2, 3, 2, 1, 2, 1, 0, 1, 2, 1, 2, 3, 2, 3, 2, 1, 2, 3, 2, 3 /)
integer :: axes(3,27) = reshape( (/ &
!   1      2      3      4      5      6      7      8      9     10     11     12     13     14
1,3,5, 1,3,0, 1,3,6, 1,5,0, 1,0,0, 1,6,0, 1,4,5, 1,4,0, 1,4,6, 3,5,0, 3,0,0, 3,6,0, 5,0,0, 0,0,0, &
!  15     16     17     18     19     20     21     22     23     24     25     26     27
6,0,0, 4,5,0, 4,0,0, 4,6,0, 2,3,5, 2,3,0, 2,3,6, 2,5,0, 2,0,0, 2,6,0, 2,4,5, 2,4,0, 2,4,6 /), (/3,27/) )

if (diagonal_jumps) then
	nax = naxes(jump)
!	if (dbug) write(nfres,*) 'fix_lastdir: jump,nax: ',nax,jump
	if (nax == 0) then
	    write(logmsg,*) 'Should not get here: fix_lastdir: nax=0'
		call logger(logmsg)
	    stop
		fix_lastdir = random_int(1,6,kpar)
	elseif (nax == 1) then
		fix_lastdir = axes(1,jump)
	else
		k = random_int(1,nax,kpar)
!		if (dbug) write(nfres,*) 'random_int: ',k
		fix_lastdir = axes(k,jump)
	endif
else
	write(logmsg,*) 'fix_lastdir: Not MOORE model: ',jump
	call logger(logmsg)
	stop
endif
end function

end module
