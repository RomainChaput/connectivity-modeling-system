!****************************************************************************
!* Extension for Connectivity Modeling System (CMS)                               *
!* File : mod_orientation.f90                                                          *
!* Created: 2020-03-25                                                *
!* Code contributors: Romain Chaput, Mohamed Iskandarani                                 *
!*                                                                          *
!****************************************************************************
! Give orientation abilities to particles in the CMS. Follows biased correlated random walk as described by Codling et al., 2004 - MEPS and Staaterman et al., 2012 - JTB
! Published in Chaput R, Sochala P, Miron P, Kourafalou V, Iskandarani M, 2022. Quantitative uncertainty estimation in biophysical models of fish larval connectivity in the Florida Keys. ICES Journal of Marine Science.

MODULE mod_orientation

USE globalvariables
USE mod_kinds
USE constants
USE mod_random
USE mod_iounits
USE mod_calendar

IMPLICIT NONE

real (kind = real_kind), allocatable, save, private :: polyCenterLons(:), &
polyCenterLats(:), polyCenterId(:), polyCenterKappa(:)

integer (kind = int_kind), save, private :: nbPoly ! length polygon file

CONTAINS

!   ==================== Load polygon files ================

! Subroutine to load the file containing the information about the centers of the polygons. Use the example of the subroutine load_reef_data in mod_reef.f90 of the open source CMS.
!Read and allocate latitudes, longitudes, and Kappa
subroutine load_center_data(centerFileName)

character (len = *), intent(in) :: centerFileName
integer (kind = int_kind) :: i, iunit
real (kind = real_kind) :: dumdum

call getSize(centerFileName,nbPoly)
! allocate names
allocate(polyCenterLons(nbPoly))
!print *,'polyCenterLons = ', size(polyCenterLons) ! for debug
allocate(polyCenterLats(nbPoly))
allocate(polyCenterId(nbPoly))
allocate(polyCenterKappa(nbPoly))

call get_unit(iunit)
open(UNIT=iunit,FILE=trim(centerFileName), STATUS='old')

! read file with centers of polygons
do i=1,nbPoly
	read(iunit,*) polyCenterLons(i), polyCenterlats(i), polyCenterId(i), polyCenterKappa(i)
end do

! close file
call release_unit(iunit)

!make sure all longitudes of the polygons are between 0 and 360: from mod_reef.f90
  do i=1,nbPoly
    do while (polyCenterLons(i) < 0.) 
      polyCenterLons(i) = polyCenterLons(i) + 360.
    end do 
    do while (polyCenterLons(i) > 360.) 
      polyCenterLons(i) = polyCenterLons(i) - 360.
    end do
  end do
  
 call random_von_Mises(polyCenterKappa(1), .TRUE., dumdum) ! Initializing the Von Mises sampler for faster computation. Only works if all the reefs have the same Kappa

end subroutine load_center_data

!  ===================== Find nearest reef ===============

! Subroutine to compute the distance between the larva and the reef to find if the orientation module can be called.
subroutine nearest_reef(lon_now, lat_now, lon_reef, lat_reef, Kappa_reef, reef_dist)

real (kind = real_kind), intent(in) :: lon_now, lat_now
real (kind = real_kind), intent(out) :: lon_reef, lat_reef, Kappa_reef
real (kind = real_kind) :: dist, reef_dist
real, dimension(nbPoly) :: distance
integer (kind = int_kind) :: i, nearest_reef_ID

do i=1,nbPoly
	lon_reef = polyCenterLons(i)
	lat_reef = polyCenterLats(i)
	call Distance_Sphere(lon_now, lat_now, lon_reef, lat_reef, dist)
	distance(i) = dist
end do

! Find smallest distance and nearest reef
reef_dist = MINVAL(distance)
nearest_reef_ID = MINLOC(distance,1)
lon_reef = polyCenterLons(nearest_reef_ID)
lat_reef = polyCenterLats(nearest_reef_ID)
Kappa_reef = polyCenterKappa(nearest_reef_ID)

end subroutine nearest_reef

! ================= Compute orientation U and V velocity =====

! Subroutine to compute uorient and vorient when the larvae are orienting toward the nearest reef. Compute theta and swimming speed.
! Based on biased correlated random walk algorithm developed in Codling et al., 2004 and orientation behavior implemented in Staaterman et al., 2012.
subroutine calc_orient(lon_now, lat_now, lon_old, lat_old, lon_reef, lat_reef, Kappa_reef, reef_dist, run_time, dtturb, vorient, uorient)

real (kind = real_kind), intent(in) :: lon_now, lat_now, lon_old, lat_old, lon_reef, lat_reef, Kappa_reef, reef_dist, dtturb
integer (kind = int8_kind), intent(in) :: run_time
real (kind = real_kind), intent(out) :: vorient, uorient
real (kind = real_kind) :: X, Y, thetaCurrent, thetaPref, mu, theta, d, ti, age, PLD, swimmingSpeed, htvelscl, normrn

if (reef_dist > maxDistance) then	
	htvelscl = sqrt((2*horDiffOrient)/dtturb)
	call random_gaussian(0.,1.,normrn)
	uorient = normrn*htvelscl
	call random_gaussian(0.,1.,normrn)
	vorient = normrn*htvelscl
	
else
	! Spatial dependence of orientation
	d = 1 - (reef_dist/maxDistance)
	! Computing direction of the reef
	thetaPref = -Haversine(lon_now, lat_now, lon_reef, lat_reef)
	! Computing direction larva is heading from previous point
	thetaCurrent = Haversine(lon_old, lat_old, lon_now, lat_now)
	! Turning angle
	mu = -d*(thetaCurrent - thetaPref)
	
	! Computing the Von Mises distribution to pick new direction
	! Biased Correlated Random Walk: Colding et al., 2004 MEPS
	! Random sampling of Von Mises
	call random_von_Mises(Kappa_reef, .FALSE., ti)
	! Apply correction
	theta = ti - thetaCurrent - mu
	
	! Compute swimming speed for the larvae following Staaterman et al., 2012 JTB and def_globalvariables.f90
	age = run_time/real(secs_in_day)
	PLD = timeMax/real(secs_in_day)
	swimmingSpeed = swimmingSpeedHatch + 10**((log10(age)/log10(PLD))*log10(swimmingSpeedSettle - swimmingSpeedHatch))
	swimmingSpeed = swimmingSpeed/100
	
	! Compute u and v orientation velocity
	uorient = swimmingSpeed * cos(theta)
	vorient = swimmingSpeed * sin(theta)
	
end if
end subroutine calc_orient

! Subroutine to compute uorient and vorient when the larvae are orienting against the currents.
subroutine calc_rheotaxis(vf, uf, run_time, vorient, uorient)

	real (kind = real_kind), intent(in) :: vf, uf
	integer (kind = int8_kind), intent(in) :: run_time
	real (kind = real_kind), intent(out) :: vorient, uorient
	real (kind = real_kind) :: age, PLD, swimmingSpeed, ti, thetaRheo, theta, uv, norm_swim

	! Compute swimming speed for the larvae following Staaterman et al., 2012 JTB
	age = run_time/real(secs_in_day)
	PLD = timeMax/real(secs_in_day)
	swimmingSpeed = swimmingSpeedHatch + 10**((log10(age)/log10(PLD))*log10(swimmingSpeedSettle - swimmingSpeedHatch))
	swimmingSpeed = swimmingSpeed/100
	norm_swim = abs(swimmingSpeed)
	
	!Compute current speed absolute value
	uv = sqrt(uf*uf+vf*vf)
	
	! Compute randomness of direction
	call random_von_Mises(5.0 , .TRUE., ti)
	! Compute rheotaxis heading
	thetaRheo = -atan2(vf, uf)
	theta = thetaRheo + ti

		if (norm_swim < uv) then
		! Compute u and v rheotaxis velocity
		uorient = swimmingSpeed * cos(theta)
		vorient = swimmingSpeed * sin(theta)
		else
		uorient = uv * cos(theta)
		vorient = uv * sin(theta)
		end if

end subroutine calc_rheotaxis

! Subroutine to compute uorient and vorient when the larvae have a fixed cardinal orientation.
subroutine calc_cardinal(Cardinal_heading, run_time, vorient, uorient)

real (kind = real_kind), intent(in) :: Cardinal_heading
integer (kind = int8_kind), intent(in) :: run_time
real (kind = real_kind), intent(out) :: vorient, uorient
real (kind = real_kind) :: thetaCard, theta, ti, age, PLD, swimmingSpeed
	
	! Computing preferred direction
	thetaCard = deg2rad*Cardinal_heading
	! Random sampling of Von Mises
	call random_von_Mises(5.0, .TRUE., ti)
	! Preferred direction plus stochastic behavior
	theta = thetaCard + ti
	
	! Compute swimming speed for the larvae following Staaterman et al., 2012 JTB
	age = run_time/real(secs_in_day)
	PLD = timeMax/real(secs_in_day)
	swimmingSpeed = swimmingSpeedHatch + 10**((log10(age)/log10(PLD))*log10(swimmingSpeedSettle - swimmingSpeedHatch))
	swimmingSpeed = swimmingSpeed/100
	
	! Compute u and v orientation velocity
	uorient = swimmingSpeed * cos(theta)
	vorient = swimmingSpeed * sin(theta)
	
end subroutine calc_cardinal

!=====================================================================

! Haversine formula to compute angles from lat and lon
real (kind = real_kind) function Haversine(lon1,lat1,lon2,lat2)
	real (kind = real_kind), intent(in) :: lon1, lat1, lon2, lat2
	real (kind = real_kind) :: X, Y, rlon1, rlon2, rlat1, rlat2
	
	rlon1 = deg2rad*lon1
	rlat1 = deg2rad*lat1
	rlon2 = deg2rad*lon2
	rlat2 = deg2rad*lat2
	X = cos(rlat2)*sin(rlon2 - rlon1)
    Y = cos(rlat1)*sin(rlat2) - sin(rlat1)*cos(rlat2)*cos(rlon2 - rlon1)
    Haversine = atan2(Y,X)
end function Haversine

!================================================
! deallocate memory of polygon files
! use the example of the subroutine dealloc_reef in mod_reef.f90 of the open source CMS
subroutine dealloc_orient

 deallocate(polyCenterLons)
 deallocate(polyCenterLats)
 deallocate(polyCenterId)
 deallocate(polyCenterKappa)

end subroutine dealloc_orient

!================================================
END MODULE mod_orientation
