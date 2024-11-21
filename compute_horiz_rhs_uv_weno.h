#if defined UV_COR || (defined CURVGRID && defined UV_ADV)
        do j=jstrV-1,jend                ! Add Coriolis terms and
          do i=istrU-1,iend              ! contribution to advection
            cff=0.5*Hz(i,j,k)*(          ! associated with curvilinear
# ifdef UV_COR
     &              fomn(i,j)            ! horizontal coordinates.
# endif
# if (defined CURVGRID && defined UV_ADV)
     &             +0.5*( (v(i,j,k,nrhs)+v(i,j+1,k,nrhs))*dndx(i,j)
     &                   -(u(i,j,k,nrhs)+u(i+1,j,k,nrhs))*dmde(i,j))
# endif
     &                                                             )
# ifdef WEC
#  if (defined CURVGRID && defined UV_ADV)
     ! Define a cff1 that is added to UFx, VFe for stokes terms 
           cff1=0.5*Hz(i,j,k)*(
     &          0.5*( dndx(i,j)*(vst(i,j,k)+vst(i,j+1,k))
     &               -dmde(i,j)*(ust(i,j,k)+ust(i+1,j,k)) ))
#  else
           cff1 = 0.0
#  endif
           UFx(i,j)=(cff+cff1)*(v(i,j,k,nrhs)+v(i,j+1,k,nrhs))
           UFe(i,j)=cff*(vst(i,j,k)+vst(i,j+1,k))
           VFe(i,j)=(cff+cff1)*(u(i,j,k,nrhs)+u(i+1,j,k,nrhs))
           VFx(i,j)=cff*(ust(i,j,k)+ust(i+1,j,k)) 
# else /* no WEC */

            UFx(i,j)=cff*(v(i,j,k,nrhs)+v(i,j+1,k,nrhs))
            VFe(i,j)=cff*(u(i,j,k,nrhs)+u(i+1,j,k,nrhs))
# endif
          enddo
        enddo
        do j=jstr,jend
          do i=istrU,iend
            ru(i,j,k)=ru(i,j,k)+0.5*(UFx(i,j)+UFx(i-1,j))
# ifdef WEC
     &               + 0.5*(UFe(i,j)+UFe(i-1,j))
# endif
          enddo
        enddo
        do j=jstrV,jend
          do i=istr,iend
            rv(i,j,k)=rv(i,j,k)-0.5*(VFe(i,j)+VFe(i,j-1))
# ifdef WEC
     &               -0.5*(VFx(i,j)+VFx(i,j-1))
# endif
          enddo
        enddo

# ifdef DIAGNOSTICS
	! JM It would be better to not do the u,v boundary points
        if (CORR_STAGE) then
          if (diag_uv.and.calc_diag) then
            Udiag(:,:,k,icori) = 0.5*(UFx(1:nx,1:ny)+UFx(0:nx-1,1:ny))*dxdyi_u
            Vdiag(:,:,k,icori) =-0.5*(VFe(1:nx,1:ny)+VFe(1:nx,0:ny-1))*dxdyi_v
          endif
        endif
# endif /* DIAGNOSTICS */

#endif  /* UV_COR */


#ifdef UV_ADV

! Add horizontal advection of momentum: compute diagonal [UFx,VFe]
! and off-diagonal [UFe,VFx] components of momentum flux tensor due
! to horizontal advection; after that add their divergence to r.h.s.

        call fluxes_uv(UFx,VFx,UFe,VFe,FlxU(:,:,k),FlxV(:,:,k),u(:,:,k,nrhs),v(:,:,k,nrhs))

        do j=jstr,jend
          do i=istrU,iend
            ru(i,j,k)=ru(i,j,k)-UFx(i+1,j) + UFx(i,j)
     &                         -UFe(i,j+1) + UFe(i,j)
          enddo
        enddo
        do j=jstrV,jend
          do i=istr,iend
            rv(i,j,k)=rv(i,j,k)-VFx(i+1,j)+VFx(i,j)
     &                         -VFe(i,j+1)+VFe(i,j)
          enddo
        enddo
#endif /* UV_ADV */
