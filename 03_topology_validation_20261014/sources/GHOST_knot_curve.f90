  subroutine knot_curve(icurve,r0,s,nc,pos,tan)
    use var
    implicit none
    integer      , intent (in) :: icurve
    real(kind=GP), intent (in) :: r0,s
    integer      , intent(out) :: nc
    real(kind=GP), intent(out) :: pos(3,2),tan(3,2)
    select case (icurve)
    case (1) ! ring in the plane z = pi
      nc = 1
      pos(:,1) = (/ pi+r0*cos(s), pi+r0*sin(s), pi /)
      tan(:,1) = (/ -r0*sin(s), r0*cos(s), 0.0_GP /)
    case (2) ! trefoil
      nc = 1
      pos(:,1) = (/ pi+(r0*(sin(s)+2.0_GP*sin(2.0_GP*s)))/3.0_GP,   &
                    pi+(r0*(cos(s)-2.0_GP*cos(2.0_GP*s)))/3.0_GP,   &
                    pi-(r0*sin(3.0_GP*s))/3.0_GP /)
      tan(:,1) = (/ (r0*(cos(s)+4.0_GP*cos(2.0_GP*s)))/3.0_GP,      &
                    (r0*(-sin(s)+4.0_GP*sin(2.0_GP*s)))/3.0_GP,     &
                    -(r0*cos(3.0_GP*s)) /)
    case default ! two linked rings
      nc = 2
      pos(:,1) = (/ pi-r0/2+r0*cos(s), pi+r0*sin(s), pi /)
      tan(:,1) = (/ -r0*sin(s), r0*cos(s), 0.0_GP /)
      pos(:,2) = (/ pi+r0/2+r0*cos(s), pi, pi+r0*sin(s) /)
      tan(:,2) = (/ -r0*sin(s), 0.0_GP, r0*cos(s) /)
    end select
  end subroutine knot_curve
  
end module ic_quantum
