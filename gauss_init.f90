subroutine gauss_init()
use variables
print*, "Calculating Gauss Values!"
wg(:,:)=0.0_prec
xg(:,:)=0.0_prec

wg(1,1)=2.0_prec

wg(2,1)=1.0_prec
wg(2,2)=1.0_prec

wg(3,1)=5.0_prec/9.0_prec
wg(3,2)=8.0_prec/9.0_prec
wg(3,3)=5.0_prec/9.0_prec

wg(4,1)=(18.0_prec - sqrt(30.0_prec)) / 36.0_prec
wg(4,2)=(18.0_prec + sqrt(30.0_prec)) / 36.0_prec
wg(4,3)=(18.0_prec + sqrt(30.0_prec)) / 36.0_prec
wg(4,4)=(18.0_prec - sqrt(30.0_prec)) / 36.0_prec

wg(5,1)=(322.0_prec - 13.0_prec * sqrt(70.0_prec)) / 900.0_prec
wg(5,2)=(322.0_prec + 13.0_prec * sqrt(70.0_prec)) / 900.0_prec
wg(5,3)=128.0_prec/225.0_prec
wg(5,4)=(322.0_prec + 13.0_prec * sqrt(70.0_prec)) / 900.0_prec
wg(5,5)=(322.0_prec - 13.0_prec * sqrt(70.0_prec)) / 900.0_prec

xg(1,1)=0.0_prec

xg(2,1)=-1.0_prec/sqrt(3.0_prec)
xg(2,2)=1.0_prec/sqrt(3.0_prec)

xg(3,1)= -sqrt(15.0_prec)/5.0_prec
xg(3,2)=0.0_prec
xg(3,3)= sqrt(15.0_prec)/5.0_prec

xg(4,1)=-sqrt(525.0_prec+70.0_prec*sqrt(30.0))/35.0_prec
xg(4,2)=-sqrt(525.0_prec-70.0_prec*sqrt(30.0))/35.0_prec
xg(4,3)=sqrt(525.0_prec-70.0_prec*sqrt(30.0))/35.0_prec
xg(4,4)=sqrt(525.0_prec+70.0_prec*sqrt(30.0))/35.0_prec

xg(5,1)= -sqrt(245.0_prec + 14.0_prec * sqrt(70.0_prec)) / 21.0_prec
xg(5,2)= -sqrt(245.0_prec - 14.0_prec * sqrt(70.0_prec)) / 21.0_prec
xg(5,3)=0.0_prec
xg(5,4)= sqrt(245.0_prec - 14.0_prec * sqrt(70.0_prec)) / 21.0_prec
xg(5,5)= sqrt(245.0_prec + 14.0_prec * sqrt(70.0_prec)) / 21.0_prec

end subroutine gauss_init
