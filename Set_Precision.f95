module Set_Precision

Implicit none

Integer,Parameter::Single=Selected_Real_Kind(p=6,r=37)
Integer,Parameter::Double=Selected_Real_Kind(p=13,r=200)
Integer,Parameter::Extended=Selected_Real_Kind(p=17,r=2000)
Integer,Parameter::Quad=Selected_Real_Kind(p=26,r=200)
Integer,Parameter:: prec=Double

end module Set_Precision

