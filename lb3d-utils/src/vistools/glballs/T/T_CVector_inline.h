// Defaultkonstruktor
/*
*/
T_CVector::T_CVector(void)
{
 setZero();
}

// Konstruktor
/*
*/
T_CVector::T_CVector(T_R3_t x,T_R3_t y,T_R3_t z)
{
 m_x=x;
 m_y=y;
 m_z=z;
}

// setzt den Vector auf den Ursprung
/*!
 Ist der schnellste Weg alle Koordinaten auf den Wert 0. zu setzen.
*/
void T_CVector::setZero(void)
{
 m_x=T_R3_ZERO;
 m_y=T_R3_ZERO;
 m_z=T_R3_ZERO;
}

// pr"uft ob der Vector auf dem Ursprung liegt
/*
*/
bool T_CVector::isZero(void)
{
 return (m_x==T_R3_ZERO) && (m_y==T_R3_ZERO) && (m_z==T_R3_ZERO);
}

// setzt den Vektor
/*
*/
void T_CVector::set(const T_R3_t x,const T_R3_t y,const T_R3_t z)
{
 m_x=x;
 m_y=y;
 m_z=z;
}

// setzt den Vektor
/*
*/
void T_CVector::set(const T_R3_t *p)
{
 m_x=p[0];
 m_y=p[1];
 m_z=p[2];
}

// setzt die X Koordinate
/*
*/
void T_CVector::setX(const T_R3_t x)
{
 m_x=x;
}

// setzt die Y Koordinate
/*
*/
void T_CVector::setY(const T_R3_t y)
{
 m_y=y;
}

// setzt die Z Koordinate
/*
*/
void T_CVector::setZ(const T_R3_t z)
{
 m_z=z;
}

// addiert Wert zur X Koordinate
/*
*/
void T_CVector::addX(const T_R3_t x)
{
 m_x+=x;
}

// addiert Wert zur Y Koordinate
/*
*/
void T_CVector::addY(const T_R3_t y)
{
 m_y+=y;
}

// addiert Wert zur Z Koordinate
/*
*/
void T_CVector::addZ(const T_R3_t z)
{
 m_z+=z;
}

// liefert alle Koordinaten
/*
*/
void T_CVector::get(T_R3_t *px,T_R3_t *py,T_R3_t *pz)
{
 *px=m_x;
 *py=m_y;
 *pz=m_z;
}

// liefert alle Koordinaten
/*
*/
void T_CVector::get(T_R3_t *p)
{
 p[0]=m_x;
 p[1]=m_y;
 p[2]=m_z;
}

// liefert die X Koordinate
/*
*/
T_R3_t T_CVector::getX(void) const
{
 return m_x;
}

// liefert die Y Koordinate
/*
*/
T_R3_t T_CVector::getY(void) const
{
 return m_y;
}

// liefert die Z Koordinate
/*
*/
T_R3_t T_CVector::getZ(void) const
{
 return m_z;
}

// addiert zwei Vektoren
/*!
 Es werden zwei Vektoren addiert und das Ergebnis als neuer Vektor
 zur"uckgegeben
*/
T_CVector T_CVector::operator+(const T_CVector &CV) const
{
 return T_CVector(
  m_x+CV.getX(),
  m_y+CV.getY(),
  m_z+CV.getZ()
 );
}

// subtrahiert Vektoren von einander
/*!
 subtrahiert den zweiten Vektor vom ersten und gibt das Ergebnis als
 neuen Vektor zur"uck.
*/
T_CVector T_CVector::operator-(const T_CVector &CV) const
{
 return T_CVector(
  m_x-CV.getX(),
  m_y-CV.getY(),
  m_z-CV.getZ()
 );
}

// multipliziert einen skalaren Wert auf eine Vektor
/*!
 Das Ergebnis wird als neuer Vektor zur"uckgegeben.
*/
T_CVector T_CVector::operator*(const T_R3_t scale) const
{
 return T_CVector(
  scale*m_x,
  scale*m_y,
  scale*m_z
 );
}

// "uber nimmt den Wert eines Vektors
/*
*/
T_CVector &T_CVector::operator=(const T_CVector &CV)
{
 m_x=CV.getX();
 m_y=CV.getY();
 m_z=CV.getZ();
 return *this;
}

// Vergleich zwei Vektoren
/*
*/
bool T_CVector::operator==(const T_CVector &CV) const
{
 return (m_x==CV.getX()) && (m_y==CV.getY()) && (m_z==CV.getZ());
}

// Vergleich zwei Vektoren
/*
*/
bool T_CVector::operator!=(const T_CVector &CV) const
{
 return (m_x!=CV.getX()) || (m_y!=CV.getY()) || (m_z!=CV.getZ());
}

// addiert einen Vektor
/*
*/
T_CVector &T_CVector::operator+=(const T_CVector &CV)
{
 m_x+=CV.getX();
 m_y+=CV.getY();
 m_z+=CV.getZ();
 return *this;
}

// subtrahiert einen Vektor
/*
*/
T_CVector &T_CVector::operator-=(const T_CVector &CV)
{
 m_x-=CV.getX();
 m_y-=CV.getY();
 m_z-=CV.getZ();
 return *this;
}

// berechnet die L"ange eines Vektors
/*!
 Die Berechnung erfolgt gem"a"s:
 \f[ l=\sqrt{x^2+y^2+z^2} \f]
*/
T_R3_t T_CVector::length(void) const
{
 return std::sqrt(m_x*m_x+m_y*m_y+m_z*m_z);
}

// berechnet des Quadrat der L"ange eines Vektors
/*!
 Ist auf den zweiten Blick klar, braucht man z.B. das Quadrat eines
 Abstandes, ist es schneller wenn man sich das Wurzelziehen spart.
 folglich gibt es hierf"ur eine zu bevorzugende Memberfunktion,
 definiert gem"a"s:
 \f[ l=x^2+y^2+z^2 \f]
*/
T_R3_t T_CVector::lengthSquare(void) const
{
 return m_x*m_x+m_y*m_y+m_z*m_z;
}

// reskaliert einen Vektor auf die L"ange eins und gibt die urspr"ungliche L"ange zur"uck
/*!
 Die Funktion f"uhrt eins zu eins folgendes aus:
 \f[ \vec{v}:=\vec{v}/{|\vec{v}|}
*/
T_R3_t T_CVector::normalize(void)
{
 T_R3_t l=length();
 T_R3_t scale=1./l;
 m_x*=scale; 
 m_y*=scale; 
 m_z*=scale; 
 return l;
}

// multipliert einen skalaren Wert auf den Vektor
/*!
 Im Gegensatz zum Operator wird der Vektor selbst skaliert.
 Es wird keine neue Instanz der Klasse T_CVektor generiert
*/
T_CVector &T_CVector::scale(T_R3_t scale)
{
 m_x*=scale; 
 m_y*=scale; 
 m_z*=scale; 
 return *this;
}

// Spieglung des Vektors am Ursprung
/*!
 Entspricht der Vorzeichen umkehr aller Koordinaten.
*/
T_CVector &T_CVector::invert(void)
{
 m_x=-m_x;
 m_y=-m_y;
 m_z=-m_z;
 return *this;
}

// gibt eine normalisierte Version des Vektors zur"uck
/*
*/
T_CVector T_CVector::norm(void) const
{
 T_R3_t scale=1./length();
 return T_CVector(m_x*scale,m_y*scale,m_z*scale); 
}

// gibt das Kreuzprodukt zweier Vektoren zur"uck
/*!
 In der Hoffnung, da"s das Kreuzprodukt so definiert ist:
 \f[\vec{u}=\vec{v}\times\vec{w}\f]
 mit
 \f[u_x=v_yw_z-v_zw_y\f]
 \f[u_y=v_zw_x-v_xw_z\f]
 \f[u_z=v_xw_y-v_yw_x\f]
 im Source sieht das ganze dann so aus:
 
 u=v.cross(w);
*/
T_CVector T_CVector::cross(const T_CVector &CV) const
{
 return T_CVector(
  (m_y*CV.getZ())-(m_z*CV.getY()),
  (m_z*CV.getX())-(m_x*CV.getZ()),
  (m_x*CV.getY())-(m_y*CV.getX())
 );
}

// berechnet das Vektorprodukt zweier Vektoren
/*!
 dies entspricht:
 \f[l=\vec{v}\vec{w}\f]
 im Source sieht das ganze dann so aus:
 
 l=v.scalar(w);
*/
T_R3_t T_CVector::scalar(const T_CVector &CV)
{
 return ((m_x*CV.getX())+(m_y*CV.getY()))+(m_z*CV.getZ());
}

// berechnet das Vektorprodukt zweier Vektoren
/*
*/
T_R3_t T_CVector::scalar(const T_R3_t x,const T_R3_t y,const T_R3_t z)
{
 return m_x*x+m_y*y+m_z*z;
}

// liefert die Projektion des Vektors auf den "ubergebenen zur"uck
/*!
 mathematisch :
 \f[ \vec{u}=(\vec{v}\vec{n})\vec{n} \f]
 im Source sieht das ganze dann so aus:
 
 u=v.projection(n);
*/
T_CVector T_CVector::projection(const T_CVector &CV)
{
 T_R3_t scale=scalar(CV); 
 return T_CVector(
  scale*CV.getX(),
  scale*CV.getY(),
  scale*CV.getZ()
 );
}

T_R3_t T_CVector::split(const T_CVector &CV,T_CVector &CVNormalPart,T_CVector &CVTangentialPart)
{
 T_R3_t NormalPartLength=scalar(CV);
 CVNormalPart.set(NormalPartLength*this->m_x,NormalPartLength*this->m_y,NormalPartLength*this->m_z);
 CVTangentialPart=CV-CVNormalPart;
 return NormalPartLength;
}
