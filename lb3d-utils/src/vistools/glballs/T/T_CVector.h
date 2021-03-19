#include "config.h"

#if !defined(T_CVECTOR_H)
#define T_CVECTOR_H

#include "T_cmath.h"

#include "T.h"

//! 3D Vektor
/*!
 Die Klasse definiert einen Vektor im 3D Raum. Die drei Koordinaten
 sind mit x, y und z bezeichnet.
*/
class T_CVector {
 public:
  //! Defaultkonstruktor
  inline T_CVector(void);
  //! Konstruktor
  inline T_CVector(T_R3_t x,T_R3_t y,T_R3_t z);
  //! setzt den Vector auf den Ursprung
  inline void setZero(void);
  //! pr"uft ob der Vector auf dem Ursprung liegt
  inline bool isZero(void);
  //! setzt den Vektor
  inline void set(const T_R3_t x,const T_R3_t y,const T_R3_t z);
  //! setzt den Vektor
  inline void set(const T_R3_t *p);
  //! setzt die X Koordinate
  inline void setX(const T_R3_t x);
  //! setzt die Y Koordinate
  inline void setY(const T_R3_t y);
  //! setzt die Z Koordinate
  inline void setZ(const T_R3_t z);
  //! addiert Wert zur X Koordinate
  inline void addX(const T_R3_t x);
  //! addiert Wert zur Y Koordinate
  inline void addY(const T_R3_t y);
  //! addiert Wert zur Z Koordinate
  inline void addZ(const T_R3_t z);
  //! liefert alle Koordinaten
  inline void get(T_R3_t *px,T_R3_t *py,T_R3_t *pz);
  //! liefert alle Koordinaten
  inline void get(T_R3_t *p);
  //! liefert die X Koordinate
  inline T_R3_t getX(void) const;
  //! liefert die Y Koordinate
  inline T_R3_t getY(void) const;
  //! liefert die Z Koordinate
  inline T_R3_t getZ(void) const;
  //! addiert zwei Vektoren
  inline T_CVector operator+(const T_CVector &CV) const;
  //! subtrahiert Vektoren von einander
  inline T_CVector operator-(const T_CVector &CV) const;
  //! multipliziert einen skalaren Wert auf eine Vektor
  inline T_CVector operator*(const T_R3_t scale) const;
  //! "uber nimmt den Wert eines Vektors
  inline T_CVector &operator=(const T_CVector &CV);
  //! Vergleich zwei Vektoren
  inline bool operator==(const T_CVector &CV) const;
  //! Vergleich zwei Vektoren
  inline bool operator!=(const T_CVector &CV) const;
  //! addiert einen Vektor
  inline T_CVector &operator+=(const T_CVector &CV);
  //! subtrahiert einen Vektor
  inline T_CVector &operator-=(const T_CVector &CV);
  //! berechnet die L"ange eines Vektors
  inline T_R3_t length(void) const;
  //! berechnet des Quadrat der L"ange eines Vektors
  inline T_R3_t lengthSquare(void) const;
  //! reskaliert einen Vektor auf die L"ange eins und gibt die urspr"ungliche L"ange zur"uck
  inline T_R3_t normalize(void);
  //! multipliert einen skalaren Wert auf den Vektor
  inline T_CVector &scale(T_R3_t length);
  //! Spieglung des Vektors am Ursprung
  inline T_CVector &invert(void);
  //! gibt eine normalisierte Version des Vektors zur"uck
  inline T_CVector norm(void) const;
  //! gibt das Kreuzprodukt zweier Vektoren zur"uck
  inline T_CVector cross(const T_CVector &pCV) const;
  //! berechnet das Vektorprodukt zweier Vektoren
  inline T_R3_t scalar(const T_CVector &pCV);
  //! berechnet das Vektorprodukt zweier Vektoren
  inline T_R3_t scalar(const T_R3_t x,const T_R3_t y,const T_R3_t z);
  //! liefert die Projektion des Vektors auf den "ubergebenen zur"uck
  inline T_CVector projection(const T_CVector &pCV);
  //
  inline T_R3_t split(const T_CVector &CV,T_CVector &CVNormalPart,T_CVector &CVTangentialPart);
  void dg(void){fprintf(stderr,"<%f, %f, %f>\n",m_x,m_y,m_z);}
 private:
  //! die X-Koordinate des Vektors
  T_R3_t m_x;
  //! die Y-Koordinate des Vektors
  T_R3_t m_y;
  //! die Z-Koordinate des Vektors
  T_R3_t m_z;
};

#include "T_CVector_inline.h"

#endif
