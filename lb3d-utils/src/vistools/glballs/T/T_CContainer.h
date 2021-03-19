#include "config.h"

#if !defined(T_CCONTAINER_H)
#define T_CCONTAINER_H

#include "T.h"

/*!
 F"uhrt eine einfache for Schleife mit dem Iterator I "uber die verkettete
 Liste in CC aus. Dabei wird jedes Kettenelement als Objekt des Types
 TYPE verstanden und deren Memberfunktion FUNC aufgerufen.
*/
#define T_CC_FOR_EACH_ITEM_CALL(I,CC,TYPE,FUNC) \
for(I.begin(CC);I.isValid();++I){\
 ((TYPE *)I.get())->FUNC;\
}


/*!
 F"uhrt eine einfache for Schleife mit dem Iterator I "uber die verkettete
 Liste in CC aus. Dabei wird jedes Kettenelement als Objekt des Types
 TYPE verstanden und deren Memberfunktion FUNC aufgerufen.
 Als R"uckgabewert wird TRUE erwartet, wenn die Function erfolgreich
 ausgef"uhrt wurde. Wird FALSE zur"uckgegeben, wird die Schleife abgebrochen
 und FLAG wird auf FALSE gesetzt. Die Schleife wird nicht ausgef"uhrt, wenn 
 FLAG zu Beginn den Wert FALSE hat.
*/
#define T_CC_FOR_EACH_ITEM_BOOLCALL(I,CC,TYPE,FUNC,FLAG) \
for(I.begin(CC);I.isValid() && FLAG;++I){\
 FLAG&=((TYPE *)I.get())->FUNC;\
}


/*!
 F"uhrt eine einfache for Schleife mit dem Iterator I "uber die verkettete
 Liste in CC aus. Dabei wird jedes Kettenelement als Objekt des Types
 TYPE verstanden und der Memberfunktion FUNC "ubergeben.
*/
#define T_CC_FOR_EACH_ITEM_EXEC(I,CC,TYPE,FUNC) \
for(I.begin(CC);I.isValid();++I){\
 FUNC((TYPE *)I.get());\
}

/*!
 F"uhrt eine einfache for Schleife mit dem Iterator I "uber die verkettete
 Liste in CC aus. Dabei wird jedes Kettenelement als Objekt des Types
 TYPE verstanden und der Memberfunktion FUNC "ubergeben.
 Als zweiter Parameter wird T "ubergeben. 
*/
#define T_CC_FOR_EACH_ITEM_EXEC_T(I,CC,TYPE,FUNC,T) \
for(I.begin(CC);I.isValid();++I){\
 FUNC((TYPE *)I.get(),T);\
}



class T_CContainer;
class T_CContainerIterator;
class T_CContainerObject;

/*************************************************************************/

//! Element einer bidirektionalen verketten Liste
/*!
 Das Element ist genau einer Liste zugeordnet, wird das Element einer Liste
 zugeordnet, wird es automatisch aus der ausgetragen, der es zu diesem 
 Zeitpunkt zugeornet ist.
*/
class T_CContainerItem {
 friend class T_CContainer;
 friend class T_CContainerIterator;
 public:
  //! Default Konstruktor
  T_CContainerItem(T_CContainerObject *pObject);
  //! Default Destruktor
  ~T_CContainerItem();
  //! setzt die Identifikaion des Listenelementes
  inline void setId(T_Id_t Id);
  //! liefert die Identifikaion des Listenelementes
  inline T_Id_t getId(void);
  //! entfernt dieses Element aus der verketteten Liste
  void unlink();
  //! f"ugt dieses Element in eine verketteten Liste ein
  void linkWith(T_CContainer *pCC);
  //! meldet, ob das Item einer verketteten Liste zugeordnet wurde
  inline bool hasContainer(void);
 private:
  //! ordnet dem Item einen neuen Container zu
  inline void setContainer(T_CContainer *pCC);
  //! liefert das referenzierte Objekt des Listenelementes
  inline T_CContainerObject *get(void);
  //! setzt den Vorg"anger dieses Elementes im Sinne einer verketten Liste
  inline void setPrev(T_CContainerItem *pCCI);
  //! setzt den Nachfolger dieses Elementes im Sinne einer verketten Liste
  inline void setNext(T_CContainerItem *pCCI);
  //! liefert den Vorg"anger in der verketteten Liste
  inline T_CContainerItem *getPrev(void) const;
  //! liefert den Nachfolger in der verketteten Liste
  inline T_CContainerItem *getNext(void) const;
 private:
  //! das Objekt, da"s von diesem Kettenglied verwaltet wird
  T_CContainerObject *m_pCCObject;
  //! das n"achste Glied der verketteten Liste
  /*!
   gibt es kein nachfolgendes Kettenglied, steht hier NULL
  */
  T_CContainerItem *m_pCCINext;
  //! das vorherige Glied der verketteten Liste
  /*!
   gibt es kein vorangehendes Kettenglied, steht hier NULL
  */
  T_CContainerItem *m_pCCIPrev;
  //! die verkettete Liste in die dieses Kettenglied eingebunden ist
  T_CContainer *m_pCC; 
  //! die Identifikation des Kettengliedes
  T_Id_t m_Id;
};

/*************************************************************************/

//! Basisklasse f"ur verkettbare Objekte
/*!
 Damit Objekte einer Klasse in eine verkettete Liste eingebunden werden
 k"onnen m"ussen die Klasse von T_CContainerObject abgeleitet sein.
 Zus"atzlich mu"s f"ur jede Einbindem"oglichkeit ein Objekt des
 Types T_CContainerItem existieren.
*/
class T_CContainerObject {
 public:
  //! Defaultkonstruktor
  T_CContainerObject(void);
  //! Defaultdekonstruktor
  virtual ~T_CContainerObject(void){};
  //! setzt die Identifikaion des Listenelementes
  inline void setId(T_Id_t Id);
  //! liefert die Identifikaion des Listenelementes
  inline T_Id_t getId(void);
  //! entfernt dieses Objekt aus der verketteten Liste
  inline void unlink();
  //! f"ugt dieses Objekt in eine verketteten Liste ein
  inline void linkWith(T_CContainer *pCC);
  //! meldet, ob das Item einer verketteten Liste zugeordnet wurde
  inline bool hasContainer(void);
 private:
  //! die Einbindung des Objektes erfolgt "uber dieses Kettenglied
  T_CContainerItem m_CCI;
};

/*************************************************************************/

//! Basisklasse f"ur Objekte, die in zwei Listen verkettbar sind

/*!
 Damit Objekte einer Klasse in eine verkettete Liste eingebunden werden
 k"onnen m"ussen die Klasse von T_CContainerObject abgeleitet sein.
 Zus"atzlich mu"s f"ur jede Einbindem"oglichkeit ein Objekt des
 Types T_CContainerItem existieren.

 Es gibt zwei Listen eine in der die Objekte verwaltet werden,
 die sogenannte Masterliste und eine zweite, "uber die lediglich
 Zugriff auf die Objekte genommen wird.

 Bei der Benutzung sollte darauf geachtet werden, da"s nur die Masterliste
 den entsprechenden Flag besitzt.
*/
class T_CContainerObjectMaster : public T_CContainerObject {
 public:
  //! Defaultkonstruktor
  T_CContainerObjectMaster(void);
  //! Defaultdekonstruktor
  virtual ~T_CContainerObjectMaster(void){};
  //! setzt die Identifikaion des Listenelementes
  inline void setIdMaster(T_Id_t Id);
  //! liefert die Identifikaion des Listenelementes
  inline T_Id_t getIdMaster(void);
  //! entfernt dieses Objekt aus der verketteten Liste
  inline void unlinkMaster();
  //! f"ugt dieses Objekt in eine verketteten Liste ein
  inline void linkWithMaster(T_CContainer *pCC);
 private:
  //! die Einbindung des Objektes in die Masterliste erfolgt "uber dieses Kettenglied
  T_CContainerItem m_CCI;
};

/*************************************************************************/


//! Iterator f"ur bidirektional verkettete Listen
class T_CContainerIterator {
 friend class T_CContainer;
 public:
  //! Default Konstruktor
  inline T_CContainerIterator(void);
  //! Konstruktor, setzt den Iterator an den Anfang der "ubergebenen Liste
  T_CContainerIterator(T_CContainer *pCC);
  //! Konstruktor, setzt den Iterator an den Anfang der "ubergebenen Liste
  T_CContainerIterator(T_CContainer &CC);
  //! setzt den Iterator an den Anfang der "ubergebenen verketteten Liste
  inline void begin(T_CContainer *pCC);
  //! setzt den Iterator an den Anfang der "ubergebenen verketteten Liste
  inline void begin(T_CContainer &CC);
  //! setzt den Iterator an das Ende der "ubergebenen verketteten Liste
  inline void end(T_CContainer *pCC);
  //! setzt den Iterator an das Ende der "ubergebenen verketteten Liste
  inline void end(T_CContainer &CC);
  //! "uberpr"uft ob sich der Iterator am Anfang oder Ende der verketten Liste befindet
  inline bool isValid(void) const;
  //! vergleicht ob zwei Iteratoren auf die selbe Stelle in der selben Liste zeigen
  inline bool operator==(const T_CContainerIterator &CCI) const;
  //! "Ubernahme der Liste und Position in der Liste
  inline T_CContainerIterator &operator=(const T_CContainerIterator &CCI);
  //! Springt in der Liste um ein Element vorw"arts
  inline void operator++(void);
  //! Springt in der Liste um ein Element r"uckw"arts
  inline void operator--(void);
  //! gibt einen Pointer des dem aktuellen Listenelement zugeordneten Objektes zur"uck
  inline T_CContainerObject *get(void);
 private:
  //! gibt das Listenelement zur"uck, auf dass der Iterator gerade zeigt.
  inline T_CContainerItem *getItem(void) const;
 private:
  //! Kettenglied auf dem der Iterator derzeit sitzt
  T_CContainerItem *m_pCCI;
};

/*************************************************************************/

//! bidirektional verkette Liste
class T_CContainer {
 friend class T_CContainerItem;
 friend class T_CContainerIterator;
 public:
  //! Konstruktor
  T_CContainer(void);
  //! Destruktor
  virtual ~T_CContainer(void);
  //! lehrt die verkettete Liste
  void clear(void);
  //! absorbiert die "ubergebene Liste
  void absorb(T_CContainer *pCC);
  //! Gibt das Listenelement mit Identifikation Id zur"uck.
  T_CContainerItem *getItem(T_Id_t Id);
  //! Gibt die Anzahl der Elemente in der Liste aus.
  T_N_t countItems(void);
  //! Gibt das letzte Objekt in der Liste aus 
  inline T_CContainerObject *getLast(void);
  //! markiert die Liste als Hauptverwalter der Listenelemente
  inline void setMaster(void);
  //! markiert die Liste als Hauptverwalter der Kettenelemente der Liste
  inline void setDeleteItems(void);
  //! dient dem debuggen, gibt Position im Source und die Zahl der enthaltenen Elemente aus
  inline void dg(void){DG;fprintf(stderr,"Container : %"T_N_PRINTF"\n",countItems());};
 private:
  //! Setzt den Beginn der verketteten Liste
  inline void setBeginAnchor(T_CContainerItem *pCCI);
  //! Setzt das Ende der verketteten Liste
  inline void setEndAnchor(T_CContainerItem *pCCI);
  //! Gibt das erste Listenelement zur"uck.
  inline T_CContainerItem *getBeginAnchor(void) const;
  //! Gibt das letzte Listenelement zur"uck.
  inline T_CContainerItem *getEndAnchor(void) const;
 private:
  //! Zeiger auf das erste Kettenglied der verwalteten verketteten Liste
  T_CContainerItem *m_pCCIBeginAnchor;
  //! Zeiger auf das letzte Kettenglied der verwalteten verketteten Liste
  T_CContainerItem *m_pCCIEndAnchor;
  //! ist dieses Flag gesetzt versteht sich die Liste als Masterliste
  /*!
   Als Masterliste wendet sie den delete Operator auf alle Kettenelemente
   und deren zugeordneten Objekte an, wenn der destruktor der Liste 
   aufgerufen wurde.
  */
  bool m_bMaster;
  //! ist dieses Flag gesetzt versteht sich die Liste als Verwalter der Kettenelemente
  /*!
   Normaler Weise sind die in der Liste verketteten Items Membervariablen.
   Es macht daher keinen Sinn einen Destruktor seitens der verketteten
   Liste anzuwenden. Ist dies jedoch nicht der Fall mu"ss der Destruktor
   der Elemente aufgerufen werden um s"amtlichen Speicher frei zu machen.
   Wird dieses Flag gesetzt, wird der Destruktor der Kettenelemente
   (T_CContainerItem) aufgerufen.
  */ 
  bool m_bDeleteItems;
};

#include "T_CContainer_inline.h"

//! 3D Array aus Pointern der Klasse T_CContainer
typedef T_CContainer* T_CContainerPtr3x3x3[3][3][3];

#endif
