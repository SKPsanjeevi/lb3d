// setzt die Identifikaion des Listenelementes
/*
*/
void T_CContainerItem::setId(T_Id_t Id)
{
 m_Id=Id;
}

// liefert die Identifikaion des Listenelementes
/*
*/
T_Id_t T_CContainerItem::getId(void)
{
 return m_Id;
}

// meldet, ob das Item einer verketteten Liste zugeordnet wurde
/*
*/
bool T_CContainerItem::hasContainer(void)
{
 return (m_pCC!=NULL);
}

// ordnet dem Item einen neuen Container zu
/*
*/
void T_CContainerItem::setContainer(T_CContainer *pCC)
{
 m_pCC=pCC;
}

// liefert das referenzierte Objekt des Listenelementes
/*
*/
T_CContainerObject *T_CContainerItem::get(void)
{
 return m_pCCObject;
}

// setzt den Vorg"anger dieses Elementes im Sinne einer verketten Liste
/*
*/
void T_CContainerItem::setPrev(T_CContainerItem *pCCI)
{
 m_pCCIPrev=pCCI;
}

// setzt den Nachfolger dieses Elementes im Sinne einer verketten Liste
/*
*/
void T_CContainerItem::setNext(T_CContainerItem *pCCI)
{
 m_pCCINext=pCCI;
}

// liefert den Vorg"anger in der verketteten Liste
/*
*/
T_CContainerItem *T_CContainerItem::getPrev(void) const
{
 return m_pCCIPrev;
}

// liefert den Nachfolger in der verketteten Liste
/*
*/
T_CContainerItem *T_CContainerItem::getNext(void) const
{
 return m_pCCINext;
}

/*************************************************************************/

// setzt die Identifikaion des Listenelementes
/*
*/
void T_CContainerObject::setId(T_Id_t Id)
{
 m_CCI.setId(Id);
}

// liefert die Identifikaion des Listenelementes
/*
*/
T_Id_t T_CContainerObject::getId(void)
{
 return m_CCI.getId();
}

// entfernt dieses Objekt aus der verketteten Liste
/*
*/
void T_CContainerObject::unlink()
{
 m_CCI.unlink();
}

// f"ugt dieses Objekt in eine verketteten Liste ein
/*
*/
void T_CContainerObject::linkWith(T_CContainer *pCC)
{
 m_CCI.linkWith(pCC);
}

// meldet, ob das Item einer verketteten Liste zugeordnet wurde
bool T_CContainerObject::hasContainer(void)
{
 return m_CCI.hasContainer();
}

/*************************************************************************/

// setzt die Identifikaion des Listenelementes
/*
*/
void T_CContainerObjectMaster::setIdMaster(T_Id_t Id)
{
 m_CCI.setId(Id);
}

// liefert die Identifikaion des Listenelementes
/*
*/
T_Id_t T_CContainerObjectMaster::getIdMaster(void)
{
 return m_CCI.getId();
}

// entfernt dieses Objekt aus der verketteten Liste
/*
*/
void T_CContainerObjectMaster::unlinkMaster()
{
 m_CCI.unlink();
}

// f"ugt dieses Objekt in eine verketteten Liste ein
/*
*/
void T_CContainerObjectMaster::linkWithMaster(T_CContainer *pCC)
{
 m_CCI.linkWith(pCC);
}


/*************************************************************************/

// Default Konstruktor
/*
*/
T_CContainerIterator::T_CContainerIterator(void)
{
 m_pCCI=NULL;
}

// setzt den Iterator an den Anfang der "ubergebenen verketteten Liste
/*
*/
void T_CContainerIterator::begin(T_CContainer *pCC)
{
 m_pCCI=pCC->getBeginAnchor();
}

// setzt den Iterator an den Anfang der "ubergebenen verketteten Liste
/*
*/
void T_CContainerIterator::begin(T_CContainer &CC)
{
 m_pCCI=CC.getBeginAnchor();
}

// setzt den Iterator an das Ende der "ubergebenen verketteten Liste
/*
*/
void T_CContainerIterator::end(T_CContainer *pCC)
{
 m_pCCI=pCC->getEndAnchor();
}

// setzt den Iterator an das Ende der "ubergebenen verketteten Liste
/*
*/
void T_CContainerIterator::end(T_CContainer &CC)
{
 m_pCCI=CC.getEndAnchor();
}


// "uberpr"uft ob sich der Iterator am Anfang oder Ende der verketten Liste befindet
/*
*/
bool T_CContainerIterator::isValid(void) const
{
 return m_pCCI;
}

// vergleicht ob zwei Iteratoren auf die selbe Stelle in der selben Liste zeigen
/*
*/
bool T_CContainerIterator::operator==(const T_CContainerIterator &CCI) const
{
 return m_pCCI==CCI.getItem();
}

// "Ubernahme der Liste und Position in der Liste
/*
*/
T_CContainerIterator &T_CContainerIterator::operator=(const T_CContainerIterator &CCI)
{
 m_pCCI=CCI.getItem();
 return *this;
}

// Springt in der Liste um ein Element vorw"arts
/*
*/
void T_CContainerIterator::operator++(void)
{
 m_pCCI=m_pCCI->getNext();
}

// Springt in der Liste um ein Element r"uckw"arts
/*
*/
void T_CContainerIterator::operator--(void)
{
 m_pCCI=m_pCCI->getPrev();
}

// gibt einen Pointer des dem aktuellen Listenelement zugeordneten Objektes zur"uck
/*
*/
T_CContainerObject *T_CContainerIterator::get(void)
{
 return m_pCCI->get();
}

// gibt das Listenelement zur"uck, auf dass der Iterator gerade zeigt.
/*
*/
T_CContainerItem *T_CContainerIterator::getItem(void) const
{
 return m_pCCI;
}

/*************************************************************************/

// Setzt den Beginn der verketteten Liste
/*
*/
void T_CContainer::setBeginAnchor(T_CContainerItem *pCCI)
{
 m_pCCIBeginAnchor=pCCI;
}

// Setzt das Ende der verketteten Liste
/*
*/
void T_CContainer::setEndAnchor(T_CContainerItem *pCCI)
{
 m_pCCIEndAnchor=pCCI;
}

// Gibt das erste Listenelement zur"uck.
/*
*/
T_CContainerItem *T_CContainer::getBeginAnchor(void) const
{
 return m_pCCIBeginAnchor;
}

// Gibt das letzte Listenelement zur"uck.
/*
*/
T_CContainerItem *T_CContainer::getEndAnchor(void) const
{
 return m_pCCIEndAnchor;
}

// Gibt das letzte Objekt in der Liste aus
/*
*/
T_CContainerObject *T_CContainer::getLast(void)
{
 if(!m_pCCIEndAnchor)return NULL;
 return m_pCCIEndAnchor->get();
}

// markiert die Liste als Hauptverwalter der Listenelemente
/*!
 Das setzen des Flags bewirkt, da"s beim Aufruf des Destruktors
 die Destruktoren der verwalteten T_CContainerObject Instancen
 aufgerufen werden. Damit r"aumt sich der Inhalt der Listen selbst auf.
*/
void T_CContainer::setMaster(void)
{
 m_bMaster=TRUE;
}

// markiert die Liste als Hauptverwalter der Kettenelemente der Liste
/*!
 Normaler Weise sind die in der Liste verketteten Items Membervariablen.
 Es macht daher keinen Sinn einen Destruktor seitens der verketteten
 Liste anzuwenden. Ist dies jedoch nicht der Fall mu"ss der Destruktor
 der Elemente aufgerufen werden um s"amtlichen Speicher frei zu machen.
 Wird dieses Flag gesetzt, wird der Destruktor der Kettenelemente
 (T_CContainerItem) aufgerufen.
*/
void T_CContainer::setDeleteItems(void)
{
 m_bDeleteItems=TRUE;
}

