#include "T_CContainer.h"

// Default Konstruktor
/*!
 Jedes Listenelement mu"s einem Objekt zugeordnet werden.
 Der Pointer auf dieses Objekt wird zur"uckgegeben, wenn
 die Liste durch einen Iterator abgeschritten wird.
*/
T_CContainerItem::T_CContainerItem(T_CContainerObject *pCCObject)
{
 m_pCCObject=pCCObject;
 m_pCCINext=NULL;
 m_pCCIPrev=NULL;
 m_pCC=NULL;
 m_Id=T_ID_DEFAULT;
}

// Default Destruktor
/*!
 Wird das Element zerst"ort, wird es automatisch aus der verketteten
 Liste ausgetragen in der es sich befindet.
*/
T_CContainerItem::~T_CContainerItem()
{
 m_pCCObject=NULL;
 unlink();
}

// entfernt dieses Element aus der verketteten Liste
/*
*/
void T_CContainerItem::unlink()
{
 if(!m_pCC)return;
 if(m_pCCIPrev)m_pCCIPrev->setNext(m_pCCINext);
 else m_pCC->setBeginAnchor(m_pCCINext);
 if(m_pCCINext)m_pCCINext->setPrev(m_pCCIPrev);
 else m_pCC->setEndAnchor(m_pCCIPrev);
 m_pCCINext=NULL;
 m_pCCIPrev=NULL;
 m_pCC=NULL;
}

// f"ugt dieses Element in eine verketteten Liste ein
/*!
 Bevor ein Element in eine verkettete Liste eingetragen wird, wird
 sie gegebenenfalls zuerst aus der alten Liste ausgetragen.
 Das Elment wird stets am Ende der verketteten Liste eingetragen.
*/
void T_CContainerItem::linkWith(T_CContainer *pCC)
{
 if(!pCC)return;
 if(m_pCC)unlink();
 m_pCCINext=NULL;
 m_pCCIPrev=pCC->getEndAnchor();
 if(m_pCCIPrev)m_pCCIPrev->setNext(this);
 else pCC->setBeginAnchor(this);
 pCC->setEndAnchor(this);
 m_pCC=pCC;
}

/*************************************************************************/

// Defaultkonstruktor
/*
*/
T_CContainerObject::T_CContainerObject(void)
:m_CCI(this)
{
}

/*************************************************************************/

// Defaultkonstruktor
/*
*/
T_CContainerObjectMaster::T_CContainerObjectMaster(void)
:m_CCI(this)
{
}


/*************************************************************************/

// Default Konstruktor
/*
*/
//T_CContainerIterator::T_CContainerIterator(void)
//{
 //m_pCC=NULL;
 //m_pCCI=NULL;
//}

//  Konstruktor, setzt den Iterator an den Anfang der "ubergebenen Liste
/*
*/
T_CContainerIterator::T_CContainerIterator(T_CContainer *pCC)
{
 this->begin(pCC);
}

//  Konstruktor, setzt den Iterator an den Anfang der "ubergebenen Liste
/*
*/
T_CContainerIterator::T_CContainerIterator(T_CContainer &CC)
{
 this->begin(CC);
}

/*************************************************************************/

// Konstruktor
/*
*/
T_CContainer::T_CContainer(void)
{
 m_pCCIBeginAnchor=NULL;
 m_pCCIEndAnchor=NULL;
 m_bMaster=FALSE;
 m_bDeleteItems=FALSE;
}

// Destruktor
/*!
 Es werden automatisch alle Listenelemente aus der verketteten Liste
 ausgetragen.
*/
T_CContainer::~T_CContainer(void)
{
 clear();
}

// lehrt die verkettete Liste
/*
*/
void T_CContainer::clear(void)
{
 T_CContainerObject *pCCO;
 T_CContainerItem *pCCI;

 if(m_bMaster){
  while(m_pCCIBeginAnchor){
   pCCI=m_pCCIBeginAnchor;
   pCCO=pCCI->get();
   if(m_bDeleteItems){
    T_DELETE(pCCI);
   } else {
    pCCI->unlink();
   }
   T_DELETE(pCCO)
  }
 } else {
  if(m_pCCIBeginAnchor){
   if(m_bDeleteItems){
    while(m_pCCIBeginAnchor){
     T_CContainerItem *pCCI=m_pCCIBeginAnchor;
     T_DELETE(pCCI);
    }
   } else {
    fprintf(stderr,"Error : T_CContainer contains %"T_N_PRINTF"\n",countItems());
    while(m_pCCIBeginAnchor)m_pCCIBeginAnchor->unlink();
   }
  }
 }
}

// absorbiert die "ubergebene Liste
/*!
 alle Elemente der "ubergebenen Liste werden "ubernommen, 
 und ans Ende der bestehenden Liste angeh"angt, die
 "ubergebene Liste ist danach leer.
*/
void T_CContainer::absorb(T_CContainer *pCC){
 T_CContainerItem *pCCI;
 pCCI=pCC->getBeginAnchor();
 if(!pCCI)return;

 do{
  pCCI->setContainer(this);
  pCCI=pCCI->getNext();
 }while(pCCI);

 pCCI=pCC->getBeginAnchor();
 if(m_pCCIEndAnchor){
  m_pCCIEndAnchor->setNext(pCCI);
 } else {
  m_pCCIBeginAnchor=pCCI;
 }
 pCCI->setPrev(m_pCCIEndAnchor);
 m_pCCIEndAnchor=pCC->getEndAnchor();

 pCC->setBeginAnchor(NULL);
 pCC->setEndAnchor(NULL);
}

// Gibt das Listenelement mit Identifikation Id zur"uck.
/*!
 Dabei wird die verketten Liste mit dem Anfang beginnend
 nach dem ersten Element mit der
 Identifikation Id durchsucht und der Pointer des diesem Element
 zugeordenten Elements zur"uckgegeben. Wird keine Element gefunden
 wird ein NULL-Pointer zur"uckgegeben
 \warning Die Routine ist laaangsam.
*/
T_CContainerItem *T_CContainer::getItem(T_Id_t Id)
{
 T_CContainerIterator I;

 for(I.begin(this);I.isValid();++I){
  if((I.getItem())->getId()==Id)return I.getItem();
 }
 return NULL;
}

// Gibt die Anzahl der Elemente in der Liste aus.
/*!
 \warning Es wird wirklich jedes mal gez"ahlt wieviele Elemente in der Liste sind
*/
T_N_t T_CContainer::countItems(void)
{
 T_CContainerIterator I;
 T_N_t counter=0;

 for(I.begin(this);I.isValid();++I){
  counter++;
 }
 return counter;
}


