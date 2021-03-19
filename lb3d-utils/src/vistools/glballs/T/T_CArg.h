#include "config.h"

#if !defined(T_CARG_H)
#define T_CARG_H

#include <string.h>
#include "T.h"
#include "T_CContainer.h"

//! Kommandozeilenparameter
/*!
 Diese Klasse verwaltet die Informationen f"ur jeweils einen
 Kommandozeilenparameter.
 \internal
*/
class T_CArgItem : public T_CContainerObject {
 public :
  //! Defaultkonstruktor
  T_CArgItem(const char *Class,const char *Name,const char *ShortName,const char *Value,const char *Description);
  //! Defaultdestruktor
  ~T_CArgItem(void);
  //! setzt den Wert dieses Parameters
  bool set(const char *Value);
  //! gibt den Wert dieses Parameters aus
  char *get(void) const;
  //! schreibt den Wert dieses Parameters in einen Stream
  void dump(FILE *out);
  //! schreibt die Kurzbeschreibung des Parameters in einen Stream
  void describe(FILE *out);
  //! pr"uft ob der Parameter auf die angegebene Klasse und den angegebenen Namen pa"st
  bool isOption(const char *Class,const char *Name) const;
  //! pr"uft ob der Parameter auf die angegebene Klasse pa"st
  bool isClass(const char *Class) const;
  //! pr"uft ob der Parameter auf den angegebenen Namen pa"st
  bool hasName(const char *Name) const;
 private :
  //! Klasse des Parameters
  char *m_Class;
  //! Name des Parameters
  char *m_Name;
  //! Kurzform des Namens des Parameters
  char *m_ShortName;
  //! der Wert auf den der Parameter gesetzt ist
  char *m_Value;
  //! Beschreibung der Bedeutung des Paramters (f"ur -?)
  char *m_Description;
};

/*************************************************************************/

//! Commandline Parser
/*! 
 Die Klasse implementiert das Einlesen von Parametern wahlweise
 "uber die Kommandozeile oder aus einer Datei. Es werden nur
 Kommandozeilenparameter akzeptiert, die vor dem Parsen durch das
 Programm registriert wurden. Jeder Parameter ist einer Klasse
 zugeordnet, so da"s Doppelbenennungen bez"uglich verschiedener
 Programmteile aufgel"ost werden k"onnen. Die Angabe eines Parameters
 folgt dem Schema ''-Klasse::Parametername=Wert''.
 Wird keine Klasse angegeben, wird statt dessen ein Defaultwert genommen.
 Bei der Registrierung werden f"unf Angaben (Strings) erwartet. 
 	- Die Klasse der der Parameter zugeordnet werden soll,
 	- der eigentliche Name des Parameters, 
 	- optional ein Kurzname,
 	- der Defaultwert des Parameters und 
 	- eine Kurzbeschreibung der Bedeutung des Parameters.

 Die Klasse T_CArg besitzt einige Parameter, die der
 Klasse ''arg'' zugeordnet sind.
 	- --help : die allseits bekannte Hilfe
 	- -dump : gibt alle Parameter am Ende des Programmes in die hier angegebene Datei aus
 	- -file : lie"st zus"atzlich zur Kommandozeile die Parameter aus einer Datei

 Es ist dabei zu beachten, da"s Werte die in der eingelesenen Datei gesetzt
 werden stehts von der Kommandozeile "uberschreiben werden.
*/
class T_CArg {
 public: 
  //! Defaultkonstruktor
  T_CArg(const char *DefaultClass);
  //! Defaultdestruktur
  ~T_CArg(void);
  //! registriert einen Parameter
  void registerParameter(const char *Class,const char *Name,const char *ShortName,const char *Value,const char *Description);
  //! parst die von der main() Funktion stammenden Kommandozeilenparameter
  bool parseArgs(int argc,char *argv[]);
  //! liefert den Wert eines Parameters als bool
  bool isValidParameter(const char *Class,const char *Name);
  //! setzt einen Parameter auf einen Wert
  bool setParameter(const char *Class,const char *Name,const char *Value);
  //! liefert den Wert eines Parameters als String
  void getParameter(const char *Class,const char *Name,char **pValue);
  //! liefert den Wert eines Parameters als T_Real32_t
  void getParameter(const char *Class,const char *Name,T_Real32_t *pValue);
  //! liefert den Wert eines Parameters als T_Real64_t
  void getParameter(const char *Class,const char *Name,T_Real64_t *pValue);
  //! liefert den Wert eines Parameters als T_Int32_t
  void getParameter(const char *Class,const char *Name,T_Int32_t *pValue);
  //! liefert den Wert eines Parameters als T_Int64_t
  void getParameter(const char *Class,const char *Name,T_Int64_t *pValue);
  //! liefert den Wert eines Parameters als bool
  void getParameter(const char *Class,const char *Name,bool *pValue);
  //! liefert den Wert eines Parameters Minimum..Maximum als T_Real64_t
  void getParameter(const char *Class,const char *Name,T_Real64_t *pMin,T_Real64_t *pMax);
  //! liefert den Wert eines Parameters Minimum..Maximum als T_Int32_t
  void getParameter(const char *Class,const char *Name,T_Int32_t *pMin,T_Int32_t *pMax);
  //! liefert den Wert eines Parameters Minimum..Maximum als T_Int64_t
  void getParameter(const char *Class,const char *Name,T_Int64_t *pMin,T_Int64_t *pMax);
 private:
  //! parst die von der main() Funktion stammenden Kommandozeilenparameter
  void parseArgsInternal(int argc,char *argv[],bool Flag);
  //! parst eine Datei
  void parseFile(bool Flag);
  //! gibt die Beschreibung aller Parameter auf die Standartausgabe
  void describe(void);
  //! schreibt alle Parameter in die "uber Kommandozeile angegebene Datei
  bool dumpArgs(void);
  //! schreibt alle Parameter in einen Stream
  void dumpArgs(FILE *out);
 private:
  //! verkettete Liste aller registrierten Parameter
  T_CContainer m_CCAI;
  //! Name der Klasse der als Defaultwert angenommen wird
  char *m_DefaultClass;
};

#endif
