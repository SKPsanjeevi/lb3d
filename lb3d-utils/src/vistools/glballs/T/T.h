#include "config.h"

#if !defined(T_H)
#define T_H

//#define _ISOC99_SOURCE
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <stdlib.h> // free()
#include <stdio.h> // stderr ...

#include "T_Types.h" 

/*!
 Datentyp f"ur Ortskoordinaten oder Distanzen
*/
typedef T_Real64_t T_R3_t; 
/*!
 Definiert den Wert Null f"ur Ortskoordinaten
*/
#define T_R3_ZERO T_Real64_ZERO
/*!
 Definiert den kleinsten darstellbaren Wert
*/
#define T_R3_MIN T_Real64_MIN
/*!
 Definiert den gr"o"sten darstellbaren Wert
*/
#define T_R3_MAX T_Real64_MAX
/*!
 Definiert das Ausgabeformat
*/
#define T_R3_SCANF T_Real64_SCANF
#define T_R3_PRINTF T_Real64_PRINTF
/*!
 Definiert welches MPI Handle benutzt werden muss
*/
#define T_R3_MPI T_Real64_MPI

/*!
 Datentyp f"ur Zeitkoordinaten oder Zeitspannen
*/
typedef T_Real64_t T_Time_t;
/*!
 Definiert den Wert Null f"ur Zeitkoordinaten
*/
#define T_TIME_ZERO T_Real64_ZERO
/*!
 Definiert den kleinsten darstellbaren Wert
*/
#define T_TIME_MIN T_Real64_MIN
/*!
 Definiert den gr"o"sten darstellbaren Wert
*/
#define T_TIME_MAX T_Real64_MAX
/*!
 Definiert das Ausgabeformat
*/
#define T_TIME_SCANF T_Real64_SCANF
#define T_TIME_PRINTF T_Real64_PRINTF
/*!
 Definiert welches MPI Handle benutzt werden muss
*/
#define T_TIME_MPI T_Real64_MPI


/*!
 Datentyp f"ur Massen
*/
typedef T_Real64_t T_Mass_t;
/*!
 Definiert den Wert Null f"ur Massen
*/
#define T_MASS_ZERO T_Real64_ZERO
/*!
 Definiert das Ausgabeformat
*/
#define T_MASS_SCANF T_Real64_SCANF
#define T_MASS_PRINTF T_Real64_PRINTF
/*!
 Definiert welches MPI Handle benutzt werden muss
*/
#define T_MASS_MPI T_Real64_MPI

/*!
 Datentyp f"ur Dichten
*/
typedef T_Real64_t T_Density_t;
/*!
 Definiert den Wert Null f"ur Massen
*/
#define T_DENSITY_ZERO T_Real64_ZERO
/*!
 Definiert das Ausgabeformat
*/
#define T_DENSITY_SCANF T_Real64_SCANF
#define T_DENSITY_PRINTF T_Real64_PRINTF
/*!
 Definiert welches MPI Handle benutzt werden muss
*/
#define T_DENSITY_MPI T_Real64_MPI

/*!
 Datentyp f"ur Energien
*/
typedef T_Real64_t T_Energy_t;
/*!
 Definiert den Wert Null f"ur Energien
*/
#define T_ENERGY_ZERO T_Real64_ZERO
/*!
 Definiert das Ausgabeformat
*/
#define T_ENERGY_SCANF T_Real64_SCANF
#define T_ENERGY_PRINTF T_Real64_PRINTF
/*!
 Definiert welches MPI Handle benutzt werden muss
*/
#define T_ENERGY_MPI T_Real64_MPI

/*!
 Datentyp f"ur die Angabe der Klasse der ein Teilchen zugeordnet ist
 innerhalb des Ablaufes der Simulation. Innerhalb der Simulation sind
 die bekannten Klassen aufsteigend nummeriert. Der kleinste erlaubte
 Wert ist eins.
*/
typedef T_Int32_t T_ParticleClassIntern_t;
/*!
 Das Teilchen ist keiner g"ultigen Klasse zugeordnet worden.
*/
#define T_PARTICLECLASSINTERN_INVALID 0

#define T_PARTICLECLASSINTERN_ZERO 0
/*!
 Definiert das Ausgabeformat
*/
#define T_PARTICLECLASSINTERN_SCANF T_Int32_SCANF
#define T_PARTICLECLASSINTERN_PRINTF T_Int32_PRINTF
/*!
 Definiert welches MPI Handle benutzt werden muss
*/
#define T_PARTICLECLASSINTERN_MPI T_Int32_MPI

/*!
 Datentyp f"ur die Angabe der Klasse der ein Teilchen zugeordnet ist
 ausserhalb des Programmes (z.B. Konfigurationsdateien). 
 Es k"onnen beliebige Werte zugeordnet werden, erst w"ahrend der Simulation
 wird gekl"art ob die f"ur die Teilchen angegebenen Klassen existieren. 
 Existiert die Klasse nicht wird das Teilchen ignoriert.
 Der kleinste erlaubte Wert ist eins. 
*/
typedef T_Int32_t T_ParticleClassExtern_t;
/*!
 Das Teilchen ist keiner g"ultigen Klasse zugeordnet worden.
*/
#define T_PARTICLECLASSEXTERN_INVALID 0
/*!
 Definiert das Ausgabeformat
*/
#define T_PARTICLECLASSEXTERN_SCANF T_Int32_SCANF
#define T_PARTICLECLASSEXTERN_PRINTF T_Int32_PRINTF
/*!
 Definiert welches MPI Handle benutzt werden muss
*/
#define T_PARTICLECLASSEXTERN_MPI T_Int32_MPI

/*!
 Datentyp f"ur die Angabe der Klasse der eine Wand zugeordnet ist
 innerhalb des Ablaufes der Simulation. Innerhalb der Simulation sind
 die bekannten Klassen aufsteigend nummeriert. Der kleinste erlaubte
 Wert ist eins.
*/
typedef T_Int32_t T_WallClassIntern_t;
/*!
 Das Teilchen ist keiner g"ultigen Klasse zugeordnet worden.
*/
#define T_WALLCLASSINTERN_INVALID 0
/*!
 Definiert das Ausgabeformat
*/
#define T_WALLCLASSINTERN_SCANF T_Int32_SCANF
#define T_WALLCLASSINTERN_PRINTF T_Int32_PRINTF
/*!
 Definiert welches MPI Handle benutzt werden muss
*/
#define T_WALLCLASSINTERN_MPI T_Int32_MPI

/*!
 Datentyp f"ur die Angabe der Klasse der ein Teilchen zugeordnet ist
 ausserhalb des Programmes (z.B. Konfigurationsdateien). 
 Es k"onnen beliebige Werte zugeordnet werden, erst w"ahrend der Simulation
 wird gekl"art ob die f"ur die Teilchen angegebenen Klassen existieren. 
 Existiert die Klasse nicht wird das Teilchen ignoriert.
 Der kleinste erlaubte Wert ist eins. 
*/
typedef T_Int32_t T_WallClassExtern_t;
/*!
 Das Teilchen ist keiner g"ultigen Klasse zugeordnet worden.
*/
#define T_WALLCLASSEXTERN_INVALID 0
/*!
 Definiert das Ausgabeformat
*/
#define T_WALLCLASSEXTERN_SCANF T_Int32_SCANF
#define T_WALLCLASSEXTERN_PRINTF T_Int32_PRINTF
/*!
 Definiert welches MPI Handle benutzt werden muss
*/
#define T_WALLCLASSEXTERN_MPI T_Int32_MPI

/*!
 Datentyp f"ur die Identifizierung einzelner Elemente einer verketteten Liste
*/
typedef T_Int32_t T_Id_t;
/*!
 Standardwert f"ur den Datentyp T_Id_t
*/
#define T_ID_DEFAULT 0
/*!
 Definiert das Ausgabeformat
*/
#define T_ID_SCANF T_Int32_SCANF
#define T_ID_PRINTF T_Int32_PRINTF
/*!
 Definiert welches MPI Handle benutzt werden muss
*/
#define T_ID_MPI T_Int32_MPI

/*!
 Datentyp f"ur die Anzahl der Elemente einer verketteten Liste
*/
typedef T_Int32_t T_N_t;
/*!
 Definiert den Wert Null f"ur Anzahlen
*/
#define T_N_ZERO 0
/*!
 Definiert das Ausgabeformat
*/
#define T_N_SCANF T_Int32_SCANF
#define T_N_PRINTF T_Int32_PRINTF
/*!
 Definiert welches MPI Handle benutzt werden muss
*/
#define T_N_MPI T_Int32_MPI


#include "T_Macros.h" 

#include "T_Debug.h"


#endif
