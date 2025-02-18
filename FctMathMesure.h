// test de sécurité pour éviter l'inclusion multiple des définitions
#ifndef _FCTMATHMESURE_H_
#define _FCTMATHMESURE_H_


// Applicable uniquement en x86_64
// Suppression de cette partie qui bug en arduino :-(
#if defined(__amd64__) || defined(__amd64) || defined(__x86_64__) || defined(__x86_64)



#include "CMesure.h"

///////////////////////////// Fonctions mathématiques //////////////////////////

// Les fonctions suivantes nécessite l'utilisation
// de la loi de propagation des incertitudes.
// cf norme NF ENV 13005 (Guide d'Utilisation de la Métrologie : GUM)

// 1 - ce groupe de fonction est à une seule variable,
// ce qui donne, pour l'équation y = f(x), la formule suivante:
// U²(y) = [df(x)/dx]² * U²(x)

	CMesure fabs (CMesure M);
	CMesure sin  (CMesure M);
	CMesure cos  (CMesure M);
	CMesure tan  (CMesure M);
	CMesure acos (CMesure M);
	CMesure asin (CMesure M);
	CMesure atan (CMesure M);
	CMesure cosh (CMesure M);
	CMesure sinh (CMesure M);
	CMesure tanh (CMesure M);
	CMesure log  (CMesure M);
	CMesure log10(CMesure M);
	CMesure log2 (CMesure M);
	CMesure exp  (CMesure M);
	CMesure sqrt (CMesure M);
	CMesure cbrt (CMesure M);

	CMesure ceil (CMesure M);
	CMesure floor(CMesure M);


// 2 - ce groupe de fonction est à deux variables,
// ce qui donne, pour l'équation z = f(x,y), la formule suivante:
// U²(z) = [df(x)/dx]² * U²(x) + [df(y)/dy]² * U²(y)

	CMesure pow  (CMesure Base, CMesure Puiss);
	CMesure pow  (double  Base, CMesure Puiss);
	CMesure pow  (CMesure Base, double   Puiss);

	CMesure atan2(CMesure X, CMesure Y);
	CMesure modf (CMesure X, CMesure* Fpart);

// 3 - ce groupe de fonction me pose problème pour déterminer
// comment calculer l'incertitude type. Je vais certainement faire
// un choix personnel unilatérale parfaitement arbitraire.

	CMesure MIN  (CMesure A, CMesure B);
	CMesure MAX  (CMesure A, CMesure B);
	CMesure SIGN (CMesure A);



#endif	// fin de protection x86_64

#endif	// fin de _FCTMATHMESURE_H_