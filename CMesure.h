// test de sécurité pour éviter l'inclusion multiple des définitions
#ifndef _CMESURE_H_
#define _CMESURE_H_

//__amd64__ __amd64 __x86_64__ __x86_64
#if defined(__amd64__) || defined(__amd64) || defined(__x86_64__) || defined(__x86_64)
	#include <iostream>
#endif
using namespace std;

///////////// Classe de gestion simultanée d'une mesure et de son incertitude ///////////

// Dans la suite, tous les calculs sont efectués en considérant les données indépendantes
// Il n'y a donc pas de gestion des covariances ... et allez vous faire foutre.


class  CMesure
{
private:

	double valeur;      // Valeur de la mesure
    double variance;	// Valeur du carré l'incertitude type
    double alpha;       // taux de confiance pour le calcul du coeff d'élargissement
                        // => par défaut à 95.45% équivalent à un coeff d'élargissement
						// de 2.00 en loi normale

public:

//////////////////////////// Constructeurs et destructeur ////////////////////////////

	// constructeurs

	// valeur = 0.0; variance = CMESURE_EPS^2.0;	alpha = 95.45;
    CMesure();
	// valeur = valeur;   variance = CMESURE_EPS^2.0;	alpha = 95.45;						   
    CMesure(double _valeur);
	// valeur = v;   variance = epsilon * epsilon;		alpha = fabs(a);
    CMesure(double _valeur, double _epsilon, double _alpha);

#if defined(__amd64__) || defined(__amd64) || defined(__x86_64__) || defined(__x86_64)
	// pour passer en paramètre la chaine ( valeur +/- IT | alpha% )
	CMesure(char* m);
#endif				   

	CMesure(double _valeur, double _it, char _loi = 'N', double _alpha = 95.45);
	// Dans le cadre de mesures effectuées dans des conditions bien identifiées,
	// il est possible d'estimer l'incertitude type directement à partir de
	// l'intervalle de tolérance à l'aide des lois suivante
	//
	//		1) 'R' : Résolution d'un indicateur numérique       : epsilon = it / rac(12.0)
	//		2) 'H' : Hystérésis tel que it = MAXI - MINI        : epsilon = it / rac(12.0)
	//		3) 'S' : évolution Sinusoïdale sur it = MAXI - MINI : epsilon = it / 1.4
	//		4) 'N' : loi Normale par défaut, K = 2              : epsilon = it / 2.0
	//		5) 'C' : appareil de classe +/- it                  : epsilon = it / rac(3.0)
	//	 	6) 'P' : appareil de classe it = p%					: epsilon = (valeur * p / 100.0) / K_alpha(95.45)
	//
	// REM : la loi N est la loi par défaut dans tout bon certificat d'étalonnage qui se respecte

	// constructeur de recopie
	CMesure(const CMesure& X); 

	// destructeur
    ~CMesure();                             



////////////////// Accesseurs aux données privées de la classe /////////////////

	// Read only
	double Val(void);		// LA mesure en cours de traitement
    double Alpha(void);		// Taux de confiance
	double Eps(void);		// Incertitude type.
	double Variance(void);	// Variance = carré de l'incertitude type	
	double IT(void);		// Intervalle de tolérance = Eps x K

    double Fx(double _input, bool _inv);	// Fonction de répartition
	double K(void);	// Coeff d'élargissement


    
    // surdéfinition des opérateurs
    
    CMesure  operator+ (const CMesure& M) const;
    CMesure& operator+=(const CMesure& M);

	CMesure  operator+ (const double&  V) const;
	CMesure& operator+=(const double&  V);

    CMesure  operator- (const CMesure& M) const;
    CMesure& operator-=(const CMesure& M);

	CMesure  operator- (const double&  V) const;
	CMesure& operator-=(const double&  V);

    CMesure  operator* (const CMesure& M) const;
    CMesure& operator*=(const CMesure& M);

	CMesure  operator* (const double&  V) const;
	CMesure& operator*=(const double&  V);

    CMesure  operator/ (const CMesure& M) const;
    CMesure& operator/=(const CMesure& M);

	CMesure  operator/ (const double&  V) const;
	CMesure& operator/=(const double&  V);

    CMesure& operator= (const CMesure& M);
	CMesure& operator= (const double&  V);

    bool operator==(const CMesure& M) const;
    bool operator!=(const CMesure& M) const;
    bool operator<=(const CMesure& M) const;
    bool operator>=(const CMesure& M) const;
    bool operator< (const CMesure& M) const;
    bool operator> (const CMesure& M) const;


////////////////////////// Fonctions amies de la classe /////////////////////////

	friend CMesure operator-(CMesure& M);
	friend CMesure operator+(double V, const CMesure& M);
	friend CMesure operator-(double V, const CMesure& M);
	friend CMesure operator*(double V, const CMesure& M);
	friend CMesure operator/(double V, const CMesure& M);

    // affichage au format (valeur +/- IG)
#if defined(__amd64__) || defined(__amd64) || defined(__x86_64__) || defined(__x86_64)
    friend ostream& operator<<(ostream& os, CMesure const &M);
#endif
	
	// fonctions mathématique perso
	friend CMesure Conjug(CMesure& M);
};

#endif	// fin de _CMESURE_H_
