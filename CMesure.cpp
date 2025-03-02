// Cette ligne est spécifique à Visual C++
// pour pouvoir utiliser le sscanf
//#define _CRT_SECURE_NO_DEPRECATE

#if defined(__amd64__) || defined(__amd64) || defined(__x86_64__) || defined(__x86_64)
	#include <cmath>
	#include <cfloat>
#else
	#include <math.h>
	#include <float.h>
#endif

#include <stdio.h>	// sscanf

/*
=> float.h est utilisée pour :
	#define DBL_EPSILON     2.2204460492503131e-016 // smallest such that 1.0+DBL_EPSILON != 1.0
	#define DBL_MAX         1.7976931348623158e+308 // MAXI value
*/

#include "CMesure.h"

using namespace std;

//#define CMESURE_MAX (sqrt(DBL_MAX)/2.0)
//#define CMESURE_EPS (sqrt(DBL_EPSILON))

/////////////////////////////////////////////////////////////////////////////////
/////////////////// Pour éviter les problémes de compilation ////////////////////
/////////////////////////////////////////////////////////////////////////////////

static inline double _CMESURE_MAX_(void) { return (sqrt(DBL_MAX)/2.0); }
#define CMESURE_MAX _CMESURE_MAX_()

static inline double _CMESURE_EPS_(void) { return (sqrt(DBL_EPSILON)); }
#define CMESURE_EPS _CMESURE_EPS_()

static inline double MINI(const double a, const double b) { return (((a)<(b)) ? (a) : (b)); }
static inline double MAXI(const double a, const double b) { return (((a)>(b)) ? (a) : (b)); }
static inline double SIGN(const double a)				  { return (((a)>=0) ? (1) : (-1)); }


////////////////////////////////////////////////////////////////////////////////
///////////////////////// Constructeurs et destructeur /////////////////////////
////////////////////////////////////////////////////////////////////////////////

// constructeur par défaut
CMesure::CMesure()
{
	this->valeur  = 0.0;
	this->variance = CMESURE_EPS*CMESURE_EPS;
	this->alpha   = 95.45;
}

// constructeur pour une constante : epsilon MINI
CMesure::CMesure(double _valeur)
{
	this->valeur  = _valeur;
	this->variance = CMESURE_EPS*CMESURE_EPS;
	this->alpha   = 95.45;
}

// construction d'une mesure en passant tous les paramètres
CMesure::CMesure(double _valeur, double _epsilon, double _alpha)
{
	this->valeur  = _valeur;
	this->variance = _epsilon*_epsilon;
	this->alpha   = MINI( MAXI(_alpha, 0.0), 100.0);
}

#if defined(__amd64__) || defined(__amd64) || defined(__x86_64__) || defined(__x86_64)
// Permet de creer une mesure à partir d'une sauvegarde texte
CMesure::CMesure(char* m)                     
{
	double v, it, a, e;
	sscanf(m, "( %lf +/- %lf | %lf%% )", &v, &it, &a);
	this->valeur  = v;
	this->alpha   = MINI( MAXI(a, 0.0), 100.0);
	e = fabs(it/this->K());
	this->variance = e*e;
}
#endif

// permet de créer une mesure à partir d'une loi de distribution connue
CMesure::CMesure(double _valeur, double _it, char _loi, double _alpha)
{
	// Dans le cadre de mesures effectuées dans des conditions bien identifiées,
	// il est possible d'estimer l'incertitude type directement à partir de
	// l'intervalle de tolérance à l'aide des lois suivante
	//
	//		0) 'D' : données directes sans traitement           : variance = it
	//		1) 'R' : Résolution d'un indicateur numérique       : epsilon = it / rac(12.0)
	//		2) 'H' : Hystérésis tel que it = MAXI - MINI        : epsilon = it / rac(12.0)
	//		3) 'S' : évolution Sinusoïdale sur it = MAXI - MINI : epsilon = it / 1.4
	//		4) 'N' : loi Normale par défaut, K = 2              : epsilon = it / 2.0
	//		5) 'C' : appareil de classe +/- it                  : epsilon = it / rac(3.0)
	//	 	6) 'P' : appareil de classe it = p%					: epsilon = (valeur * p / 100.0) / K_alpha(95.45)

	this->valeur = _valeur;
	this->alpha  = _alpha;

	switch(_loi)
	{
		case 'D' : this->variance = fabs(_it);			break;
		case 'R' : this->variance = (_it * _it) / 12.0;	break;
		case 'H' : this->variance = (_it * _it) / 12.0;	break;
		case 'S' : this->variance = (_it * _it) / 1.96;	break;
		case 'C' : this->variance = (_it * _it) / 3.0;	break;
		case 'P' : 
			this->variance = (_valeur * fabs(_it) / 100.0) / (this->K_alpha(_alpha)); 
			this->variance = this->variance * this->variance;
			break;
		case 'N' : // c'est la loi par défaut dans tout bon certificat d'étalonnage qui se respecte
		default  : this->variance = (_it * _it) / 4.0;	break;	
	}

}


// construceteur de recopie
CMesure::CMesure(const CMesure& X)
{
	this->valeur    = X.valeur;
	this->variance  = X.variance;
	this->alpha     = X.alpha;
}

// destructeur => vide car il n'y a aucune allocations dynamiques 
CMesure::~CMesure() 
{
	// Pas d'allocations dynamiques ... donc rien à détruire.
}


////////////////////////////////////////////////////////////////////////////////
////////////////// Accesseurs aux données privées de la classe /////////////////
////////////////////////////////////////////////////////////////////////////////

double CMesure::Val     (void) { return this->valeur;            }
double CMesure::Alpha   (void) { return this->alpha;             }
double CMesure::Variance(void) { return this->variance;          }
double CMesure::Eps     (void) { return sqrt(this->Variance());  }
double CMesure::IT      (void) { return this->Eps() * this->K(); }

// Coeff d'élargissement calculé en fct de alpha
double CMesure::K_alpha(double _alpha)
{
	// Calcul par interpolation du coeff d'élargissement à l'aide
	// des valeurs décrites dans la norme "NF ENV 13005"
    double p[13] = { 99.95 , 99.73 , 99.30, 99.00 , 98.76 , 95.45 , 95.00 , 90.00 , 86.64 , 68.27 , 50.000 , 38.29 , 0.000 };
    double k[13] = { 3.500 , 3.000 , 2.698, 2.576 , 2.500 , 2.000 , 1.960 , 1.645 , 1.500 , 1.000 , 0.6745 , 0.500 , 0.000 };

    double a, b;
	int i;
    
	if     (_alpha >= p[0]) return k[0];
	else if(_alpha <= p[12]) return k[12];
	else
	{
		// Recherche du cadran dans lequel on se situe
		for(i=1; i<13; i++) if(_alpha >= p[i]) break;

		// Interpolation de la valeur du coefficient d'élargissement
		a = (k[i] - k[i-1]) / (p[i] - p[i-1]);
		b = k[i-1] - (a * p[i-1]);

		return (a * _alpha + b);
	}
}

double CMesure::K(void) { return this->K_alpha(this->Alpha()); }
    
////////////////////////////////////////////////////////////////////////////////
//////////////////////// surdéfinition des opérateurs //////////////////////////
////////////////////////////////////////////////////////////////////////////////

CMesure  CMesure::operator+ (const CMesure& M) const
{
	// U²(this + M) = U²(this) + U²(M)
	
	CMesure R;

	R.valeur   = this->valeur + M.valeur;
	R.variance = this->variance + M.variance;
	R.alpha    = MAXI(this->alpha, M.alpha);

    return R;
}

CMesure& CMesure::operator+=(const CMesure& M)
{
	(*this) = (*this) + M;
    return (*this);
}

CMesure  CMesure::operator+ (const double&  V) const
{
	// une constante est une mesure dont l'incertitude est minimale
	return (*this) + CMesure(V);
}

CMesure& CMesure::operator+=(const double&  V)
{
	// une constante est une mesure dont l'incertitude est minimale
	(*this) = ((*this) + CMesure(V));
	return (*this);
}

CMesure  CMesure::operator- (const CMesure& M) const
{
    // U²(this - M) = U²(this) + U²(M)
	
	CMesure R;

	R.valeur   = this->valeur - M.valeur;
	R.variance = this->variance + M.variance;
	R.alpha    = MAXI(this->alpha, M.alpha);

    return R;
}
CMesure& CMesure::operator-=(const CMesure& M)
{
	(*this) = (*this) - M;
    return (*this);    
}

CMesure  CMesure::operator- (const double&  V) const
{
	// une constante est une mesure dont l'incertitude est minimale
	return (*this) - CMesure(V);
}

CMesure& CMesure::operator-=(const double&  V)
{
	// une constante est une mesure dont l'incertitude est minimale
	(*this) = ((*this) - CMesure(V));
	return (*this);
}

CMesure  CMesure::operator* (const CMesure& M) const
{
    // U²(R) = U²(this) * M² + this² * U²(M)

	CMesure R;
	
	R.valeur   = this->valeur * M.valeur;
	R.variance = (this->valeur * this->valeur * M.variance) + (this->variance * M.valeur * M.valeur);
	R.alpha    = MAXI(this->alpha, M.alpha);

    return R;
}

CMesure& CMesure::operator*=(const CMesure& M)
{
 	(*this) = (*this) * M;
    return (*this);
}

CMesure  CMesure::operator* (const double&  V) const
{
	// une constante est une mesure dont l'incertitude est minimale
	return (*this) * CMesure(V);
}

CMesure& CMesure::operator*=(const double&  V)
{
	// une constante est une mesure dont l'incertitude est minimale
	(*this) = (*this) * CMesure(V);
	return (*this);
}

CMesure  CMesure::operator/ (const CMesure& M) const
{
	// U²(R) = (  U(this)² * M²) + (this² * U(M)²)  ) * (1 / M²) 
	// CAS DE LA DIVISION DE/PAR ZERO !!! (traite l'infinie comme une valeur)
	//		R.valeur = +/-inf si dénominateur nul
	//		eps = +inf si dénom est nul

	CMesure R;
	
	if(M == CMesure(0.0))	
	{
		R.valeur   = (double)SIGN(this->valeur) * CMESURE_MAX;
		R.variance = CMESURE_MAX;
		R.alpha    = MAXI(this->alpha, M.alpha);
	}
	else
	{					
		R.valeur   = (this->valeur / M.valeur);
		//((self.Val().powf(2.0_f64) * RMesure_rhs.Variance()) + (RMesure_rhs.Val().powf(2.0_f64) * self.Variance())) / RMesure_rhs.Val().powf(4.0_f64)
		//R.variance = ((pow(this->valeur, 2.0) * M.variance) + (pow(M.valeur, 2.0) * this->variance)) / (pow(M.valeur, 4.0));
		R.variance = ((this->valeur * this->valeur * M.variance) + (M.valeur * M.valeur * this->variance)) / (M.valeur * M.valeur * M.valeur * M.valeur);
		R.alpha    = MAXI(this->alpha, M.alpha);
	}

    return R;
}

CMesure& CMesure::operator/=(const CMesure& M)
{
	(*this) = (*this) / M;
    return (*this);
}

CMesure  CMesure::operator/ (const double&  V) const
{
	// une constante est une mesure dont l'incertitude vaut z�ro
	return (*this) / CMesure(V);
}

CMesure& CMesure::operator/=(const double&  V)
{
	// une constante est une mesure dont l'incertitude vaut z�ro
	(*this) = (*this) / CMesure(V);
	return (*this);
}

CMesure& CMesure::operator= (const CMesure& M)
{
    this->valeur   = M.valeur;
    this->variance = M.variance;
    this->alpha    = M.alpha;
    return (*this);
}
CMesure& CMesure::operator= (const double&  V)
{
	// une constante est une mesure dont l'incertitude est minimale
	(*this) = CMesure(V);
	return (*this);
}

/////////////////////////////////////////////////////////////////////////////////
//////////////////////// Fonctions de tests pour les VA /////////////////////////
/////////////////////////////////////////////////////////////////////////////////

//
// Dans cette classe, les mesures sont considérées comme des
// variables aléatoires, effectuer un test d'ordonnancement
// entre deux mesures revient à effectuer un test statistique
// entre deux variables aléatoires.
//

/////////////////////////////////////////////////////////////////////////////////
//                                                                             //
//           PRINCIPE DE RESOLUTION DES TESTS D'ORDONNANCEMENT                 //
//           -------------------------------------------------                 //
//                                                                             //
//   -> Calcul de R = A - B                                                    //
//   -> contrôle de la position du résultat par rapport à son propre IT        //
//                                                                             //
//                                                                             //
//                     -IT(A-B)       0      +IT(A-B)                          //
//  -inf ------------------+----------+----------+-----------------> (A - B)   //
//                                                                             //
//              (A!=B)             (A==B)                (A!=B)                //
//  -inf ------------------[----------+----------]------------------ +inf      //
//                                                                             //
//                      (A<=B)                             (A>B)               //
//  -inf ----------------------------------------]------------------ +inf      //
//                                                                             //
//               (A<B)                         (A>=B)                          //
//  -inf ------------------[---------------------------------------- +inf      //
//                                                                             //
//                                                                             //
/////////////////////////////////////////////////////////////////////////////////


bool CMesure::operator==(const CMesure& M) const
{ 
	//
	// Dans cette classe, les mesures sont considérées comme des
	// variables aléatoires : tester si A == B revient à effectuer,
	// ce que l'on nomme en statistique, un test bilatéral.
	//
	// Cela consiste à calculer la VA équivalente à la différence des
	// deux VA testées et à vérifier que sa moyenne est comprise dans
	// son propre intervalle de tolérance centré en zéro.
	//
	CMesure D;
	D = (*this) - M;
	return ( fabs(D.valeur) <= D.IT() );
}

bool CMesure::operator!=(const CMesure& M) const { return !((*this) == M); }
bool CMesure::operator<=(const CMesure& M) const { return !((*this) >  M); }
bool CMesure::operator>=(const CMesure& M) const { return !((*this) <  M); }

bool CMesure::operator< (const CMesure& M) const 
{
	// tester est-ce que (A < B) consiste à vérifier si
	// 1) A et B sont significativement diffentes
	// 2) A.valeur < B.valeur

	CMesure D;
	D = (*this) - M;
	return ( D.valeur < (-1.0 * D.IT()) );
}

bool CMesure::operator>( const CMesure& M) const
{
	// tester est-ce que (A > B) consiste à vérifier si
	// 1) A et B sont significativement diffentes
	// 2) A.valeur > B.valeur

	CMesure D;
	D = (*this) - M;
	return ( D.valeur > D.IT() );
}

/////////////////////////////////////////////////////////////////////////////////
////////////////////////// Fonctions amies de la classe /////////////////////////
/////////////////////////////////////////////////////////////////////////////////

CMesure operator- (CMesure& M) { return (CMesure(0.0) - M); }

CMesure operator+ (double V, const CMesure& M) { return CMesure(V) + M; }
CMesure operator- (double V, const CMesure& M) { return CMesure(V) - M; }
CMesure operator* (double V, const CMesure& M) { return CMesure(V) * M; }
CMesure operator/ (double V, const CMesure& M) { return CMesure(V) / M; }

#if defined(__amd64__) || defined(__amd64) || defined(__x86_64__) || defined(__x86_64)
ostream& operator<<(ostream& os, CMesure const &M)
{
	CMesure tmp = M;

	os << "( " << tmp.Val() << " +/- " << tmp.IT() << " | " << tmp.Alpha() << "% )";
	//os << "( " << tmp.valeur << " +/- " << (tmp.valeur * tmp.K()) << " | " << tmp.alpha << "% )";
	return os;
}
#endif

CMesure Conjug(CMesure& M) 
{
	CMesure tmp = CMesure(M.Val(), M.Eps(), M.Alpha()); 
	return tmp;
}


