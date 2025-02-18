#define _USE_MATH_DEFINES // M_LN10 M_PI

// Applicable uniquement en x86_64
// Suppression de cette partie qui bug en arduino :-(
#if defined(__amd64__) || defined(__amd64) || defined(__x86_64__) || defined(__x86_64)


#include <cmath>
#include <cfloat>


#include "CMesure.h"

#include "FctMathMesure.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
///////////////////////////// Fonctions mathématiques //////////////////////////
////////////////////////////////////////////////////////////////////////////////

// Les fonctions suivantes nécessite l'utilisation
// de la loi de propagation des incertitudes.
// cf norme NF ENV 13005 (Guide d'Utilisation de la Métrologie : GUM)


// 1 - ce groupe de fonction est à une seule variable,
// ce qui donne, pour l'équation y = f(x), la formule suivante:
// U²(y) = [df(x)/dx]² * U²(x)

CMesure fabs (CMesure M) 
{
	// si M.Val() <  0 => fabs(M.Val()) = -1.0 * M.Val() et df = -1
	// si M.Val() >= 0 => fabs(M.Val()) =  1.0 * M.Val() et df =  1
	// dans tous les cas, lorsque df sera élevé au carré, df² = 1
	// => U²(x) = sqrt(M.Eps() * M.Eps()) = fabs(M. Eps())
	
	return CMesure(
		(double)fabs(M.Val()), 
		M.Variance(), 
		'D', 
		M.Alpha());
}

CMesure sin  (CMesure M)
{
	// d[sin(x)] = cos(x)
	
	double df;
	df = (double)cos(M.Val());
	
	return CMesure(
		(double)sin(M.Val()),
		df * df * M.Variance(),
		'D',
		M.Alpha());
}

CMesure cos  (CMesure M)
{
	// d[cos(x)] = -sin(x)
	
	double df;
	df = (-1.0) * (double)sin(M.Val());
	
	return CMesure(
		(double)cos(M.Val()),
		df * df * M.Variance(),
		'D',
		M.Alpha());
}

CMesure tan  (CMesure M)
{
	// d[tan(x)] = 1 + tan²(x)
	
	double df, val;
	val = (double)tan(M.Val());
	df = 1.0 + (val * val);

	return CMesure(
		val,
		df * df * M.Variance(),
		'D',
		M.Alpha());
}

CMesure acos (CMesure M)
{
	// d[acos(x)] = -1 / sqrt(1 - x²) 
	// df² = (-1)² / sqrt(1 - x²)²
	// df² = 1 / (1 - x²)
	
	double df2;
	df2 = (double) ( 1.0 / (1.0 - (M.Val() * M.Val())));

	return CMesure(
		(double)acos(M.Val()),
		df2 * M.Variance(),
		'D',
		M.Alpha());
}

CMesure asin (CMesure M)
{
	// d[asin(x)] = 1 / sqrt(1 - x²)
	// df² = 1² / sqrt(1 - x²)²
	// df² = 1 / (1 - x²)
	
	double df2;
	df2 = (double) ( 1.0 / ( 1.0 - (M.Val() * M.Val()) ));

	return CMesure(
		(double)asin(M.Val()),
		df2 * M.Variance(),
		'D',
		M.Alpha());
}

CMesure atan (CMesure M)
{
	// d[atan(x)] = 1 / (1 - x²)
	
	double df;
	df = (double) ( (1.0) / (1.0 - M.Val() * M.Val()) );

	return CMesure(
		(double)atan(M.Val()),
		df * df * M.Variance(),
		'D',
		M.Alpha());
}

CMesure cosh (CMesure M)
{
	// d[acos(x)] = sinh(x)
	
	double df;
	df = (double) sinh(M.Val());

	return CMesure(
		(double)cosh(M.Val()), 
		df * df * M.Variance(),
		'D',
		M.Alpha());
}

CMesure sinh (CMesure M)
{
	// d[sinh(x)] = cosh(x)
	
	double df;
	df = (double) cosh(M.Val());

	return CMesure(
		(double)sinh(M.Val()), 
		df * df * M.Variance(),
		'D',
		M.Alpha());
}

CMesure tanh (CMesure M)
{
	// d[tanh(x)] = 1 + tanh²(x)
	
	double df, val;
	val = (double)tanh(M.Val());
	df = 1.0 + (val * val);

	return CMesure(
		val,
		df * df * M.Variance(),
		'D',
		M.Alpha());
}

CMesure log  (CMesure M)
{
	// d[log(x)] = 1 / x
	
	double df;
	df = 1.0 / M.Val();

	return CMesure(
		(double)log(M.Val()), 
		df * df * M.Variance(),
		'D',
		M.Alpha());
}

CMesure log10(CMesure M)
{
	// d[log10(x)] = d[log(x)/log(10)]
	//             = 1/log(10) * d[log(x)]
	//             = 1/log(10) * 1/x

	double df;
	//df = (1.0 / log(10.0)) * (1.0 / M.Val());
	df = (1.0 / M_LN10) * (1.0 / M.Val());

	return CMesure(
		(double)log10(M.Val()), 
		df * df * M.Variance(),
		'D',
		M.Alpha());
}

CMesure log2(CMesure M)
{
	// d[log10(x)] = d[log(x)/log(2)]
	//             = 1/log(2) * d[log(x)]
	//             = 1/log(2) * 1/x

	double df;
	//df = (1.0 / log(10.0)) * (1.0 / M.Val());
	df = (1.0 / M_LN2) * (1.0 / M.Val());

	return CMesure(
		(double)log2(M.Val()), 
		df * df * M.Variance(),
		'D',
		M.Alpha());
}

CMesure exp  (CMesure M)
{
	// d[exp(x)] = exp(x)

	double df;
	df = exp(M.Val());

	return CMesure(
		(double)exp(M.Val()), 
		df * df * M.Variance(),
		'D',
		M.Alpha());
}

CMesure sqrt (CMesure M)
{
	// d[rac(x)] = 1 / (2 * rac(x))
	// df² = 1² / (2 * rac(x))²
	// df² = 1 / (4 * |x|)

	double df2;
	df2 = 1.0 / (4.0 * fabs(M.Val()));
	
	return CMesure(
		(double)sqrt(M.Val()), 
		df2 * M.Variance(),
		'D',
		M.Alpha());
}

CMesure cbrt (CMesure M) { return pow(M, 1.0/3.0); }



// 2 - ce groupe de fonction est à deux variables,
// ce qui donne, pour l'équation z = f(x,y), la formule suivante:
// U²(z) = [df(x)/dx]² * U²(x) + [df(y)/dy]² * U²(y)

CMesure pow(CMesure M, CMesure Puiss)
{
	// d[pow(m,p)] = d[m^p]/dm + d[m^p]/dp
	//
	// d[m^p]/dm = p * m^(p-1)
	//
	// d[m^p]/dp = d[e^(p*ln(m))]
	//           = m^(p-1) * (m*1*ln(m) + 1*m)
	//           = m*(1+ln(m))*m^(p-1)
	//           = (1+ln(m))*m^p
	
	double dfM, dfP, val;
	val = (double)pow(M.Val(), Puiss.Val());
	dfM = Puiss.Val() * pow(M.Val(), Puiss.Val() - 1.0);
	dfP = (1.0 + log(M.Val())) * val;

	return CMesure(
		val, 
		(double)((dfM * dfM * M.Variance()) + (dfP * dfP * Puiss.Variance())), 
		'D',
		(M.Alpha() > Puiss.Alpha() ? M.Alpha() : Puiss.Alpha()));
}

CMesure pow  (double Base, CMesure Puiss) { return pow((CMesure) CMesure(Base), (CMesure) Puiss);          }
CMesure pow  (CMesure Base, double Puiss) { return pow((CMesure) Base,          (CMesure) CMesure(Puiss)); }

CMesure atan2(CMesure X, CMesure Y)      
{
	//return atan((CMesure) (Y / X) );

	CMesure Mzero(0.0);

	if ((X == Mzero) && (Y == Mzero)) 	{ return Mzero; 			 }
	else if (X >  Mzero)				{ return atan(Y / X); 		 }
	else if (Y >= Mzero)				{ return atan(Y / X) + M_PI; }
	else								{ return atan(Y / X) - M_PI; }
}

CMesure modf (CMesure X, CMesure* Fpart)
{
	(*Fpart) = X - floor((CMesure) X);
	return floor((CMesure) X);
}


// 3 - ce groupe de fonction me pose problème pour déterminer
// comment calculer l'incertitude type. Je vais certainement faire
// un choix personnel unilatérale parfaitement arbitraire.


// Le calcul d'epsilon pour floor et ceil a posé plusieurs questions:
//		1) Application d'un coeff de proportionnalité newVal/oldVal ?
//		2) Considérer qu'une valeur seuillée possède un epsilon nul 
//		   car c'est une valeure certaine ?
//		3) Augmenter epsilon de suffisament pour conserver l'ancien IT
//		   dans un nouveau centré en newVal ?
//
// Résultat de mes réflexions : 
//		Solution 1 : problème identifé si floor ou ceil retourne 0.0 car 
//		cela provoque un epsilon infini => REFUSEE
//		Solution 2 : supprimer de façon artificielle une incertitude sur
//		une mesure est contraire à la philosophie de cette classe => REFUSEE
//		Solution 3 : c'est la moins pire, selon moi, elle offre un compromis 
//		acceptable entre la philosophie de cette classe et les loies mathématiques
//		sous jacente

CMesure ceil (CMesure M)
{ 
	double valeur, epsilon, alpha;
	valeur  = ceil(M.Val());
	epsilon = (M.IT() + fabs(valeur - M.Val())) / M.K();
	alpha   = M.Alpha();
	return CMesure(valeur, epsilon, alpha);
}

CMesure floor(CMesure M)
{
	double valeur, epsilon, alpha;
	valeur  = floor(M.Val());
	epsilon = (M.IT() + fabs(M.Val() - valeur)) / M.K();
	alpha   = M.Alpha();
	return CMesure(valeur, epsilon, alpha);
}

// Pour une incertitude de mesure les MIN et MAX ne sont pas symétriques
// Etant donné que les inégalités font intervenir les IT, il y a trois cas :
// -> A<B
// -> A>B
// -> A potentiellement égal à B : donc A ni inff ni supp à B mais A.val != B.val
// Dans ce troisième cas, faut-il
// 1) retourner A ? (puisque indécidable donc potentiellement équivalent)
// 2) retourner B ? (puisque indécidable donc potentiellement équivalent)
// 3) une autre valeur => choix unilatéral personnel parfaitement arbitraire de retourner la moyenne des deux ... et je vous emmerde

CMesure MIN (CMesure A, CMesure B)
{ 
	CMesure R;

	if (A == B)	    { R = (A + B) / 2.0; }	// cas indécidable
	else if (A < B) { R = A; }				// A < B
	else            { R = B; }				// A > B
	
	return R; 
}

CMesure MAX (CMesure A, CMesure B) 
{ 
	CMesure R;

	if (A == B)	    { R = (A + B) / 2.0; }	// cas indécidable
	else if (A < B) { R = B; }				// A < B
	else            { R = A; }				// A > B
	
	return R; 
}

CMesure SIGN(CMesure A)
{ 
	return ((A.Val() >= 0) ? CMesure(1.0) : CMesure(-1.0)); 
}


#endif // fin protestion x86_64
