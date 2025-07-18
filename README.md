# CMesure
Uncertainties calculation C++ class


## Variance propagation

| Function | Value | Uncertaincy |
| :------------: | :-------------: | :-------------: |
| $$F = K$$         | $$K$$ | $$0$$ |
| $$F = K.x$$         | $$K.\bar x$$ | $$\sigma_x$$ |
| $$F = x + y$$       | $$\bar x + \bar y$$ | $$\sigma_F^2 = \sigma_x^2 + \sigma_y^2$$ |
| $$F = x - y$$       | $$\bar x - \bar y$$ | $$\sigma_F^2 = \sigma_x^2 + \sigma_y^2$$ |
| $$F = x.y$$         | $$\bar x . \bar y$$ | $$\sigma_F^2 = \sigma_x^2 . y^2 + x^2 . \sigma_y^2$$ |
| $$F = \frac{x}{y}$$ |  $$\frac{\bar x}{\bar y}$$| $$\sigma_F^2 = (  (\sigma_x^2 . y^2) + (x^2 . \sigma_y^2)  ) . \frac{1}{y^4}$$ |

## Extended uncertaincy (a.k.a "Intervalle de Tolérence" in french)

$$IT = \bar \sigma_V.K$$
Wherre K is the number of standard deviation to take in account regarding the "Normal distribution"

## Comparisons and ordering

	-> Evaluate R = A - B
	-> Check where is the extented uncertaincy regarding the result


			    -IT(A-B)      0      +IT(A-B)
	-inf ------------------+----------+----------+-----------------> (A - B)

		(A!=B)			(A==B)			(A!=B) 
	-inf ------------------[----------+----------]------------------ +inf

		(A<=B)						(A>B)
	-inf ----------------------------------------]------------------ +inf

		(A<B)						(A>=B) 
	-inf ------------------[---------------------------------------- +inf


## Variance propagation for fucntions
| Function | Uncertaincy |
| :------------: | :-------------: |
| $$y = f(x)$$    | $$\sigma_y^2 = \left(\frac{\mathrm{d}f(x)}{\mathrm{d}x}\right)^2 \sigma_x^2$$ |
| $$z = f(x,y)$$  | $$\sigma_z^2 = \left(\frac{\mathrm{d}f(x)}{\mathrm{d}x}\right)^2 \sigma_x^2 + \left(\frac{\mathrm{d}f(y)}{\mathrm{d}y}\right)^2 \sigma_y^2$$ |
| $$y = \left\| x \right\|$$         | $$\sigma_y^2=\sigma_x^2$$ |
| $$y = \sin(x) $$    | $$\sigma_y^2 = \cos(x)^2  \sigma_x^2$$ |
| $$y = \cos(x) $$    | $$\sigma_y^2 = \left(-\sin(x) \right)^2  \sigma_x^2$$ |
| $$y = \tan(x) $$    | $$\sigma_y^2 = \left[1.0 + \tan^2(x) \right]^2  \sigma_x^2$$ |
| $$y = \arccos(x) $$ | $$\sigma_y^2 = \frac{1}{1-x^2} . \sigma_x^2$$ |
| $$y = \arcsin(x) $$ | $$\sigma_y^2 = \frac{1}{1-x^2} . \sigma_x^2$$ |
| $$y = \arctan(x) $$ | $$\sigma_y^2 = \left(\frac{1}{1-x^2} \right)^2 . \sigma_x^2$$ |
| $$y = cosh(x) $$    | $$\sigma_y^2 = sinh(x)^2  \sigma_x^2$$ |
| $$y = sinh(x) $$    | $$\sigma_y^2 = cosh(x)^2  \sigma_x^2$$ |
| $$y = tanh(x) $$    | $$\sigma_y^2 = \left(1.0 + tanh(x)^2 \right)^2  \sigma_x^2$$ |
| $$y = \log(x) $$ | $$\sigma_y^2 = \left(\frac{1}{x} \right)^2  \sigma_x^2$$ |
| $$y = \log_{10}(x) $$ | $$\sigma_y^2 = \left[\left(\frac{1}{log(10)} \right)\left(\frac{1}{x} \right)\right]^2  \sigma_x^2$$ |
| $$y = \log_{2}(x) $$ | $$\sigma_y^2 = \left[\left(\frac{1}{log(2)} \right)\left(\frac{1}{x} \right)\right]^2  \sigma_x^2$$ |
| $$y = e^x $$ | $$\sigma_y^2 = \left( e^x \right)^2  \sigma_x^2$$ |
| $$y = \sqrt x $$ | $$\sigma_y^2 = \left( \frac{1}{4\left\| x \right\|} \right)  \sigma_x^2$$ |
| $$y = x^y$$         | $$\sigma_z^2 = \left(y.x^{y-1}\right)^2 \sigma_x^2 + \left[(1+ln(x)).x^y \right]^2 \sigma_y^2$$ |

## Code example

```cpp - C++
#include <CMesure.h>

int main(void)
{

    CMesure fc(250000.0);	// 250 kHz
    CMesure Ct(0.000000001, 10.0,'P'); // 1nF +/- 10%

    CMesure Rt = 1.0 / (fc * Ct);
	
    cout << "fc =" << fc << endl;
	cout << "Ct =" << Ct << endl;
	cout << "Rt =" << Rt << endl;

	return 0;
}
```

```
fc =( 250000 +/- 2.98023e-08 | 95.45% )
Ct =( 1e-09 +/- 1e-10 | 95.45% )
Rt =( 4000 +/- 400 | 95.45% )
```

