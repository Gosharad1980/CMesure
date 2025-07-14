# CMesure
Uncertaicies calculation C++ class


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
	-> Check where is the extented uncertaicy regarding the result


						-IT(A-B)      0      +IT(A-B)
	-inf ------------------+----------+----------+-----------------> (A - B)

				(A!=B)			    (A==B)                (A!=B) 
	-inf ------------------[----------+----------]------------------ +inf

						(A<=B)                             (A>B)
	-inf ----------------------------------------]------------------ +inf

				(A<B)                         (A>=B) 
	-inf ------------------[---------------------------------------- +inf

---
