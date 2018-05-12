#ex1
p = 1:20
p_inv = rev(p)
soma = p + p_inv
add3 = p + 3
vetor1p = 1/p

#ex2
gasolina = c(60, 55, 38, 87, 65, 63, 43, 44, 45, 50, 78, 67)
soma = sum(gasolina)
names(gasolina) = month.name

minValue = min(gasolina)
minMes = min(names(gasolina))
maxValue = max(gasolina)
maxMes = max(names(gasolina))
media = mean(gasolina)