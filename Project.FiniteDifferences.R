	#Grupo 10
	#Integrantes:
		#Anderson Rodrigues 11096215
		#André Schittini 21016814
		#Daniel Badajoz 11056607
		#Gustavo S. de Paula 21085415
		#Icaro Olinda 11063713
		#Luane Souza 11051413
		#Marcos V. F. de Oliveira 11067212
		#Mariana Assumpção 11109910

		#teste

options(digits=22)
Main = function(){
	#Extremos
	a = 1
	b = 1.5
	#Número de partições
	n = 10
	#Solução no extremo inicial
	Alpha = sol(a)
	#Solução no extremo final
	Beta = sol(b)
	#Tamanho total do intervalo
	h = h_creator(n,a,b)
	#Tolerância ao erro para os métodos de resolução dos sistemas matriciais
	tol=10^(-6)
	#Numero máximo de interações que serão executadas pelos metodos de resolução de sitemas matriciais
	N_max=10^10
	#Estimando os valores da solução y(t) nos pontos igualmente espaçados
	x = x_creator(n,h,a)
	#Calcula a o valor da função na malha de pontos
	ysol=Solucao(n,x)
	#Cria e preenche a matriz com os termos dependentes
	Matriz = Gerar.Matriz(n,h,x)
	#Cria e preenche o vetor com os termos independentes
	Vetor=Gerar.Vetor.Term.Inde(n,h,x,Alpha,Beta)
	#Cria uma matriz com os termos dependentes e independentes
	MatrizMerged = Matriz
	MatrizMerged<-cbind(MatrizMerged,Vetor)[,c(1:(length(MatrizMerged[1,])+1))]

	# Merged Matrix debug:
	 print("Exibindo a Matriz")
	 print(Matriz)
	 print("Exibindo o Vetor")
	 print(Vetor)
	# print("MatrizMerged")
	# print(MatrizMerged)

	#Solução utilizando metodos internos do R
	#yR = solve(Matriz,Vetor)

	#Solução utilizando metodo Jacobi
	yJacobi= Jacobi(n, h, x, Matriz, Vetor, tol, N_max)
	#Solução utilizando metodo de Eliminação de Gauss
	yEliminGauss = EliminGauss(MatrizMerged)
	#Solução utilizando metodo de Gauss Seidel
	yGaussSeidel = Gauss.Seidel(n, h, x, Matriz, Vetor, tol, N_max)

	#LAGRANGE
	#N_Lagrange = n-2
	#p_Lagrange = array(0,c(N_Lagrange))
	#for(i in 1:(n-1)){
	#		p_Lagrange[i] = Interpolador.Lagrange(N_Lagrange,n,x[i],x,ysol)
	#	}

	#Transposição do vetor coluna
	#dim(yR) = NULL

	#Soluções encontradas
	print("A solução da equação nos pontos x ysol:")
	print(ysol)

	#print("Solução R: ")
	#print(yR)
	print("Solução aproximada porJacobi: ")
	print(yJacobi)

	print("Solução aproximada por Eliminação de Gauss: ")
	print(yEliminGauss)

	print("Solução aproximada por Gauss Seidel: ")
	print(yGaussSeidel)

	#print("Solução Lagrange: ")
	#print(p_Lagrange)



	#Exibição em gráfico cartesiano sem os extremos 'a' e 'b'
	#plot(x,ysol, type="l", pch=1, lty=1, col="black",ann = FALSE)
	#lines(x,yR, type="o", lty=1,col="black")
	#lines(x,yJacobi, type="p", pch=3, col="blue")
	#lines(x,yGaussSeidel, type="p", pch=1, col="green")
	#lines(x,yEliminGauss, type="p", pch=5, col="black")
	#lines(x,p_Lagrange, type="o", lty=1, col="black")
	#Nomeia os eixos x e y na cor preta
	title(xlab="X", col.lab=rgb(0,0,0))
	title(ylab="Y", col.lab=rgb(0,0,0))


}

Qx = function(x){
	return(16*x^2+28+((6)/x^2))
	}

Rx = function(x){
	return(0)
	}
#Solução
	sol = function(x){return((x^3)*exp(2*x^2))}

h_creator=function(n,a,b){
	h = (b-a)/n
	return (h)
}
x_creator=function(n,h,a){
	x = array(0,c(n - 1))
	for(i in 1:(n - 1)){
		x[i] = a + (i*h)
	}
	return (x)
}
Solucao=function(n,x){
	ysol = array(0,c(n - 1))
	for(i in 1:(n - 1)){
		ysol[i] = sol(x[i])		#Qx(x[i])*sol(x[i]) + Rx(x[i])
	}
	return (ysol)
}
#Método para gerar a Matriz
Gerar.Matriz = function(n,h,x){
		#Definindo a matriz do sistema de equações
		M = array(0,c(n-1,n-1))

		#Preenchendo a diagonal principal da matriz com:
		for (i in 1 : (n-1)){
			M[i,i] = 2 + (h^2 * Qx(x[i]))
		}

		#Preenchendo as diagonais adjacentes com -1
		for(i in 1 : (n-2)){
			M[i+1,i] = -1
			M[i,i+1] = -1
		}
		return (M)
	}

#Método gerar vetor INDEPENDENTES
Gerar.Vetor.Term.Inde=function(n,h,x,Alpha,Beta){
		#CONSTRUINDO O VETOR DE TERMOS INDEPENDENTES
		#Definindo o vetor de soluções do sistema
		v = array(0,c(n-1,1))

		#Preenchendo o vetor solução
		for(i in 1:(n-1)){
			v[i,1] = -h^2*Rx(x[i])
		}
		v[1,1] = v[1,1] + Alpha
		v[n-1,1] = v[n-1,1] + Beta
		return (v)
	}

#Metodo de resolução de sistemas matriciais por Jacobi
Jacobi=function(n, h, x, Matriz, Vetor, tol, N_max){
		w_0 = array(0,c(n-1))	#Criação do vetor w_0 que é a solução de "chute"
		for(i in 1:(n-1)){
			w_0[i]= 0	#"Chute" inicial
		}
		#Criação dos vetores w
		w_old = array(0,c(n-1))
		w_new = array(0,c(n-1))

		for(i in 1:(n-1)){
			w_old[i] = w_0[i] + (2*tol)
			w_new[i] = w_0[i]
		}
		norma= abs(max(w_new - w_old))		#modulo da maior diferença entre os valores dos vetores w
		iter=0

		#Laço para cálculos do método
		while(norma>tol && iter <N_max){
			iter=iter+1
			w_old=w_new
			w_new[1]=(1/(2+h^2*Qx(x[1])))*(Vetor[1]+w_old[2])
			w_new[n-1]=(1/(2+h^2*Qx(x[n-1])))*(Vetor[n-1]+w_old[n-2])
			for(i in 2:(n-2)){
				w_new[i]=(1/(2+h^2*Qx(x[i])))*(Vetor[i]+w_old[i-1]+w_old[i+1])
			}
			norma= abs(max(w_new-w_old))
		}
		return(w_new)
	}

#Metodo de resolução de sistemas matriciais por Eliminação de Gauss
EliminGauss=function(MatrixtoSolve){
	#turn the inputed matrix in a triangular matrix
	triangularMatrix = EliminGauss.toTriangular(MatrixtoSolve)
	# # -- debug
	# print("triangularMatrix: ")
	# print(triangularMatrix)

	# solves the linear system
	solvedSystem = EliminGauss.solveSystem(triangularMatrix)
	# # -- debug
	# print("solvedSystem: ")
	# print(solvedSystem)
	dim(solvedSystem) = NULL

	return(solvedSystem)
}

# Transforma a matriz em uma matriz triangular
EliminGauss.toTriangular=function(M){

	matrixLines= length(M[,1])    # length of first column in M
	matrixColumns = length(M[1,]) #length of first line in M

	# -- lines iteration
	for (i in 1:matrixLines){

		if(i==1){   # do nothing in first line
		next
		}

		# -- columns iteration
		for(j in 1:matrixColumns){

			# It's not a pivot
			if(i<=j){
				next
			}

			# It's a pivot:
			pivot = M[i,j]/M[j,j]

			# Calculate each column for actual pivot
			for(k in j:length(M[1,])){
				M[i,k]= M[i,k]-pivot*M[j,k]
			}

			# iteration debug:
			# strLine = "line"
			# strCol = "column"
			# strAux = paste(strLine,i,",",strCol,j)
			# print(strAux)
			# print(paste("pivot",pivot))
			# print(M);

		}
	}

	return(M)
}


#Encontra os valores de cada variável
EliminGauss.solveSystem=function(M){
 	# length of first column in M
	matrixLines= length(M[,1])
	# length of first line in M
	matrixColumns = length(M[1,])
	# number of variables on system
	numVar = matrixLines
	# column of independent terms
	indepTerm = matrixColumns

	# Create and start fill each element in output vector with zero
	outVector= c(1: numVar)
	for(el in outVector){
		outVector[el]=0
	}

	# Iterate lines in matrix from down to up
	for(var in numVar:1){

		# store the independent term
		outVector[var] = M[var,indepTerm]

		# solve the last line
		if(var>=numVar){
			outVector[var] = outVector[var]/M[var,var]
		next
		}

		# Iterate column in matrix from right to left
		for(i in numVar:(var+1)){

			outVector[var] = outVector[var] - M[var,i]*outVector[i]

			# iteration debug:
			# strLine = "line"
			# strCol = "column"
			# strAux = paste(strLine,var,",",strCol,i)
			# print(strAux)
		}

		outVector[var] = outVector[var]/M[var,var]

		# iteration output debug:
		# print(outVector)
	}

	return(matrix(outVector,numVar,1,FALSE))

}

#Metodo de resolução de sistemas matriciais por Gauss Seidel
Gauss.Seidel=function(n, h, x, Matriz, Vetor, tol, N_max){
		w_0 = array(0,c(n-1))
		for(i in 1:(n-1)){
			w_0[i]= 0
		}
		w_old = array(0,c(n-1))
		w_new = array(0,c(n-1))
		for(i in 1:(n-1)){
			w_old[i] = w_0[i] + (2*tol)
			w_new[i] = w_0[i]
		}
		norma= abs(max(w_new - w_old))
		iter=0
		while(norma>tol && iter <N_max){
			iter=iter+1
			w_old=w_new
			w_new[1]=(1/(2+h^2*Qx(x[1])))*(Vetor[1]+w_old[2])
			for(i in 2:(n-2)){
				w_new[i]=(1/(2+h^2*Qx(x[i])))*(Vetor[i]+w_old[i-1]+w_old[i+1])
			}
			w_new[n-1]=(1/(2+h^2*Qx(x[n-1])))*(Vetor[n-1]+w_old[n-2])
			norma= abs(max(w_new-w_old))
		}
		return(w_new)
}

#Método para criação do polinomio de Lagrange
#sendo N_Lagrange o grau do polinomio, x o ponto, X vetor de pontos, y solução aproximada
Interpolador.Lagrange = function(N_Lagrange,n, x, X, y){
    N_Lagrange=N_Lagrange+1
	P_L_final=0
	grau=2
	#Laço para realizar a somatória dos L(x) calculados em P_L_parcial
    for(i in 1:(N_Lagrange)){
		P_L_parcial=1.0

		#Laço para realizar a produtoria e gerar os L(x)
	   for(j in 1:(N_Lagrange)){
          if(j!=i){
				P_L_parcial= P_L_parcial * (x-X[j])/(X[i]-X[j])
			}
       }

	   P_L_final=P_L_final + P_L_parcial*y[i]
	   grau=grau*2
	}
	#print("grau")
	#print(grau)
	return (P_L_final)
}


#ERROS | Análise em função de X
erros_x=function(a,b){
	n=array(0,c(8))
	for(i in 1:(8)){
		n[i]= 40 + 5*i
		}

	z=1
	#Tolerância ao erro para os métodos de resolução dos sistemas matriciais
	tol=10^(-9)
	#Numero máximo de interações que serão executadas pelos metodos de resolução de sitemas matriciais
	N_max=10^5
	while(z <= 8){
		En = array(0,c(n[z]-1))
		#Solução no extremo inicial
		Alpha = sol(a)
		#Solução no extremo final
		Beta = sol(b)
		#Tamanho total do intervalo
		h = h_creator(n[z],a,b)
		#Estimando os valores da solução y(t) nos pontos igualmente espaçados
		x = x_creator(n[z],h,a)
		ysol=Solucao(n[z],x)
		#Cria e preenche a matriz com os termos dependentes
		Matriz = Gerar.Matriz(n[z],h,x)
		#Cria e preenche o vetor com os termos independentes
		Vetor=Gerar.Vetor.Term.Inde(n[z],h,x,Alpha,Beta)
		#Cria uma matriz com os termos dependentes e independentes
		MatrizMerged = Matriz
		MatrizMerged<-cbind(MatrizMerged,Vetor)[,c(1:(length(MatrizMerged[1,])+1))]
		#Solução utilizando metodo Jacobi
		#yJacobi= Jacobi(n[z], h, x, Matriz, Vetor, tol, N_max)
		#Solução utilizando metodo de Eliminação de Gauss
		#yEliminGauss = EliminGauss(MatrizMerged)
		#Solução utilizando metodo de Gauss Seidel
		yGaussSeidel = Gauss.Seidel(n[z], h, x, Matriz, Vetor, tol, N_max)
		for(i in 1:n[z]-1){
			En[i] = abs(ysol[i] - yEliminGauss[i])
		}
		#En = log10(En)
		if(z==1){plot((x),((En)), type="p", pch=1, lty=2, col="black",ylim=c(0,4*10^-2),xlim=c(1,1.6),ann = FALSE)}
		if(z>1){lines(x,En, type="p", lty=1, col="black")}
		#En[z] = abs(ysol - yGaussSeidel)
		#En[z] = abs(ysol - yJacobi)
		z=z+1
	}

	#n = log10(n)
	#Exibição em gráfico cartesiano sem os extremos 'a' e 'b'
	#plot((n),((En)), type="p", pch=1, lty=2, col="black",ann = FALSE)
	#lines(x,p_Lagrange, type="o", lty=1, col="black")
	#Nomeia os eixos x e y na cor preta
	#title(xlab="X", col.lab=rgb(0,0,0))
	#title(ylab="ERRO", col.lab=rgb(0,0,0))
	#Minimos_Quadrados(En,n)
	#Minimos_Quadrados_log(En,n)
}
#ERROS | Análise EM FUNÇÃO DE N
erros_N=function(a,b){
	n=array(0,c(7))
	for(i in 1:(7)){
		n[i]=5 + 5*i
		}
	z=1
	En = array(0,c(7))
	#Tolerância ao erro para os métodos de resolução dos sistemas matriciais
	tol=10^(-9)
	#Numero máximo de interações que serão executadas pelos metodos de resolução de sitemas matriciais
	N_max=10^5
	while(z <= 7){
		#Solução no extremo inicial
		Alpha = sol(a)
		#Solução no extremo final
		Beta = sol(b)
		#Tamanho total do intervalo
		h = h_creator(n[z],a,b)
		#Estimando os valores da solução y(t) nos pontos igualmente espaçados
		x = x_creator(n[z],h,a)
		ysol=Solucao(n[z],x)
		#Cria e preenche a matriz com os termos dependentes
		Matriz = Gerar.Matriz(n[z],h,x)
		#Cria e preenche o vetor com os termos independentes
		Vetor=Gerar.Vetor.Term.Inde(n[z],h,x,Alpha,Beta)
		#Cria uma matriz com os termos dependentes e independentes
		MatrizMerged = Matriz
		MatrizMerged<-cbind(MatrizMerged,Vetor)[,c(1:(length(MatrizMerged[1,])+1))]
		#Solução utilizando metodo Jacobi
		#yJacobi= Jacobi(n[z], h, x, Matriz, Vetor, tol, N_max)
		#Solução utilizando metodo de Eliminação de Gauss
		yEliminGauss = EliminGauss(MatrizMerged)
		#Solução utilizando metodo de Gauss Seidel
		#yGaussSeidel = Gauss.Seidel(n[z], h, x, Matriz, Vetor, tol, N_max)

		#ERRO ESTAVA AQUI
		En[z] = max(abs(yEliminGauss - ysol)) #OBS: O ERRO ESTAVA NESSA SUBTRAÇÃO, SENDO QUE A MESMA ESTAVA INVERTIDA
		#print("yEliminGauss - ysol")
		#print(yEliminGauss - ysol)
		#print(abs(yEliminGauss - ysol))
		#print(max(abs(yEliminGauss - ysol)))
		z=z+1
	}
	print("Valores do Erro:")
	for(d in 1:(7)){
		print(En[d])
		}

	#En = log10(En)
	#n = log10(n)

	#Exibição em gráfico cartesiano sem os extremos 'a' e 'b'
	#plot((n),((En)), type="p", pch=1, lty=2, col="black",ann = FALSE)
	#lines(x,p_Lagrange, type="o", lty=1, col="black")
	#Nomeia os eixos x e y na cor preta
	#title(xlab="N", col.lab=rgb(0,0,0))
	#title(ylab="ERRO", col.lab=rgb(0,0,0))

	#title(xlab="Log10 N", col.lab=rgb(0,0,0))
	#title(ylab="Log10 ERRO", col.lab=rgb(0,0,0))

	Minimos_Quadrados(En,n)
	#Minimos_Quadrados_log(En,n)
}


Minimos_Quadrados=function(En,Nerro){
	f_1=function(x){return (1/x)}
	f_2=function(x){return (1/x^2)}
	f_3=function(x){return (1/x^3)}
	tab = array(0,c(length(En),3))
	for(i in 1: length(En)){
		tab[i,1]=f_1(Nerro[i])
		tab[i,2]=f_2(Nerro[i])
		tab[i,3]=f_3(Nerro[i])
	}
	M=array(0,c(3,3))
	w=array(0,c(3,1))
	for(i in 1:3){
		for(j in 1:3){
			for(k in 1:length(En)){
				M[i,j] = M[i,j] + tab[k,i]*tab[k,j]
			}
		}
	}
	for(i in 1:3){
		for(k in 1:length(En)){
			w[i] = w[i] + tab[k,i]*En[k]
		}
	}
	coef = solve(M,w)

	p=function(x){return(coef[1,1]/x + coef[2,1]/x^2 + coef[3,1]/x^3)}
	#plot((Nerro),(En),type="p", pch=1, lty=2, col="red",ann = FALSE)
	par(new=TRUE)
	curve(p, from=10, to=40, col="blue")
	#title(xlab="N", col.lab=rgb(0,0,0))
	#title(ylab="ERRO", col.lab=rgb(0,0,0))

	#Valores dos Coeficientes da equação fi
	print("coef[1,1]")
	print(coef[1,1])
	print("coef[2,1]")
	print(coef[2,1])
	print("coef[3,1]")
	print(coef[3,1])

	#Calculando o valor de N para certos valores de tolerância do Erro
	TOL = 10^-5
	print("Raiz de n para TOL = 10^-5. Usando Metodo de Newton")
	print(Newton(10^4,10^-4,1,TOL,coef))
	print("Raiz para TOL = 10^-5. Usando Metodo da Bisseccao")
	print(Bisseccao(1,10^3,10^-6,TOL,coef))
}

Minimos_Quadrados_log=function(En,Nerro){
	f_1=function(x){return (1)}
	f_2=function(x){return (1/x^2)}
	tab = array(0,c(length(En),2))
	for(i in 1: length(En)){
		tab[i,1]=f_1(Nerro[i])
		tab[i,2]=f_2(Nerro[i])
	}
	M=array(0,c(2,2))
	w=array(0,c(2,1))
	for(i in 1:2){
		for(j in 1:2){
			for(k in 1:length(En)){
				M[i,j] = M[i,j] + tab[k,i]*tab[k,j]
			}
		}
	}
	for(i in 1:2){
		for(k in 1:length(En)){
			w[i] = w[i] + tab[k,i]*En[k]
		}
	}
	coef = solve(M,w)

	p=function(x){return(coef[1,1] - coef[2,1]*x)}
	par(new=TRUE)
	curve(p, from=10, to=40, col="blue")
	title(xlab="N", col.lab=rgb(0,0,0))
	title(ylab="ERRO", col.lab=rgb(0,0,0))


	#Valores do coeficiente da reta
	print("coef[1,1]")
	print(coef[1,1])
	print("coef[2,1]")
	print(coef[2,1])

	#Calculando o valor de N para certos valores de tolerância do Erro
	TOL = 10^-4
	print("Raiz para TOL = 10^-4. Usando Metodo de Newton")
	print(Newton(10^4,10^-6,1,TOL,coef))
	print("Raiz para TOL = 10^-4. Usando Metodo da Bisseccao")
	print(Bisseccao(1,10^4,10^-6,TOL,coef))
}


Newton=function(N_max,prec,X_0,TOL,coef){
	f=function(x){return(coef[1,1]/x + coef[2,1]/x^2 + coef[3,1]/x^3 - TOL)}
	derivF = function(x){return((-1)*coef[1,1]/x^2 + (-2)*coef[2,1]/x^3 + (-3)*coef[3,1]/x^4)}
	i=1
	x=X_0
	while(f(x-prec)*f(x+prec)> 0 && i < N_max){
		x = x-f(x)/derivF(x)
		i=i+1
	}
	return (x)
}

Bisseccao=function(m,M,tol,TOL,coef){
	f=function(x){return(coef[1,1]/x + coef[2,1]/x^2 + coef[3,1]/x^3 - TOL)}
	i = 0      # variável que conta o número de passos do algoritmo

	raiz = 0.5*(m+M)
	while(f(raiz-tol)*f(raiz+tol) > 0) # laço de aplicação do método. O método para quando o tamanho do intervalo de busca for
	{                           # inferior à tolerância tol desejada. Isso garante que a distância entre a raíz da equação
                             # e a aproximação raiz seja no máximo tol.
	raiz = 0.5*(m+M)
	if(f(raiz)*f(m) < 0)  # verifica se a raíz está no intervalo [m, raiz]
	{
		M = raiz
	}
	if(f(raiz)*f(m) > 0 || f(raiz)*f(m) == 0)  # verifica se a raíz está no intervalo [raiz, M]
	{
		m = raiz
	}

	i = i + 1
	return(raiz)
	}
}

#n=10
#Main(n)
a=1
b=1.5
erros_N(a,b)
