####################################
##
## Simulando dados com a função grf()
## (simulação não-condicional)
##
#####################################
##
## 1. Preparativos iniciais
## ------------------------
##
##
## 1.1 Carregando o pacote geoR
##
## se necessário atualize a versão do geoR com o comando:
#install.packages("geoR", contrib="http://www.est.ufpr.br/~paulojus/geoR/windows")
#
require(geoR)
##
## 1.2 Iniciando a ajuda no formato html (opcional)
##
#options(htmlhelp=TRUE)
#help.start()
##
## 1.3 Salvando padrões gráficos iniciais e inicializando o gerador de números aleatórios 
##
par.ori <- par(no.readonly = TRUE)
if(!exists(".Random.seed")) set.seed(123)
##
## 1.4 Inspecionando os argumentos da função grf()
##
args(grf)
#help(grf)
##
## 2. Uma primeira simulação com argumentos mínimos usando diversas opções "default"
## ---------------------------------------------------------------------------------
##
## 2.1 simula dados:
##
sim01 <- grf(50, cov.pars=c(1, .25))
##
## 2.2 inspeciona objeto com dados simulados:
##
sim01
##
## 2.3 inspeciona variogramas teórico e empírico
##
plot(sim01)
##
## 2.4 visualizando as localizações simuladas
##
## pode-se ainda visualizar as localizações dos dados simulados
par(mfrow=c(1,2))
plot(sim01, plot.loc=T)
par(mfrow=c(1,1))
## ou usar funções padrao de visualização de dados
plot.geodata(sim01)
points.geodata(sim01)
##
## Nota: veja ainda opções destas funções com:
##        args(plot.geodata)
##        args(points.geodata)
##
##
## 3. Uma segunda simulação em grid regular
## ----------------------------------------
##
## 3.1 Simulando dados
##
sim02 <- grf(256, cov.pars=c(1, .3), grid="reg") 
##
## 3.2 Visualizando as simulações em grid regular
##
image(sim02)
persp(sim02)
## e, se desejar mude o angulo de visualização 
persp(sim02, theta=30, phi=20)
##
## Exercicio 01: Explorando as opções graficas
##
##  - inspecione outras opções em help(image.grf) e help(persp.grf)
##    e use estas opções para midoficar os graficos
##
##  - use as outras funções de visualização mencionadas anteriormente
##    (plot, plot.geodata, points.geodata)
##
##
## Exercicio 02: Explorando outros argumentos da função grf()
##
##  - gere novas simulações usando os argumentos: nx, ny, xlims e ylims
##
##
## 4. Usando simulações em 1D para compreender o modelo e seus parametros
## ----------------------------------------------------------------
##
## 4.1 Simulando com diferentes alcances ("range") definidos pelo parametro $\phi$ 
##
sim11 <- grf(100, ylims=c(0,0), cov.pars=c(1, 0.25))
image(sim11)
image(sim11, type="l")
image(sim11, type="b", cex=0.5)
##
sim12 <- grf(100, ylims=c(0,0), cov.pars=c(1, 0.0))
image(sim12, type="l")
##
sim13 <- grf(100, ylims=c(0,0), cov.pars=c(1, 0.5))
image(sim13, type="l")
##
par(mfrow=c(3,1), mar=c(2,2,1,1))
yl <- range(c(sim11$data, sim12$data, sim13$data))
image(sim11, type="l", ylim=yl)
image(sim12, type="l", ylim=yl)
image(sim13, type="l", ylim=yl)
par(par.ori)
##
## Exercicio 03: 
##
##  - gere diferentes simulações variando o parametro $\sigma^2$ no argumento cov.pars
##    e discuta os resultados
##
##
## 4.2 Simulações com diferentes modelos de variograma
##
sim21 <- grf(100, ylims=c(0,0), cov.pars=c(1, 0.75), cov.model="sph")
image(sim21, type="l")
##
sim22 <- grf(50, ylims=c(0,0), cov.pars=c(1, 0.75/sqrt(3)), cov.model="mat", kappa=2)
image(sim22, type="l")
##
par(mfrow=c(3,1), mar=c(2,2,1,1))
yl <- range(c(sim11$data, sim21$data, sim22$data))
image(sim11, type="l", ylim=yl)
image(sim21, type="l", ylim=yl)
image(sim22, type="l", ylim=yl)
par(par.ori)
##
##
## Exercicio 04: Explorando os modelos de variograma
##
##  - gere e compare simulações com o modelo de Matern
##    (cov.model="mat") com diferentes valores
##    no argumento kappa (tente, por exemplo valores 0.5, 2 e 5)
##
##  - tente tambem com cov.model="gau" e observe se ha'
##    mensagens de alerta
##
##  - gere e compare simulações com outros modelos de variogramas
##
##
## 4.3 Simulações com ruido
##
sim31 <- grf(50, ylims=c(0,0), cov.pars=c(1, 0.75/sqrt(3)), nugget = 0.1, cov.model="gau")
##
## repare novamente se ha' alguma mensagem de alerta (warning) !
##
image(sim31, type="l")
##
## Agora gerando simulações com a mesma "semente" e com
## diferentes  niveis de ruido
##
set.seed(12)
sim32 <- grf(100, ylims=c(0,0), cov.pars=c(1, 0.25))
##
set.seed(12)
sim33 <- grf(100, ylims=c(0,0), cov.pars=c(1, 0.25), nug=0.25)
##
set.seed(12)
sim34 <- grf(100, ylims=c(0,0), cov.pars=c(1, 0.25), nug=0.5)
##
par(mfrow=c(3,1), mar=c(2,2,1,1))
yl <- range(c(sim32$data, sim33$data, sim34$data))
image(sim32, type="l", ylim=yl)
image(sim33, type="l", ylim=yl)
image(sim34, type="l", ylim=yl)
par(par.ori)
##
##
## Exercicio 05: Simulações com ruido
##
##  - gere em compare mais simulações com differentes valores no argumento nugget
##
##
## 5. Modelos com anisotropia
## --------------------------
##
## simulação sem estrutura espacial
sim39 <- grf(961, grid="reg", cov.pars=c(1, 0))
image(sim39)
## simulação com estrutura espacial isotropica
sim40 <- grf(961, grid="reg", cov.pars=c(1, .15))
image(sim40)
## simulação com estrutura espacial anisotropica
sim41 <- grf(961, grid="reg", cov.pars=c(1, .15), aniso.pars=c(pi/4, 3))
image(sim41)
##
## Exercicio 06: Simulações com anisotropia
##
##  - gere em compare simulações com differentes valores no argumento aniso.pars
##
##
## 6. Simulando modelos com tendencia
## ----------------------------------
##  A função grf gera simulações com media constante e igual a zero.
##  Para gerar simulações com media diferente de zero basta somar
##  quantidades desejadas aos dados simulados como no exemplo a seguir
##
sim51 <- grf(100, ylims=c(0,0), cov.pars=c(1, 0.25))
sim51$data <- 10 - 5*sim51$coords[,1] + sim51$data
image(sim51, ty="l")
##
##
## 7. Simulações de variaveis nao-Gaussianas
##
## Simulação de variavel log-Normal
##
sim61 <- grf(961, grid="reg", cov.pars=c(1, 0.25), lambda = 0)
image(sim61)
hist(sim61$data)
##
## Simulação de variavel Poisson
##
sim71 <- grf(121, grid="reg", cov.pars=c(1, 0.25))
sim71$data <- rpois(121, exp(0.5 + sim71$data))
image(sim71, col=gray(seq(1,0.2,l=11)))
text(sim71$coords[,1], sim71$coords[,2], sim71$data)
barplot(table(sim71$data))
##
## Exercicio 07: Simulações de variaveis nao-Gaussianas
##
##  - gere, explore e compare simulações com differentes valores no argumento lambda
##  - gere simulações com variavel resposta Binomial
##
##
## 8. Gerando simulações em grids muito finos (malha com muitos pontos)
##
sim81 <- grf(40401, grid = "reg", cov.pars = c(10, .2), met = "circ") 
image(sim81)
##
## =======================================================
## Para mais sobre simulações veja o pacote RandomFields.
## =======================================================
##


