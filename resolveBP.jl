#Pkg.add("MATLAB")
#Pkg.add("MAT")

using MATLAB, MAT
using JuMP,GLPK
using SparseArrays
#Pkg.add("GLPK")

#------------------------------------------------------------------------------#
#Função para extrair variáveis de um DIct
function extract(d)
           expr = quote end
           for (k, v) in d
               push!(expr.args, :($(Symbol(k)) = $v))
           end
           eval(expr)
           return
end

#Função para mostrar todo ranges no console
#Ex: myshowall(stdout, ReservaGeracao[1,1,iper,pat], false)
function myshowall(io, x, limit = false)
  println(io, summary(x), ":")
  Base.print_matrix(IOContext(io, :limit => limit), x)
end
#------------------------------------------------------------------------------#
#Extraindo variáveis MATLAB:
Path = dirname(Base.source_path());
matvariables = matread(Path * "\\test2.mat")
extract(matvariables)

#Transformando Arrays Float em Int
nusi = trunc(Int,nusi)
nsis = trunc(Int,nsis)
nper = trunc(Int,nper)
nHoras = trunc(Int,nHoras)
nCenarios = trunc(Int,nCenariosNSimul)
nfic = trunc(Int,nfic)
nusiTerm = trunc(Int,nusiTerm)
nCenariosNSimul = trunc(Int,nCenariosNSimul)
nCenarioHidro = trunc(Int,nCenarioHidro)
numAgrint = trunc(Int,numAgrint)
serie = trunc(Int,serie)
ApontadorPatamarCarga = round.(Int, ApontadorPatamarCarga)
ApontadorIntercambio = round.(Int, ApontadorIntercambio)
UHEsubsis = round.(Int,UHEsubsis)
UHEsubsis = vec(UHEsubsis)
UHEsubsis = replace(UHEsubsis, 14=>6, 15=>7)
UTEsubsis = vec(round.(Int,UTE["subsistema"]))

CVU = UTE["CVU"]
nInter = size(Intercambio["origem"],2)

#Definindo penalidades:
PenInterc = 0.01;
PenDeficit = 3000;
PenSobra = 0.001;
#PenGH = [4 2 3 2 1 1 1 1 1 1 1]; # Prioriza a geração dos subsistemas a fio d'agua (Mad., B. Monte e T. Pires)
PenGT = 10;
PenRG = maximum(CVU)*1.02;
PenRO = maximum(CVU)*1.01;


#Deficit = Float64[]

#Criando matrizes de saída
Deficit = zeros(nCenariosNSimul,nCenarioHidro,nsis,nper,nHoras);
Excesso = zeros(nCenariosNSimul,nCenarioHidro,nsis,nper,nHoras);
Interc = zeros(nCenariosNSimul,nCenarioHidro,nsis+nfic,nsis+nfic,nper,nHoras);
ViolaRNE = zeros(nCenariosNSimul,nCenarioHidro,nper,nHoras);
GH = zeros(nCenariosNSimul,nCenarioHidro,nsis,nper,nHoras);
GHusina = zeros(nCenariosNSimul,nCenarioHidro,nusi,nper,nHoras);
ViolaReservaOP = zeros(nCenariosNSimul,nCenarioHidro,nsis,nper,nHoras);
ViolaReservaGer = zeros(nCenariosNSimul,nCenarioHidro,nusi,nper,nHoras);
GT = zeros(nCenariosNSimul,nCenarioHidro,nsis,nper,nHoras);
GTusina = zeros(nCenariosNSimul,nCenarioHidro,nusiTerm,nper,nHoras);
Sobra = zeros(nCenariosNSimul,nCenarioHidro,nsis,nper,nHoras);           # Inicializa a matriz de sobras
SobraMin = ones(nCenariosNSimul,nCenarioHidro,nsis,nper);                # Inicializa a matriz de sobras
SobraMinSin = ones(nCenariosNSimul,nCenarioHidro,nper);                  # Inicializa a matriz de sobras
println("Calculando balanço de ponta determinístico...")

#varre todos os meses
for iper = mesi:nper
    #obtem o mes de simulação
    mes=rem(iper,12);
    if mes == 0
        mes = 12;
    end
    anoSimulacao = floor((iper-1)/12)+1;
	#varre todos os cenarios de nao simuladas
	cenario=0;
    for iCenNS = 1:nCenariosNSimul
        #se a probabilidade do cenario for 0, vai para o próximo
        if probCenarioNS(iCenNS,iper) > 0
            #varre todos os cenarios hidrologicos
            for iCenH = 1:nCenarioHidro
                iSerieH=SerieHidro(iCenH,anoSimulacao);
                cenario=cenario+1;
                println("   Período $(iper) Cenario $(cenario)");

                #varre todos as horas
                for ihora = 1:nHoras
                    #varre todos os cenarios de nao simuladas
                    # nesta etapa a Disponibilidade ainda esta agregada por subsistema, e
                    # sua informacao esta nas variaveis GHmax e GTmax
                    #dados de saida do PL
                    [Deficit(iCenNS,iCenH,:,iper,ihora), ...
                     Excesso(iCenNS,iCenH,:,iper,ihora), ...
                     Interc(iCenNS,iCenH,:,:,iper,ihora), ...
                     ViolaRNE(iCenNS,iCenH,iper,ihora), ...
                     GH(iCenNS,iCenH,:,iper,ihora), ...
                     GHusina(iCenNS,iCenH,:,iper,ihora), ...
                     ViolaReservaOP(iCenNS,iCenH,:,iper,ihora), ...
                     ViolaReservaGer(iCenNS,iCenH,:,iper,ihora), ...
                     GT(iCenNS,iCenH,:,iper,ihora), ...
                     GTusina(iCenNS,iCenH,:,iper,ihora)] = ... #dados de entrada do PL
                      resolveBP(iper, ihora, ApontadorPatamarCarga(mes,ihora), nsis, nfic, nusi, nusiTerm, ...
                          ApontadorSistema, Sistema, Mercado(:,iper).*PUMercado(:,mes,ihora), ...
                          DispNS(iCenNS,:,iper,ihora), PenGH(iSerieH,:,iper,ApontadorPatamarCarga(mes,ihora)), ...
                          UHEsubsis, GHmax(iSerieH,:,iper,ApontadorPatamarCarga(mes,ihora)), ROmax(:,iper,ihora), ...
                          RGmax(iSerieH,:,iper,ApontadorPatamarCarga(mes,ihora)), ...
                          GHmin(iSerieH,:,iper,ApontadorPatamarCarga(mes,ihora)), ...
                          UTE.CVU(:,iper), UTE.disp(:,iper), UTE.despacho(iSerieH,:,iper,ApontadorPatamarCarga(mes,ihora)), ...
                          UTE.subsistema, LimInterc(:,iper, ApontadorPatamarCarga(mes,ihora)), RORNEmax(iCenNS,iper,ihora), ...
                          numAgrint, Agrint, Restricao, ApontadorConfHd, ApontadorConfT);


#using Cbc # Open source solver. Must support integer programming.


# x = range(-10.0, stop=10.0, length=500)
# mat"plot($x, sin($x))"  # evaluate a MATLAB function
#
# y = range(2.0, stop=3.0, length=500)
# mat"""
#     $u = $x + $y
# 	$v = $x - $y
# """
# @show u v



# files = ["Deficit", "Excesso", "Interc","GH","GT","GTusina","Sobra"]
# saidas = Dict()
# for (n, f) in enumerate(files)
#     if f == "Interc"
#         saidas[f] = Dict("DE->PARA" =>Dict("mes" => Dict("hora" => Dict("scen" => []))), "PARA->DE" =>Dict("mes" => Dict("hora" => Dict("scen" => []))))
#     elseif f == "GTusina"
#         saidas[f] = Dict("ute" => Dict("mes" => Dict("hora" => Dict("scen" => []))))
#     else
#         saidas[f] = Dict("sistema" => Dict("mes" => Dict("hora" => Dict("scen" => []))))
#     end
# end

#Convertendo de Float para Int


# Deficit = zeros(nsis,nper,nHoras,nCenarios,nsim);
# Excesso = zeros(nsis,nper,nHoras,nCenarios,nsim);
# Interc = zeros(nsis+nfic,nsis+nfic,nper,nHoras,nCenarios,nsim);
# GH = zeros(nsis,nper,nHoras,nCenarios,nsim);
# GT = zeros(nsis,nper,nHoras,nCenarios,nsim);
# GTusina = zeros(nusiTerm,nper,nHoras,nCenarios,nsim);
# #fprint('Calculando balanço de ponta determinístico...')
# Sobra = zeros(nsis,nper,nHoras,nCenarios);

# LimInterc = Intercambio["capacidade"][:][mes,ApontadorPatamarCarga[mes,iper]]=
iper = 5
hora = 1
Inflex = UTE["inflex"][1:nusiTerm,iper][1:nusiTerm]
DispTerm = UTE["disp"][1:nusiTerm,iper][1:nusiTerm]
DispNSim = DispNS[serie,:,iper,hora]
pat = ApontadorPatamarCarga[iper,hora]


function resolveBP(iper, ihora, PatamarCarga, nsis, nfic, nUsiHidro, nUsiTerm, ApontadorSistema, Sistema, Mercado, DispNS, PenGH, ...
        #UHEsubsis, DispHidro, ReservaOperativa, ReservaGeracao, GHMin, CVU, DispTerm, Inflex, ...
        #UTEsis, LimInterc, ReservaOperativaRNE, numAgrint, Agrint, Restricao, ApontadorConfHd, ApontadorConfT)

# # Resolve o seguinte problema:
# # min sum(def) + PenInterc * sum(interc) + PenSobra * sum(sobra)
# #     + PenGH*sum(GHusina) + PenReservPotencia*sum(violaReservaOp)
# #     + PenReservGeracao*sum(violaResevaGer) + PenGT*sum(GT)
# #     + PenReservPotencia*sum(violaRNE)
# # s.a.
# # sum(GHusina) + sum(GTusina) + Deficit - Sobra + Imp - Exp = Mercado - DispNS
# # Imp < LimInterc
# # Exp < LimInterc
# # GHmin < GH - violaResevaGer < GHmax - ReservaGermax
# # sum(GH) - sum(violaResevaGer) - violaReservPotencia < sum(GHmax) - sum(ReservaGermax) - PenReservPotenciamax
# # violaReservaOp < ReservaOPmax
# # violaReservaGer < ReservaGermax
# # Inflex < GT < DispTerm
# # sum(interc) < Agrint
# # RNE - violaRNE < RNEmax - RORNEmax
# # violaRNE < RORNEmax

    #---> Declaração de variáveis


#IntercAponta = zeros(size(LimInterc,1),size(LimInterc,1));

#--------------------------------------------------------------------------
# Montagem do problema
#--------------------------------------------------------------------------
#Escolhe o solver:
BPontaProb = Model(with_optimizer(GLPK.Optimizer));

#---> Variáveis do problema
# Geração Hidroeletrica
# ghid_ub = zeros(nsis)
# reservaP_ub = zeros(nsis)
# for isis=1:nsis
#     if isis == 3 || ReservaIntercambioNE == 1
#         #se a reserva de potencia do nordeste estiver no intercambio nao
#         #retira da hidraulica
#         ghid_ub[isis] = DispHidro[isis];
#         reservaP_ub[isis] = 0;
#         # Limite Superior
#         else
#         ghid_ub[isis] = DispHidro[isis] - Mercado[isis]*reserva/100;
#         reservaP_ub[isis] = Mercado[isis]*reserva/100;   # Limite Superior
#     end
# end

# Geração Hidráulica e Reserva de Geracao
@variable(BPontaProb, ghid[1:nusi] >= 0);
@variable(BPontaProb, violaReservaGer[1:nusi] >= 0);

if ReservaIntercambioNE[iper] == 1
    # se a reserva de potencia do nordeste estiver no intercambio
    # nao retira da hidraulica
    @constraint(BPontaProb,[i=1:nsis], ghid[i] >= GHmin[1,i,iper,pat])
    @constraint(BPontaProb,[i=1:nsis], ghid[i] <= GHmax[1,i,iper,pat])
else
    @constraint(BPontaProb,[i=1:nusi], ghid[i] - violaReservaGer[i]>= GHmin[1,i,iper,pat])
    @constraint(BPontaProb,[i=1:nusi], ghid[i] - violaReservaGer[i]<= GHmax[1,i,iper,pat]
        - ReservaGeracao[1,i,iper,pat])
    @constraint(BPontaProb, [i=1:nusi], violaReservaGer[i] <= ReservaGeracao[1,i,iper,pat])
end


# Reserva operativa
@variable(BPontaProb, violaReservPotencia[1:nsis] >= 0);
@constraint(BPontaProb, [i=1:nsis], violaReservPotencia[i] <= ROmax[i,iper,hora])

# Geração Termoelétrica
@variable(BPontaProb, gterm[1:nusiTerm] >= 0);
@constraint(BPontaProb, [i=1:nusiTerm], gterm[i] >= Inflex[i])
@constraint(BPontaProb, [i=1:nusiTerm], gterm[i] <= DispTerm[i])

# Deficit
@variable(BPontaProb, def[1:nsis] >= 0)

# Sobras
@variable(BPontaProb, sob[1:nsis] >= 0);

# Intercambios
#Falta declarar penalidades dos intercâmbios e apontador adequadamente...

#Intercambio
#inter_ub = zeros(size(Intercambio["origem"],2))
#for i=1:size(Intercambio["origem"],2)
        #inter_ub[i] = Intercambio["capacidade"][:][i][5,4]; #[5,4 = mes e patamar], serão variaveis de entrada
#end

Liminter = zeros(nInter)
for i=1:nInter
    Liminter[i] = Intercambio["capacidade"][i][iper,hora]
end

@variable(BPontaProb, inter[1:nInter] >= 0)

#---> Restrições do problema
@constraint(BPontaProb, Rest_Inter[i=1:nInter], inter[i] <= Liminter[i])

# Atendimento ao mercado
@constraint(BPontaProb, Atend_Demanda[i=1:nsis],sum(ghid[j] for j in findall(x -> x==i,UHEsubsis)) +
    sum(gterm[j] for j in findall(x -> x==i,UTEsubsis)) + def[i] - sob[i] +
    sum(inter[j] for j in ApontadorIntercambio[findall(x -> x>0, ApontadorIntercambio[:,i]),i]) -
    sum(inter[j] for j in ApontadorIntercambio[i,findall(x -> x>0, ApontadorIntercambio[i,:])])
    == Mercado[i,iper] - DispNSim[i])

@constraint(BPontaProb, [i=1:nsis], sum(ghid[j] for j in findall(x -> x==i,UHEsubsis)) -
    sum(violaReservaGer[j] for j in findall(x -> x==i,UHEsubsis)) - violaReservPotencia[i]
    <= sum(GHmax[1,j,iper,pat] for j in findall(x -> x==i,UHEsubsis)) -
    ReservaGeracao_sist[1,i,iper,pat] - ROmax[i,iper,hora])

#Restrição dos nós
@constraint(BPontaProb, Rest_nosInter[i=1+nsis:nsis+nfic], sum(inter[j] for j in ApontadorIntercambio[findall(x -> x>0, ApontadorIntercambio[:,i]),i]) -
    sum(inter[j] for j in ApontadorIntercambio[i,findall(x -> x>0, ApontadorIntercambio[i,:])]) == 0)

#---> Restrições de Agrupamento de Intercambio
Agrint_reg = vec(round.(Int,Agrint["nReg"][:]))
AgrintAponta = zeros(size(Agrint["nReg"],2),nsis+nfic,nsis+nfic)

#Criando a matriz DE-PARA dos AGRINTS
for iagr = 1:round(Int,numAgrint)
    for j = 1:Agrint_reg[iagr]
       isis = trunc(Int64,Agrint["De"][iagr][j]);
       jsis = trunc(Int64,Agrint["Para"][iagr][j]);
       AgrintAponta[iagr,isis,jsis] = iagr*trunc(Int64,Agrint["Fator"][iagr][j])
    end
end

# Reserva de potencia no intercambio do NE
@variable(BPontaProb, violaRNE >= 0);
@constraint(BPontaProb, RNE_const,  sum(inter[j] for j in ApontadorIntercambio[findall(x -> x>0, AgrintAponta[1,:,:]),:]) -
    sum(inter[j] for j in ApontadorIntercambio[findall(x -> x<0, AgrintAponta[1,:,:]),:])
    - violaRNE <= Agrint["Limite"][1][iper,pat] - RORNEmax[serie,iper,hora])
@constraint(BPontaProb, violaRNE <= RORNEmax[serie,iper,hora])

#Demais AGRINTS
@constraint(BPontaProb, Agrint_const[i=1:numAgrint-1],  sum(inter[j] for j in ApontadorIntercambio[findall(x -> x>0, AgrintAponta[i+1,:,:]),:]) -
    sum(inter[j] for j in ApontadorIntercambio[findall(x -> x<0, AgrintAponta[i+1,:,:]),:]) <= Agrint["Limite"][i+1][iper,pat])


# # s.a.
# # sum(GHusina) + sum(GTusina) + Deficit - Sobra + Imp - Exp = Mercado - DispNS
# # Imp < LimInterc
# # Exp < LimInterc
# # GHmin < GH - violaResevaGer < GHmax - ReservaGermax
# # sum(GH) - sum(violaResevaGer) - violaReservPotencia < sum(GHmax) - sum(ReservaGermax) - PenReservPotenciamax
# # violaReservaOp < ReservaOPmax
# # violaReservaGer < ReservaGermax
# # Inflex < GT < DispTerm
# # sum(interc) < Agrint
# # RNE - violaRNE < RNEmax - RORNEmax
# # violaRNE < RORNEmax


# irest = 0;
#
# nrest = irest;             # Número de restrições
#
# A = zeros(nUsiHidro+nsis+numAgrint,nvar);
#
#
# #---> Demais Restricoes
# nRestricoes = 0;
# for iRest = 1:size(Restricao,2)
#     #verifica se deve incluir a restricao
#     incluir = 0;
#     #se for genérico, inclui a restricao
#     if strcmp(Restricao(iRest).VarContr,'NULL') == 1
#         incluir = 1;
#     #se for incluir dependendo do valor da carga
#     elseif strcmp(Restricao(iRest).VarContr,'CARGA') == 1
#         #verifica de qual subsistema
#         if Restricao(iRest).VarContrSubsis ~= 0
#             carga = Mercado(ApontadorSistema(Restricao(iRest).VarContrSubsis));
#         else
#             #calcula a carga do SIN
#             carga = sum(Mercado);
#         end
#         #verifica se deve incluir
#         if ((Restricao(iRest).VarContrValor(1) == -1 || carga <= Restricao(iRest).VarContrValor(1)) && ...
#             (Restricao(iRest).VarContrValor(2) == -1 || carga >= Restricao(iRest).VarContrValor(2)))
#             incluir = 1;
#         end
#     end
#     #inclui a restricao no PL
#     if incluir == 1
#         #verifica se possui limite Maximo
#         if Restricao(iRest).LimiteMax(iper,ihora) ~= 999999999
#             nRestricoes = nRestricoes + 1;
#             #varre todas as UHEs envolvidas na restricao
#             for indice = 1: size(Restricao(iRest).UHE,2)
#                 A(nUsiHidro+nsis+numAgrint+nRestricoes,ApontadorConfHd(Restricao(iRest).UHE(indice))) = ...
#                     Restricao(iRest).FatorUHE(indice);
#             end
#             #varre todas as UTEs envolvidas na restricao
#             for indice = 1: size(Restricao(iRest).UTE,2)
#                 A(nUsiHidro+nsis+numAgrint+nRestricoes,2*nUsiHidro+nsis+ApontadorConfT(Restricao(iRest).UTE(indice))) = ...
#                     Restricao(iRest).FatorUTE(indice);
#             end
#             #varre todas as interligações envolvidas na restricao
#             for indice = 1: size(Restricao(iRest).DE,2)
#                 A(nUsiHidro+nsis+numAgrint+nRestricoes,IntercAponta(Restricao(iRest).DE(indice),Restricao(iRest).PARA(indice))) = ...
#                     Restricao(iRest).FatorInter(indice);
#             end
#             #verifica se é uma restricao que deve incluir a reserva no
#             #limite de intercambio
#             if Restricao(iRest).Codigo<20
#                 A(nUsiHidro+nsis+numAgrint+nRestricoes,ApontadorReservaRNE) = -1;
#                 b(nUsiHidro+nsis+numAgrint+nRestricoes,1) = Restricao(iRest).LimiteMax(iper,ihora) - ...
#                     ReservaOperativaRNE;
#             else
#                 b(nUsiHidro+nsis+numAgrint+nRestricoes,1) = Restricao(iRest).LimiteMax(iper,ihora);
#             end
#         end
#         #verifica se possui limite Minimo
#         if Restricao(iRest).LimiteMin(iper,ihora) ~= 999999999
#             nRestricoes = nRestricoes + 1;
#             #varre todas as UHEs envolvidas na restricao
#             for indice = 1: size(Restricao(iRest).UHE,2)
#                 A(nUsiHidro+nsis+numAgrint+nRestricoes,ApontadorConfHd(Restricao(iRest).UHE(indice))) = ...
#                     -1*Restricao(iRest).FatorUHE(indice);
#             end
#             #varre todas as UTEs envolvidas na restricao
#             for indice = 1: size(Restricao(iRest).UTE,2)
#                 A(nUsiHidro+nsis+numAgrint+nRestricoes,2*nUsiHidro+nsis+ApontadorConfT(Restricao(iRest).UTE(indice))) = ...
#                     -1*Restricao(iRest).FatorUTE(indice);
#             end
#             #varre todas as interligações envolvidas na restricao
#             for indice = 1: size(Restricao(iRest).DE,2)
#                 A(nUsiHidro+nsis+numAgrint+nRestricoes,IntercAponta(Restricao(iRest).DE(indice),Restricao(iRest).PARA(indice))) = ...
#                     -1*Restricao(iRest).FatorInter(indice);
#             end
#             b(nUsiHidro+nsis+numAgrint+nRestricoes,1) = -1*Restricao(iRest).LimiteMin(iper,ihora);
#         end
#     end
# end
#    xlswrite('matrizA', A(:,1:250),'1','A1:IP320');
#    xlswrite('matrizA', A(:,251:size(A,2)),'2','A1:IP320');
#  xlswrite('matrizb', b);
#--------------------------------------------------------------------------
# Solução do problema
#--------------------------------------------------------------------------
@objective(BPontaProb,Min,sum(def[j] for j in 1:nsis)*PenDeficit + PenSobra*sum(sob[j] for j in 1:nsis) +
    sum(PenGH[serie_h,i,iper,pat]*sum(ghid[k] for k in findall(x -> x==i,UHEsubsis)) for i in 1:nsis)
    + PenRO*sum(violaReservPotencia[j] for j in 1:nsis) +
    PenRG*sum(violaReservaGer[k] for k in 1:nusi)+ PenGT*sum(gterm[l] for l in 1:nusiTerm) + PenRO*violaRNE)

status = optimize!(BPontaProb)
termination_status(BPontaProb)
FALTA PENALIZAR INTERCAMBIO
serie_h = 1
# # min sum(def) + PenInterc * sum(interc) + PenSobra * sum(sobra)
# #     + PenGH*sum(GHusina) + PenReservPotencia*sum(violaReservaOp)
# #     + PenReservGeracao*sum(violaResevaGer) + PenGT*sum(GT)
# #     + PenReservPotencia*sum(violaRNE)

[x,fval,exitflag,output] = linprog(f,A,b,Aeq,beq,lb,ub,x0,opcoes);
if (exitflag ~= 1)
   disp('')
   disp(output)
   error('Não foi encontrada uma solução ótima.')
end

#--------------------------------------------------------------------------
# Carrega a solução do problema
#--------------------------------------------------------------------------
ivar = 0;
# Geração Hidroelétrica
GH = zeros(nsis,1);
GHusina = zeros(nUsiHidro,1);
for iusina = 1:nUsiHidro
   ivar = ivar + 1;
   GHusina(iusina) = x(ivar);
   GH(ApontadorSistema(UHEsubsis(iusina))) = GH(ApontadorSistema(UHEsubsis(iusina))) + x(ivar);
end
# Reserva Geracao
ViolaReservaGer = zeros(nUsiHidro,1);
for iusina = 1:nUsiHidro
   ivar = ivar + 1;
   ViolaReservaGer(iusina) = x(ivar);
end
# Reserva Operativa
ViolaReservaOP = zeros(nsis,1);
for isis = 1:nsis
   ivar = ivar + 1;
   ViolaReservaOP(isis) = x(ivar);
end
# Geração Termoelétrica
GT = zeros(nsis,1);
GTusina = zeros(size(UTEsis,2),1);
for iusina = 1:size(UTEsis,2)
   ivar = ivar + 1;
   GTusina(iusina,1) = x(ivar);
   GT(ApontadorSistema(UTEsis(iusina)),1) = GT(ApontadorSistema(UTEsis(iusina)),1) + x(ivar);
end
# Deficit
Deficit = zeros(nsis,1);
for isis = 1:nsis
   ivar = ivar + 1;
   Deficit(isis,1) = round(x(ivar));
end
# Sobras
Sobra = zeros(nsis,1);
for isis = 1:nsis
   ivar = ivar + 1;
   Sobra(isis,1) = x(ivar);
end
# Intercambios
Interc = zeros(nsis+nfic,isis-1);
ivar1 = 0;
inicio = ivar+1;
 for isis = 1:nsis+nfic
    for jsis = 1:isis-1
       ivar1 = ivar1 + 1;
       #só inclui no PL as linhas existentes
       if LimInterc(ivar1) > 0
           ivar = ivar + 1;
           Interc(isis,jsis) = x(ivar);
       end
    end
    for jsis = isis+1:nsis+nfic
       ivar1 = ivar1 + 1;
       #só inclui no PL as linhas existentes
       if LimInterc(ivar1) > 0
           ivar = ivar + 1;
           Interc(isis,jsis) = x(ivar);
       end
    end
 end
 # subtrai o intercambio x1/x2 do x2/x1 para evitar fluxo nos dois sentido
 for isis = 1:nsis+nfic
     for jsis = isis+1:nsis+nfic
         if(Interc(isis,jsis)>=Interc(jsis,isis))
             Interc(isis,jsis) = Interc(isis,jsis)-Interc(jsis,isis);
             Interc(jsis,isis) = 0;
         else
             Interc(jsis,isis) = Interc(jsis,isis)-Interc(isis,jsis);
             Interc(isis,jsis) = 0;
         end
     end
 end
 fim = ivar;
# GHusina(32)+GHusina(33)
# Interc(11,12)
# Interc(13,1)
# xlswrite('GHusina', GHusina);
# xlswrite('Interc', x(inicio:fim));
# pause
 ivar = ivar + 1;
 ViolaRNE = x(ivar);

 #end
