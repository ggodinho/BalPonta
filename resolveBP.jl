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
ncenNS = trunc(Int,serie)
ApontadorPatamarCarga = round.(Int, ApontadorPatamarCarga)
ApontadorIntercambio = round.(Int, ApontadorIntercambio)
UHEsubsis = round.(Int,UHEsubsis)
UHEsubsis = vec(UHEsubsis)
UHEsubsis = replace(UHEsubsis, 14=>6, 15=>7)
UTEsubsis = vec(round.(Int,UTE["subsistema"]))
ApontadorInterc_FIC_de = ApontadorIntercambio[nsis+1:nsis+nfic,1:nsis]
ApontadorInterc_FIC_para = ApontadorIntercambio[1:nsis,nsis+1:nsis+nfic]

CVU = UTE["CVU"]
LimInter = Intercambio["capacidade"]
nInter = size(Intercambio["origem"],2)

#Definindo penalidades:
PenInterc = 0.01;
PenDeficit = 3000;
PenSobra = 0.001;
#PenGH = [4 2 3 2 1 1 1 1 1 1 1]; # Prioriza a geração dos subsistemas a fio d'agua (Mad., B. Monte e T. Pires)
PenRG = maximum(CVU)*1.02;
PenRO = maximum(CVU)*1.01;

Ano_H = [1955; 2016]
nCenarioHidro = size(Ano_H,1)
Ano_h_count = Ano_H.-(1931+1)
serie_h = Ano_h_count[iCenH]

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
                serie_h=SerieHidro(iCenH,anoSimulacao);
                cenario=cenario+1;
                println("   Período $(iper) Cenario $(cenario)");

                #varre todos as horas
                for hora = 1:nHoras
                    #varre todos os cenarios de nao simuladas
                    # nesta etapa a Disponibilidade ainda esta agregada por subsistema, e
                    # sua informacao esta nas variaveis GHmax e GTmax
                    #dados de saida do PL

Inflex = UTE["inflex"][1:nusiTerm,iper][1:nusiTerm]
DispTerm = UTE["disp"][1:nusiTerm,iper][1:nusiTerm]
pat = ApontadorPatamarCarga[iper,hora]

mes = rem(iper,12);
if mes == 0
    mes = 12;
end

# Deficit(iCenNS,iCenH,:,iper,hora),
#     Excesso(iCenNS,iCenH,:,iper,hora),
#     Interc(iCenNS,iCenH,:,:,iper,hora),
#     ViolaRNE(iCenNS,iCenH,iper,hora),
#     GH(iCenNS,iCenH,:,iper,hora),
#     GHusina(iCenNS,iCenH,:,iper,hora),
#     ViolaReservaOP(iCenNS,iCenH,:,iper,hora),
#     ViolaReservaGer(iCenNS,iCenH,:,iper,hora),
#     GT(iCenNS,iCenH,:,iper,hora),
#     GTusina(iCenNS,iCenH,:,iper,hora)] =
    resolveBP(iper, hora, serie_h, pat, nsis, nfic, nusi, nusiTerm, nInter, Mercado[:,iper].*PUMercado[:,mes,hora],
        DispNS[icenNS,:,iper,hora], PenGH[serie_h,:,iper,pat], UHEsubsis, GHmax[serie_h,:,iper,pat],
        ROmax[:,iper,hora], RGmax[serie_h,:,iper,pat], GHmin[serie_h,:,iper,pat], CVU[:,iper],
        DispTerm, UTE["despacho"][serie_h,:,iper,pat], UTEsubsis, LimInter, RORNEmax[icenNS,iper,hora],
        numAgrint, Agrint, ApontadorIntercambio, ApontadorInterc_FIC_de, ApontadorInterc_FIC_para)

function resolveBP(iper, hora, serie_h, pat, nsis, nfic, nusi, nusiTerm, nInter, Mercado,
        DispNSim, PenGH, UHEsubsis, GHMax, ROmax, RGmax, GHMin, CVU, DispTerm, Inflex,
        UTEsubsis, LimInter, RORNEmax, numAgrint, Agrint, ApontadorIntercambio,
        ApontadorInterc_FIC_de, ApontadorInterc_FIC_para)

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

    #--------------------------------------------------------------------------
    # Montagem do problema
    #--------------------------------------------------------------------------
    #Escolhe o solver:
    BPontaProb = Model(with_optimizer(GLPK.Optimizer));

    #--------------------------------------------------------------------------
    # Variáveis do problema
    #--------------------------------------------------------------------------
    # Geração Hidráulica e Reserva de Geracao
    @variable(BPontaProb, ghid[1:nusi] >= 0);
    @constraint(BPontaProb,[i=1:nusi], ghid[i] >= GHmin[i])
    @constraint(BPontaProb,[i=1:nusi], ghid[i] <= GHmax[i])

    #Reserva de Geração
    @variable(BPontaProb, violaReservaGer[1:nusi] >= 0);
    @constraint(BPontaProb,[i=1:nusi], violaReservaGer[i] <= RGmax[i])

    # Reserva operativa
    @variable(BPontaProb, violaReservPotencia[1:nsis] >= 0);
    @constraint(BPontaProb, [i=1:nsis], violaReservPotencia[i] <= ROmax[i])

    # Geração Termoelétrica
    @variable(BPontaProb, gterm[1:nusiTerm]);
    @constraint(BPontaProb, [i=1:nusiTerm], gterm[i] >= Inflex[i])
    @constraint(BPontaProb, [i=1:nusiTerm], gterm[i] <= DispTerm[i])

    # Deficit
    @variable(BPontaProb, def[1:nsis] >= 0)

    # Sobras
    @variable(BPontaProb, sob[1:nsis] >= 0);

    # Intercambios
    @variable(BPontaProb, inter[1:nInter] >= 0)
    @constraint(BPontaProb, Rest_Inter[i=1:nInter], inter[i] <= LimInter[i][iper,hora])

    #--------------------------------------------------------------------------
    # Restrições
    #--------------------------------------------------------------------------

    # Atendimento ao mercado
    @constraint(BPontaProb, Atend_Demanda[i=1:nsis],
        sum(ghid[j] for j in findall(x -> x==i,UHEsubsis)) +
        sum(gterm[j] for j in findall(x -> x==i,UTEsubsis)) + def[i] - sob[i] +
        sum(inter[j] for j in ApontadorIntercambio[findall(x -> x>0, ApontadorIntercambio[:,i]),i]) -
        sum(inter[j] for j in ApontadorIntercambio[i,findall(x -> x>0, ApontadorIntercambio[i,:])])
        == Mercado[i] - DispNSim[i])

    #Restrição de Geração Hidráulica
    @constraint(BPontaProb, [i=1:nusi], ghid[i] - violaReservaGer[i]
        <= GHmax[i] - RGmax[i])

    #Restrição de Potência Operativa
    @constraint(BPontaProb, [i=1:nsis],
        sum(ghid[j] for j in findall(x -> x==i,UHEsubsis))
        - sum(violaReservaGer[j] for j in findall(x -> x==i,UHEsubsis))
        - violaReservPotencia[i]
        <= sum(GHmax[j] for j in findall(x -> x==i,UHEsubsis))
        - sum(RGmax[j] for j in findall(x -> x==i,UHEsubsis))
        - ROmax[i])

    #Restrição dos nós
    @constraint(BPontaProb, Rest_nosInter[i=1+nsis:nsis+nfic],
        sum(inter[j] for j in ApontadorIntercambio[findall(x -> x>0, ApontadorIntercambio[:,i]),i]) -
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
    @constraint(BPontaProb, RNE_const,
        sum(inter[j] for j in ApontadorIntercambio[findall(x -> x>0, AgrintAponta[1,:,:]),:])
        - sum(inter[j] for j in ApontadorIntercambio[findall(x -> x<0, AgrintAponta[1,:,:]),:])
        - violaRNE <= Agrint["Limite"][1][iper,pat] - RORNEmax)
    @constraint(BPontaProb, violaRNE <= RORNEmax)

    #Demais AGRINTS
    @constraint(BPontaProb, Agrint_const[i=1:numAgrint-1],
        sum(inter[j] for j in ApontadorIntercambio[findall(x -> x>0, AgrintAponta[i+1,:,:]),:]) -
        sum(inter[j] for j in ApontadorIntercambio[findall(x -> x<0, AgrintAponta[i+1,:,:]),:])
        <= Agrint["Limite"][i+1][iper,pat])

    #--------------------------------------------------------------------------
    # Função objetivo
    #--------------------------------------------------------------------------
    @objective(BPontaProb, Min,
        sum(def[j] for j in 1:nsis)*PenDeficit
        + PenSobra*sum(sob[j] for j in 1:nsis)
        + sum(PenGH[i]*sum(ghid[k] for k in findall(x -> x==i,UHEsubsis)) for i in 1:nsis)
        + PenRO*sum(violaReservPotencia[j] for j in 1:nsis)
        + PenRG*sum(violaReservaGer[k] for k in 1:nusi)
        + sum(CVU[l]*gterm[l] for l in 1:nusiTerm)
        + PenRO*violaRNE
        + PenInterc*sum(inter[j] for j in ApontadorIntercambio[findall(x -> x>0, ApontadorIntercambio[1:nsis,1:nsis])])
        + (PenInterc/2)*(sum(inter[j] for j in ApontadorInterc_FIC_de[findall(x -> x>0, ApontadorInterc_FIC_de[:,:])])
            + sum(inter[j] for j in ApontadorInterc_FIC_para[findall(x -> x>0, ApontadorInterc_FIC_para[:,:])])))

    #--------------------------------------------------------------------------
    # Solução do problema
    #--------------------------------------------------------------------------
    optimize!(BPontaProb)
    status = termination_status(BPontaProb)
    solution_obj = objective_value(BPontaProb)
    solution_def = value.(def)
    solution_exc = value.(sob)
    solution_interc = value.(inter)
    solution_violaRNE = value.(violaRNE)
    solution_GH = zeros(nsis)
    for i=1:nsis
        solution_GH[i] = sum(value.(ghid)[j] for j in findall(x -> x==i,UHEsubsis))
    end
    solution_GHusina = value.(ghid)
    solution_violaReservaOP = value.(violaReservPotencia)
    solution_violaReservaG = value.(violaReservaGer)
    solution_GT = zeros(nsis)
    for i=1:maximum(UTEsubsis)
        solution_GT[i] = sum(value.(gterm)[j] for j in findall(x -> x==i,UTEsubsis))
    end
    solution_GTusina = value.(gterm)

    if status == MOI.OPTIMAL
        println("Problema convergido!")
        return solution_def, solution_exc, solution_interc, solution_violaRNE,
            solution_GH, solution_GHusina, solution_violaReservaOP,
            solution_violaReservaG, solution_GT, solution_GTusina, solution_obj
    else
        println("Problema não convergido. Não foi encontrada uma solução ótima.")
        #break
    end

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


 # #---> Demais Restricoes não implementadas
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
