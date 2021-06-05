%cls :- write('\33\[2J'). -> para limpar o terminal  

%Declaracoes iniciaias

:- set_prolog_flag(discontiguous_warnings,off).
:- set_prolog_flag(single_var_warnings,off).

%Definicoes iniciais
:- op( 900,xfy,'::' ).
:- use_module(library(lists)).

:- include('pontos_recolha.pl').
:- include('arcos.pl').

% --------------------------------------------------- Predicados das Aulas -----------------------------------------------------------

nao( Questao ) :-
    Questao, !, fail.
nao( _ ).

membro(X, [X|_]).
membro(X, [_|Xs]):-
    membro(X, Xs).

comprimento(S,N) :- length(S,N).

estaVazia(L,V) :- comprimento(L,V),nao(V>0).

imprime([]).
imprime([X|T]) :- write(X), nl, imprime(T).

getArco(Origem,Destino,Dist) :- arco(Origem,Destino,Dist).

adjacente(X,Y,Z) :- getArco(X,Y,Z).
adjacente(X,Y,Z) :- getArco(Y,X,Z).

%reverte uma lista
reverseL(Ds,Es) :- reverseL(Ds, [], Es).
reverseL([],Ds,Ds).
reverseL([D|Ds],Es,Fs) :- reverseL(Ds, [D|Es], Fs).

getLixo(Local, Tipo) :- pontos_recolha(_,_,Local,Tipo,_,_).

seleciona(E, [E|Xs],Xs).
seleciona(E, [X|Xs], [X|Ys]) :- seleciona(E, Xs, Ys).

% --------------------------------------------------- Trabalho Individual -----------------------------------------------------------

inicial('R do Alecrim').

final('Av 24 de Julho').

% --------------------------------------------------- *Procura Não Informada* ----------------------------------------------------------- 

% --------------------------------------------------- DFS - Profundidade ----------------------------------------------------------- 

solveDepthFirst([Nodo|Caminho]):-
    inicial(Nodo),
    depthFirst(Nodo,[Nodo],Caminho).

depthFirst(Nodo,_, []):-
    final(Nodo).

depthFirst(Nodo, Historico, [NodoProx|Caminho]):-
    adjacente(Nodo, NodoProx, _),
    nao(membro(NodoProx, Historico)),
    depthFirst(NodoProx, [NodoProx|Historico], Caminho).

solveDepthFirstAll(L):- findall((S,C), (solveDepthFirst(S), comprimento(S,C)),L).


% --------------------------------------------------- DFS - Profundidade com Custos ----------------------------------------------------------- 

resolveDPdistancia(L):- findall((S,C), (solveDepthFirstC(S,C)), L).

solveDepthFirstC([Nodo|Caminho], Custo):-
    inicial(Nodo),
    depthFirstC(Nodo,[Nodo], Caminho, Custo).

depthFirstC(Nodo, _, [], 0):-
    final(Nodo).

depthFirstC(Nodo, Historico, [NodoProx|Caminho], Custo):-
    adjacente(Nodo, NodoProx, CustoArco),
    nao(membro(NodoProx, Historico)),
    depthFirstC(NodoProx, [NodoProx|Historico], Caminho, CustoAux),
    Custo is CustoArco + CustoAux.


% ------------------------------------------------ DFS Tipos -------------------------------------------------------------------

solveDepthFirstT([Nodo|Caminho], Tipo):-
    inicial(Nodo),
    depthFirstT(Nodo,[Nodo],Caminho, Tipo).

depthFirstT(Nodo,_,[], Tipo):-
    final(Nodo).

depthFirstT(Nodo, Historico, [NodoProx|Caminho], Tipo):-
    getLixo(NodoProx, Tipo),
    adjacente(Nodo, NodoProx,_),
    nao(membro(NodoProx, Historico)),
    depthFirstT(NodoProx, [NodoProx|Historico], Caminho, Tipo).

solveDepthFirstTAll(L,T):- findall((S,C), (solveDepthFirstT(S, T), comprimento(S,C)),L).

% --------------------------------------------------- BFS - Largura ----------------------------------------------------------- 


resolveBFS(No,Solucao) :-
    inicial(No),
    breadthFirst([[No]],Sol),
    reverseL(Sol,Solucao).


breadthFirst([[No|Caminho]|_],[No|Caminho]) :-
    final(No).

breadthFirst([Caminho|Caminhos],Solucao) :-
    estender(Caminho,NovosCaminhos),
    append(Caminhos,NovosCaminhos,Caminhos1),
    breadthFirst(Caminhos1,Solucao).

estender([No|Caminho],NovosCaminhos) :-
    findall([NovoNo,No|Caminho],
        (adjacente(No,NovoNo,_), \+ membro(NovoNo,[No|Caminho])),
        NovosCaminhos).

%resolveBFSAll(L) :- findall((S,C), (resolveBFS(S,C)), L).


% --------------------------------------------------- BFS - Largura Custos ----------------------------------------------------------- 

%NÃO ESTÁ A FUNCIONAR

resolveBFScusto(No,Solucao, Custo) :-
    inicial(No),
    breadthFirstC([[No]],Sol, Custo),
    reverseL(Sol,Solucao).

breadthFirstC([[No|Caminho]|_],[No|Caminho], 0) :-
    final(No).

breadthFirstC([Caminho|Caminhos],Solucao, Custo) :-
    extender(Caminho,NovosCaminhos,CustoArco),
    append(Caminhos,NovosCaminhos,Caminhos1),
    breadthFirstC(Caminhos1,Solucao, CustoAux),
    Custo is CustoArco + CustoAux.

extender([No|Caminho],NovosCaminhos, CustoArco) :-
    findall([NovoNo,No|Caminho],
    (adjacente(No,NovoNo,CustoArco), + membro(NovoNo,[No|Caminho])),
    NovosCaminhos).


% ----------------------------------------- Busca Iterativa Limitada em Profundidade ----------------------------------------------------------- 


resolveDPlimitada(Solucao,L) :-
    inicial(No),
    depthFirstLimited([],No,Sol,L),
    reverseL(Sol,Solucao).

depthFirstLimited(Caminho,No,[No|Caminho],_) :-
    final(No),!.

depthFirstLimited(Caminho,No,S,L) :-
    L > 0,
    adjacente(No,No1,_),
    \+ membro(No1,Caminho),
    L1 is L - 1,
    depthFirstLimited([No|Caminho],No1,S,L1).

resolveDPlimitadaAll(R,Size):- findall((S,C), (resolveDPlimitada(S,Size), comprimento(S,C)),R).

% -------------------------------------- Busca Iterativa Limitada em Profundidade com Custos ----------------------------------------------------------- 

resolveDPlimitadaC(Solucao,L,Custo) :-
    inicial(No),
    depthFirstLimitedC([],No,Sol,L,Custo),
    reverseL(Sol,Solucao).

depthFirstLimitedC(Caminho,No,[No|Caminho],L,0) :-
    final(No),!.

depthFirstLimitedC(Caminho,No,S,L, Custo) :-
    L > 0,
    adjacente(No,No1,CustoArco),
    \+ membro(No1,Caminho),
    L1 is L - 1,
    depthFirstLimitedC([No|Caminho],No1,S,L1, CustoAux),
    Custo is CustoArco + CustoAux.

solveDepthFirstLimitedCostAll(L,Size,Custo) :- 
    findall((S,C,Cost), (resolveDPlimitadaC(S,Size,Cost), comprimento(S,C)), L).


% -------------------------------------- Busca Iterativa Limitada em Profundidade com Tipos ----------------------------------------------------------- 

solveDPlimitadaT(Solucao,L, Tipo) :-
    inicial(No),
    depthFirstLimitedT([],No,Sol,L,Tipo),
    reverseL(Sol,Solucao).

depthFirstLimitedT(Caminho,No,[No|Caminho],L, Tipo) :-
    final(No),!.

depthFirstLimitedT(Caminho,No,S,L,Tipo) :-
    L > 0,
    getLixo(No1, Tipo),
    adjacente(No,No1,_),
    \+ membro(No1,Caminho),
    L1 is L - 1,
    depthFirstLimitedT([No|Caminho],No1,S,L1,Tipo).

% ----------------------------------------------------- Gulosa ----------------------------------------------------------- 

solveGulosa(Nodo, Caminho/Custo) :-
    estimativa(Nodo, Estimativa),
    agulosa([[Nodo]/0/Estimativa], CaminhoInverso/Custo/_),
    reverseL(CaminhoInverso, Caminho).

agulosa(Caminhos, Caminho) :-
    obtem_melhor_g(Caminhos, Caminho),
    Caminho = [Nodo|_]/_/_,final(Nodo).

agulosa(Caminhos, SolucaoCaminho) :-
    obtem_melhor_g(Caminhos, MelhorCaminho),
    seleciona(MelhorCaminho, Caminhos, OutrosCaminhos),
    expandeGulosa(MelhorCaminho, ExpCaminhos),
    append(OutrosCaminhos, ExpCaminhos, NovoCaminhos),
    agulosa(NovoCaminhos, SolucaoCaminho).

obtem_melhor_g([Caminho], Caminho) :- !.

obtem_melhor_g([Caminho1/Custo1/Est1,_/Custo2/Est2|Caminhos], MelhorCaminho) :-
    Est1 =< Est2, !,
    obtem_melhor_g([Caminho1/Custo1/Est1|Caminhos], MelhorCaminho).

obtem_melhor_g([_|Caminhos], MelhorCaminho) :- 
    obtem_melhor_g(Caminhos, MelhorCaminho).

expandeGulosa(Caminho, ExpCaminhos) :-
    findall(NovoCaminho, adjacenteG(Caminho,NovoCaminho), ExpCaminhos).

adjacenteG([Nodo|Caminho]/Custo/_, [ProxNodo,Nodo|Caminho]/NovoCusto/Est) :-
    adjacente(Nodo, ProxNodo, PassoCusto),
    nao(membro(ProxNodo, Caminho)),
    NovoCusto is Custo + PassoCusto,
    estimativa(ProxNodo, Est).

estimativa(Nodo,Est):-
    distance(Nodo, Est).

distance(Origem,Dis):-
    pontos_recolha(Lat1,Lon1,Origem,_,_,_),
    final(Destino),
    pontos_recolha(Lat2,Lon2,Destino,_,_,_),
    P is 0.017453292519943295,
    A is (0.5 - cos((Lat2 - Lat1) * P) / 2 + cos(Lat1 * P) * cos(Lat2 * P) * (1 - cos((Lon2 - Lon1) * P)) / 2),
    Dis is (12742 * asin(sqrt(A))).


% ---------------------------------------------- Gulosa com tipos----------------------------------------------------------- 

resolveGulosaT(Nodo,Caminho/Custo,Tipo):-
    estimativa(Nodo,Estimativa),
    agulosaT([[Nodo]/0/Estimativa],CaminhoInverso/Custo/_,Tipo),
    reverse(CaminhoInverso,Caminho).

agulosaT(Caminhos,Caminho,Tipo):-
    obtem_melhor_g(Caminhos,Caminho),
    Caminho = [Nodo|_]/_/_,final(Nodo).

agulosaT(Caminhos,SolucaoCaminho,Tipo):-
    obtem_melhor_g(Caminhos, MelhorCaminho),
    seleciona(MelhorCaminho, Caminhos, OutrosCaminhos),
    expandeGulosaT(MelhorCaminho,ExpCaminhos,Tipo),
    append(OutrosCaminhos,ExpCaminhos,NovoCaminhos),
    agulosaT(NovoCaminhos,SolucaoCaminho,Tipo).

expandeGulosaT(Caminho,ExpCaminhos,Tipo):-
    findall(NovoCaminho,adjacenteGT(Caminho,NovoCaminho,Tipo),ExpCaminhos).

adjacenteGT([Nodo|Caminho]/Custo/_,[ProxNodo,Nodo|Caminho]/NovoCusto/Est,Tipo):-
    getLixo(NodoProx,Tipo),
    adjacente(Nodo,ProxNodo,PassoCusto),
    nao(membro(ProxNodo,Caminho)),
    NovoCusto is Custo+PassoCusto,
    estimativa(ProxNodo,Est).


% ----------------------------------------------------- A* ----------------------------------------------------------- 

resolveAEstrela(Nodo, Caminho/Custo) :-
    estimativa(Nodo, Estimativa),
    aestrela([[Nodo]/0/Estima], CaminhoInverso/Custo/_),
    reverseL(CaminhoInverso, Caminho),!.

aestrela(Caminhos, Caminho) :-
    obtem_melhor(Caminhos, Caminho),
    Caminho = [Nodo|_]/_/_,final(Nodo).

aestrela(Caminhos, SolucaoCaminho) :-
    obtem_melhor(Caminhos, MelhorCaminho),
    seleciona(MelhorCaminho, Caminhos, OutrosCaminhos),
    expandeAEstrela(MelhorCaminho, ExpCaminhos),
    append(OutrosCaminhos, ExpCaminhos, NovoCaminhos),
    aestrela(NovoCaminhos, SolucaoCaminho).

obtem_melhor([Caminho], Caminho) :- !.

obtem_melhor([Caminho1/Custo1/Est1,_/Custo2/Est2|Caminhos], MelhorCaminho) :-
    Custo1 + Est1 =< Custo2 + Est2, !,
    obtem_melhor([Caminho1/Custo1/Est1|Caminhos], MelhorCaminho).

obtem_melhor([_|Caminhos], MelhorCaminho) :- 
    obtem_melhor(Caminhos, MelhorCaminho).

expandeAEstrela(Caminho, ExpCaminhos) :-
findall(NovoCaminho, adjacenteG(Caminho,NovoCaminho), ExpCaminhos). %adjacenteG/2 da Alínea anterior

