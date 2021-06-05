import pandas as pd
import numpy as np
from math import sin, cos, sqrt, atan2
import re

# Informação sobre o dataset e as suas colunas
# Latitude | Longitude | Ponto Recolha | Contentor Resíduo | Contentor Litros/Quantidade

def calculaDistancia(lat1, long1, lat2, long2):
    R = 6373.0

    dlon = long2 - long1
    dlat = lat2 - lat1

    a = (sin(dlat/2))**2 + cos(lat1) * cos(lat2) * (sin(dlon/2))**2
    c = 2 * atan2(sqrt(a), sqrt(1-a))
    
    return R*c

dataset = pd.read_excel('C:\\Users\\titoa\\OneDrive\\Documents\\Universidade 20-21\\2Semestre\\SRCR\\Projeto_Individual\\dataset.xlsx')

dadosDatasetUniversal = {}

dadosDataset = {}
dadosDataset['Lixos'] = {}
dadosDataset['Vidro'] = {}
dadosDataset['Papel e Cartão'] = {}
dadosDataset['Embalagens'] = {}
dadosDataset['Organicos'] = {}


filePontos = open("pontos_recolha.pl", "w", encoding="UTF-8")
filePontos.write("%%pontos_recolha(latitude, longitude, local, tipo, quantidade, [lista]).\n\n")

fileArcos = open("arcos.pl", "w", encoding="UTF-8")
fileArcos.write("%%arco(rua 1, rua 2, distancia).\n\n")


for line in dataset.values:
    
    pontoDeRecolha = line[2]
    res = re.search(r'([\w, -\/]+)(\[(.*)\])?', pontoDeRecolha)
    rua = re.split(r',', res[1])[0]
    if rua[-1] == " ":
        rua = rua[:-1]
    
    if rua not in dadosDataset[line[3]]:
         dadosDataset[line[3]][rua] = {
            'latitude':line[0],
            'longitude':line[1],
            'destinos':[],
            'quantidade':0
            }
    dadosDataset[line[3]][rua]['quantidade'] += line[4] 
    
    if rua not in dadosDatasetUniversal:
        dadosDatasetUniversal[rua] = {
                'latitude':line[0],
                'longitude':line[1],
                'destinos':[],
                'quantidade':0
        }
    dadosDatasetUniversal[rua]['quantidade'] += line[4] 
    
    caminhos = ""
    if res[3] is not None:
        destinos = re.split(r' ?, ?', res[3])
        for destino in destinos:
            caminhos += "'" + destino + "',"
            dadosDataset[line[3]][rua]['destinos'].append(destino)
            if destino not in dadosDatasetUniversal[rua]['destinos']:
                dadosDatasetUniversal[rua]['destinos'].append(destino)
        caminhos = caminhos[:-1]
    

    filePontos.write("pontos_recolha("  + str(line[0])
                                + ", "  + str(line[1])
                                + ", '" + rua + "'"
                                + ", '" + line[3] + "'"
                                + ", "  + str(line[4])
                                + ", [" + str(caminhos)
                                + "]).\n")

dadosDatasetUniversal['Bqr dos Ferreiros'] = {'latitude':-9.149144, 'longitude':38.708209, 'destinos':[], 'quantidade':0}
dadosDatasetUniversal['Av Dom Carlos I']   = {'latitude':-9.153078, 'longitude':38.709049, 'destinos':[], 'quantidade':0}

percorridos = []

for (rua,valor) in dadosDatasetUniversal.items():
    for destino in valor['destinos']:
        if destino in dadosDatasetUniversal:
            par = (rua,destino)
            par2 = (destino,rua)
            if par not in percorridos and par2 not in percorridos and rua != destino:
                percorridos.append(par)
                fileArcos.write("arco(" + "'" + rua + "'" +
                            ", " + "'" + destino + "'" +
                            ", " + str(calculaDistancia(valor['latitude'],valor['longitude'],
                                            dadosDatasetUniversal[destino]['latitude'],
                                            dadosDatasetUniversal[destino]['longitude'])) +
                            ").\n")

                            
#fileArcos.write("%%arco(rua 1, rua 2, distancia, tipo).\n\n")




#print(dadosDataset)
# Ponto de Recolha : X
# Dicionario :
#   - Nome da Rua: resto da informação
# Arcos :
#   - Rua A para Rua B : coordenadas A e coordenadas B