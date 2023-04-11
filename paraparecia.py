def load_table ():
    reconocer = open ("./terminos.txt") 
    leer = reconocer.readlines () 
    tipos = {}
    for i in leer:
        line = i.strip().split("\t")
        experimento = f"{line [0]}.csv" 
        muestra = f"{line [1]}_sample_table.txt" 
        tipo = line [2]
        if experimento not in tipos: 
            tipos[experimento] = {"CASO": [], "CONTROL" : []} 
            if tipo == "CASO":
                tipos [experimento]["CASO"]+= [muestra] 
            else:
                tipos [experimento]["CONTROL"]+= [muestra]
        else: 
            if tipo == "CASO":
                tipos [experimento]["CASO"]+= [muestra]   
            else:
                tipos [experimento]["CONTROL"]+= [muestra]
    return (tipos)

def symbols ():
    file = open ("./symbols.txt")
    rf= file.readlines()
    sim =[]
    for i in rf:
        s = i.replace (",","\t")
        s = s.replace (" ","")
        s = s.strip().split ("\t")
        gen = s[-1]
        sim.append(gen)
    return (sim)

def getnames (tabla,list_genes, tipo):
    gens_dic ={}
    for gen in list_genes:
        gens_dic[gen]= ["not found"]
    print (f"corriendo{tabla}")
    genes = symbols ()
    soft = tabla.split()[1].replace(".csv","_family.soft")
    soft = f"GSE{soft}"
    file = open(f"./GSE/{soft}")
    rf = file.readline ()
    out = open (f"genes{tipo}_{tabla}","w")
    while rf != "":
        id = rf.replace("\t"," ").split()[0]
        if id in list_genes:
            allinfo = rf.replace("\t"," ").split()
            for keygen in genes:
                if keygen in allinfo:
                    if keygen not in gens_dic[id]:
                        gens_dic[id]= [keygen]              
        rf = file.readline ()
    for gen in gens_dic:
        info = ""
        for name in gens_dic[gen]:
            info+= f";{name}"
        out.write(f"{gen}{info}\n")
    out.close ()
        
def reads_csv (tipos): 
   
    tablas = os.listdir ("./Todos") 
    for tabla in tablas: 
        if "GEOD" in tabla: 
            df = pd.read_csv(f"./Todos/{tabla}",sep = "\t") 
            df.set_index ("id",inplace= True)
            Casos = tipos [tabla]["CASO"] 
            Controles = tipos [tabla]["CONTROL"]
            dfControl = df.loc[:,Controles]
            dfCasos = df.loc [:,Casos]

            promedioca = dfCasos.mean(axis=1) 
            dfCasos.loc[:,("Average")] = promedioca
            promedioco = dfControl.mean(axis=1) 
            dfControl.loc[:,("Average")] = promedioco
            dfCasos ["Fold_Change"] = dfCasos ["Average"] /(dfControl ["Average"] - 1)
            dfCasos ["Logfc"] = np.log(dfCasos ["Fold_Change"])
            dfCasos2 = dfCasos [~dfCasos.isnull().any(axis =1)]
            q95 = dfCasos2 ["Fold_Change"]. quantile (.95)
            q05= dfCasos2 ["Fold_Change"]. quantile (.05)
            dfhigh = dfCasos2 [dfCasos2 ["Fold_Change"] >= q95]
            dflow = dfCasos2 [dfCasos2 ["Fold_Change"] <= q05]
            list_high = list(dfhigh.index.values)
            list_low = list(dflow.index.values)

            dfhigh.to_csv (f"./Fold_Change/{tabla}_high.csv")
            dflow.to_csv (f"./Fold_Change/{tabla}_low.csv")

            getnames (tabla, list_high,"high")
            getnames (tabla, list_low,"low")

            print (f"{tabla}, {dfhigh.shape}, {dflow.shape}")
         
def makematrix (dicmatrix, folder):
    df = pd.DataFrame.from_dict(dicmatrix) 
    df.to_csv (f"./Todos/{folder}.csv", sep = "\t", index=False) 

def cargar ():
    folders = os.listdir ("./") 
    for folder in folders: 
        if "GEO" in folder:
            subfolders = os.listdir (f"./{folder}") 
            dicmatrix = {"id": []}
            for subfolder in subfolders:
                if "Proccesed" in subfolder: 
                    files = os.listdir (f"./{folder}/{subfolder}") 
                    address = f"./{folder}/{subfolder}" 
                    for file in files: 
                        f = open (f"{address}/{file}")
                        rf = f.readlines() 
                        if file not in dicmatrix:
                            dicmatrix [file] = []
                        for i in rf:
                            if i != "":
                                line = i.strip().split() 
                                valor = line [1]

                                dicmatrix [file] += [valor]
                                if len(dicmatrix ["id"]) < len (rf):
                                    dicmatrix ["id"] += [line [0]]
            dicmatrix2 = {}
            for i in dicmatrix:
                dicmatrix2 [i] = dicmatrix [i][1:]
            makematrix (dicmatrix2,folder)
                       
import os
import pandas as pd
import matplotlib
import numpy as np 

cargar ()
table = load_table () 
reads_csv (table)