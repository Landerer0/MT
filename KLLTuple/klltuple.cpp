#include <iostream>
#include <cmath>
#include <algorithm>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <map>
#include <chrono>

#include "klltuple.hpp"
#include "minheap.hpp"

void debugLetra(string string){
    cout << string << endl;
}

using namespace std;

void KLLTuple::iniciarHeap(int numNiveles){
    heap = MinHeap(0,numNiveles);
}   

void KLLTuple::resetHeap(int numNiveles){
    heap = MinHeap(0,numNiveles);
}   

void KLLTuple::setH_pp(int newLevel){
    int nivelesARetirar = newLevel-H_pp;
    if(nivelesARetirar < 0){        // old H'' > new H'': no se realizara ningun proceso, pues si se necesitan mas niveles estos seran agregados
        nivelesARetirar = 0;
    }
    //else if(nivelesARetirar > 0) // old H'' < new H'': eliminar niveles pq menos elementos van a ser procesados (ocupara menos espacio)
    for(int i=0;i<nivelesARetirar;i++){
        if(sketch.size()<=1){
            newLevel--;    
        } 
        else {
            sketch.pop_back();
            numArreglos--;
        }
    }

    H_pp =  newLevel;
    wH_pp = pow(2,H_pp);
    mascara = pow(2,H_pp);
    numArreglos = sketch.size();
    H = H_pp+numArreglos-1;
}

void KLLTuple::setupKLL(uint64_t nP, double epsilonP, double deltaP, double cP, int minkP){
    iniciarHeap(7);
    isMrl = false;
    n = nP;
    epsilon = epsilonP;
    delta = deltaP;
    c = cP;
    minK = minkP;
    metodoCompactacion = "ordenada"; // por defecto es ordenada

    minElement = make_pair(std::numeric_limits<double>::max(),-2);
    maxElement = make_pair(-1*minElement.first,-2);
    sampleElement=make_pair(0,0);
    sampleWeight=0;

    numElementosRevisados = 0;
    numTotalElementos = 0;

    H = 1.0*log2(epsilon*n); // O(log(epsilon*n))
    k = 1.0*(double)(1/epsilon) * log2(log2((1/delta)));
    s = 1;

    H_p = H-s;
    //H_pp = H - ceil((double)log(k) / (double)log((1/c)) ); // H_pp = H - log(k)
    H_pp = H - ceil((double)log2((1.0/epsilon)) ); // H_pp = H - log(k)
    if(H_pp<0) H_pp = 0;

    wH_pp = pow(2,H_pp);
    mascara = pow(2,H_pp);
    
    numArreglos = (H - H_pp+1);
}

void KLLTuple::setupMRL(int minkP){
    iniciarHeap(7);
    isMrl = true;
    n = 0;
    minK = minkP;
    k = minK;
    epsilon = -1;
    delta = -1;
    c = 1;
    metodoCompactacion = "ordenada"; // por defecto es ordenada

    minElement = make_pair(std::numeric_limits<double>::max(),-2);
    maxElement = make_pair(-1*minElement.first,-2);
    sampleElement=make_pair(0,0);
    sampleWeight=0;

    numElementosRevisados = 0;
    numTotalElementos = 0;

    numArreglos=1;
    H=0;
    H_pp=0;
    s=1;
    
    wH_pp = pow(2,0);
    mascara = pow(2,0);
}

void KLLTuple::startSketch(){ // necesita haberse inicializado k, mink, c, H y H_pp
    vector<long> sizeTemp;
    for(int i=H_pp;i<=H;i++){
        uint64_t cantElementos;
        //cout << i << " " << H-s << endl;
        if(isMrl){
            cantElementos = k;
        }
        else{
            if(i>(H-s)) cantElementos = max(k,(long)minK);
            else cantElementos = max((int)(k*pow(c,H-i-1)),(int)minK);
        }
        //if(i>(H-s)) cantElementos = max(k,(long)minK);
        //else cantElementos = max((int)(k*pow(c,H-i-1)),(int)minK);
        // cout << i << " El: " << cantElementos << endl;
        //cout << "cantElementos en arreglo " << i+1 << ": " << cantElementos << endl;
        sizeTemp.push_back(cantElementos);
    }

    // inicializar los vectores de tam k*c^lvl
    for(int i=0;i<sizeTemp.size();i++){
        // el valor por defecto es -1, que indica "vacio"
        uint64_t cantElementos = sizeTemp.at(i);
        //cerr << "Cantidad Elementos arreglo " << i << " (nivel " << i+H_pp+1 <<") :" << cantElementos<< endl;
        pair<int64_t,int64_t> valorElemento = make_pair(-2,-2);
        vector<pair<int64_t,int64_t>> vectorAtLvlI(cantElementos,valorElemento); 
        pair<vector<pair<int64_t,int64_t>>,long> toInsert;
        toInsert.first=vectorAtLvlI;
        toInsert.second=0; // representa el num de elementos ocupados en el arreglo
        sketch.push_back(toInsert);
    }
}

void KLLTuple::startLimitedSketchMRL(uint64_t np, uint64_t espacioMaximo){
    uint64_t numNivelesNecesarios = -1;
    uint64_t minKNecesario = 0;
    uint64_t numElementosP = np;
    uint64_t espacioOcupado = espacioMaximo;
    uint64_t cantElementosPorNivel;
    bool parametrosEncontrados = false;

    espacioOcupado-=sizeInBytes();
    uint64_t espacioIteracion = 0;
    uint64_t espacioPorNivel;
    while(!parametrosEncontrados){
        numNivelesNecesarios++;

        uint64_t cantElementosIteracion;
        cantElementosIteracion = espacioOcupado-sizeof(long)*(numNivelesNecesarios+1);
        cantElementosIteracion /= ((numNivelesNecesarios+1) * sizeof(double));

        if(cantElementosIteracion < 2){
            cerr << "ERROR EN LA CREACION DEL MRL" << endl;
            return;
        }

        espacioPorNivel = (cantElementosIteracion*sizeof(double))+sizeof(long);
        uint64_t numElementosAbarcados = 0;

        espacioIteracion=0;
        for(int i = 0; i<= numNivelesNecesarios;i++){
            espacioIteracion+=espacioPorNivel;
        }
        numElementosAbarcados = cantElementosIteracion*pow(2,numNivelesNecesarios)-1;
        cerr << "NumNiveles: " << numNivelesNecesarios+1 << " NumElementosAbarcados: " << numElementosAbarcados << "/" << np<<endl;
        cerr << "\t Espacio ocupado: " << espacioIteracion  << "/" << espacioOcupado << " Cantidad de elementos por nivel: " << cantElementosIteracion << " Espacio por nivel: " << espacioPorNivel << " ... " << espacioPorNivel*(numNivelesNecesarios+1) << endl;

        cantElementosPorNivel = cantElementosIteracion;
        if(numElementosAbarcados>=np){
            parametrosEncontrados = true;
            n = numElementosAbarcados;
        }
    }

    cerr << "Num niveles: " << numNivelesNecesarios+1 << endl;
    cerr << "cant de elementos por nivel: " << cantElementosPorNivel << endl;

    minK = cantElementosPorNivel;
    k = minK;

    vector<long> sizeTemp;
    unsigned long long cantElementos;
    cantElementos = minK;
    for(int i=0;i<=numNivelesNecesarios;i++) sizeTemp.push_back(cantElementos);
    
    // inicializar los vectores de tam k*c^lvl
    for(int i=0;i<sizeTemp.size();i++){
        // el valor por defecto es -1, que indica "vacio"
        unsigned long long cantElementos = sizeTemp.at(i);
        espacioOcupado += cantElementos;
        //cerr << "Cantidad Elementos arreglo " << i << " (nivel " << i+H_pp+1 <<") :" << cantElementos<< endl;
        pair<int64_t,int64_t> valorElemento = make_pair(-2,-2);
        vector<pair<int64_t,int64_t>> vectorAtLvlI(cantElementos,valorElemento); 
        pair<vector<pair<int64_t,int64_t>>,long> toInsert;
        toInsert.first=vectorAtLvlI;
        toInsert.second=0; // representa el num de elementos ocupados en el arreglo
        sketch.push_back(toInsert);
    }

    H = numNivelesNecesarios; 
    numArreglos = H+1;
}

void KLLTuple::startLimitedSketchKLL(uint64_t np, uint64_t espacioMaximo){
    vector<long> sizeTemp;
    for(int i=H_pp;i<=H;i++){
        unsigned long long cantElementos;
        if(i>(H-s)) cantElementos = max(k,(long)minK);
        else cantElementos = max((int)(k*pow(c,H-i-1)),(int)minK);
        cout << i << " El: " << cantElementos << endl;
        //cout << "cantElementos en arreglo " << i+1 << ": " << cantElementos << endl;
        sizeTemp.push_back(cantElementos);
    }

    // HACER CALCULO AQUI DEL NUMERO DE ELEMENTOS ABARCADOS, en caso que no abarque aumentar H y H_pp en 1 
    uint64_t numElementosAbarcados = 0;
    uint64_t numVecesLlenaCompactor = 1;
    for(int i=sizeTemp.size()-2;i>=0;i--){
        numVecesLlenaCompactor = ceil( ((double)sizeTemp.at(i+1)*(double)numVecesLlenaCompactor) / ceil((double)sizeTemp.at(i)/2.0) );
    }
    cerr << "numVeces que se llena nivel inferior: " << numVecesLlenaCompactor << endl;
    uint64_t numElementosTotaleskll = pow(2,H_pp)*numVecesLlenaCompactor;
    while(numElementosTotaleskll < np){
        cerr << "H: " << H << " H_pp: " << H_pp << " numElementosAbarcados: " << numElementosTotaleskll << "/" << np<< endl;
        H_pp++;
        H++;
        wH_pp = pow(2,H_pp);
        mascara = pow(2,H_pp);
        numElementosTotaleskll *= 2;
    }
    //cerr << "H: " << H << " H_pp: " << H_pp << " numElementosAbarcados: " << numElementosTotaleskll << "/" << numElementsParam<< endl;
    
    wH_pp = pow(2,H_pp);
    mascara = pow(2,H_pp);

    // inicializar los vectores de tam k*c^lvl
    for(int i=0;i<sizeTemp.size();i++){
        unsigned long long cantElementos = sizeTemp.at(i);
        pair<int64_t,int64_t> valorElemento = make_pair(-2,-2);
        vector<pair<int64_t,int64_t>> vectorAtLvlI(cantElementos,valorElemento); 
        pair<vector<pair<int64_t,int64_t>>,long> toInsert;
        toInsert.first=vectorAtLvlI;
        toInsert.second=0; // representa el num de elementos ocupados en el arreglo
        sketch.push_back(toInsert);
    }

    cout << endl;
}
        
// KLLTuple sin espacio limitado
KLLTuple::KLLTuple(uint64_t numElementsParam, double epsilonParam, double deltaParam, double cParam, int minKp){
    //cout << "KLLTuple Basico" << endl;
    espacioCte = false;
    espacioLimitado = false;

    setupKLL(numElementsParam, epsilonParam, deltaParam, cParam, minKp);
    startSketch();
}

KLLTuple::KLLTuple(vector<int64_t> sizeOfCompactors, int H_ppParam){
    //cout << "KLLTuple Determinado" << endl;
    isMrl = false;
    espacioCte = false;
    espacioLimitado = false;
    H_pp = H_ppParam;

    vector<int64_t> sizeTemp = sizeOfCompactors;
    for(int i=0;i<sizeTemp.size();i++){
        uint64_t cantElementos = sizeTemp.at(i);
        pair<int64_t,int64_t> valorElemento = make_pair(-2,-2);
        vector<pair<int64_t,int64_t>> vectorAtLvlI(cantElementos,valorElemento); 
        pair<vector<pair<int64_t,int64_t>>,long> toInsert;
        toInsert.first=vectorAtLvlI;
        toInsert.second=0; // representa el num de elementos ocupados en el arreglo
        sketch.push_back(toInsert);
    }
    k = sizeTemp.at(sizeTemp.size()-1);
    wH_pp = pow(2,H_pp);
    mascara = pow(2,H_pp);
    numArreglos = sketch.size();
    H = H_pp+numArreglos-1;
    
}


KLLTuple::KLLTuple(uint64_t numElementsParam, double epsilonParam, double cParam, int minKp){ 
    // KLLTuple Tradicional con espacio determinado por epsilon, delta=0.01
    //cout << "KLLTuple proporcionando epsilon y con delta 0.01" << endl;
    espacioCte = false;
    espacioLimitado = true;

    setupKLL(numElementsParam, epsilonParam, 0.01, cParam, minKp);

    //! posteriormente cambiar comportamiento con startSketch y 
    startLimitedSketchKLL(numElementsParam, 15000);
}

KLLTuple::KLLTuple(uint64_t numElementsParam, double epsilonParam, double deltaParam, double cParam, int minKp, bool samplerEnMenorNivel){
    espacioCte = true;
    espacioLimitado = true;
    
    setupKLL(numElementsParam, epsilonParam, deltaParam, cParam, minKp);
    if(samplerEnMenorNivel){
        H -= H_pp;
        H_pp = 0;
    }
    wH_pp = pow(2,H_pp);
    mascara = pow(2,H_pp);

    startSketch();
} // KLLTuple con espacio constante, variable "samplerEnMenorNivel" indica que H'' debe ser 0
        
KLLTuple::KLLTuple(uint64_t minKP){ // MRL sin espacioLimitado
    //cout << "MRLTuple Basico" << endl;
    espacioLimitado = false;
    espacioCte = false;
    setupMRL(minKP);

    startSketch();
}

KLLTuple::KLLTuple(uint64_t espacioMaximo, uint64_t numElementsParam){ // MRL que se le entrega el espacio
    espacioLimitado = true;
    espacioCte = false;
    cout << "Espacio Maximo: " << espacioMaximo << " NumElementosAAbarcar: " << numElementsParam << endl;
    
    setupMRL(2);
    startLimitedSketchMRL(numElementsParam, espacioMaximo);
}

KLLTuple::KLLTuple(double epsilonParam ,uint64_t numElementsParam){
    espacioCte = false;

    int minkP = (1.0/epsilonParam)*ceil(log2(epsilonParam*numElementsParam))+1;

    //double cotaEspacio = 1.0/epsilonParam * log2(1.0/epsilonParam) * log2(1.0/epsilonParam);
    //uint64_t espacioMaximo = ceil(cotaEspacio);
    //cout << "Epsilon Entregado: " << epsilonParam << " Espacio Maximo: " << espacioMaximo << " NumElementosAAbarcar: " << numElementsParam << endl;
    //desde aqui es el mismo constructor que  KLLTuple(uint64_t espacioMaximo, uint64_t numElementsParam)
    espacioLimitado = false;
    isMrl = true;
    epsilon = epsilonParam;
    
    setupMRL(minkP);
    startSketch();

    //startLimitedSketchMRL(numElementsParam,espacioMaximo);

    sizeInBytes();
    print();
}

// creacion KLLTuple con parametros indicados
KLLTuple::KLLTuple(uint64_t HParam,uint64_t sParam,uint64_t H_ppParam,vector<int> sizeNivelesParam, bool isAnMrlParam, bool hasLimitedSpaceParam){
    isMrl = isAnMrlParam;
    espacioLimitado = hasLimitedSpaceParam;
    espacioCte = false;
    H = HParam;
    s = sParam;
    H_p = H-sParam;
    H_pp = H_ppParam;
    wH_pp = pow(2,H_pp);

    mascara = pow(2,H_pp);
    sampleElement=make_pair(0,0);
    sampleWeight=0;
    minElement = make_pair(std::numeric_limits<double>::max(),-2);
    maxElement = make_pair(-1*minElement.first,-2);

    numArreglos = (H - H_pp+1);
    numElementosRevisados = 0;
    numTotalElementos = 0;

    uint64_t espacioOcupado = 0;

    minK = sizeNivelesParam.at(0);
    k = sizeNivelesParam.at(sizeNivelesParam.size()-1);

    // inicializar los vectores de tam k*c^lvl
    for(int i=0;i<sizeNivelesParam.size();i++){
        unsigned long long cantElementos = sizeNivelesParam.at(i);
        espacioOcupado += cantElementos;
        pair<int64_t,int64_t> valorElemento = make_pair(-2,-2);
        vector<pair<int64_t,int64_t>> vectorAtLvlI(cantElementos,valorElemento); 
        pair<vector<pair<int64_t,int64_t>>,long> toInsert;
        toInsert.first=vectorAtLvlI;
        toInsert.second=0; // representa el num de elementos ocupados en el arreglo
        sketch.push_back(toInsert);
    }
} // MRL/KLLTuple que se le indican los parametros que va a ocupar

uint64_t KLLTuple::saveData(string outputFileName){
    //FORMATO:
    // KLLTuple TRADICIONAL: isMRL, minElement, maxElement, metodoCompactacion, n, epsilon, delta, c, minK, numTotalElementos, 
    //                       , sampleWeight, sampleElement
    //                       , numNiveles, numElementosRevisados, vector<posInicialNivel>, vector<elementos> 
    // MRL: isMrl, ,minElement, maxElement, metodoCompactacion minK, numNiveles, numElementosRevisados, vector<posInicialNivel>,
    //                       , vector<elementos> 


    // Abrir el archivo en modo binario
    std::ofstream archivo(outputFileName+".bin", std::ios::binary);
    uint64_t numBytes = 0;

    if (archivo) {
        // Almacenamos el tipo de estructura asociado
        archivo.write(reinterpret_cast<const char *>(&isMrl), sizeof(isMrl));
        numBytes += sizeof(isMrl);
        archivo.write(reinterpret_cast<const char *>(&espacioLimitado), sizeof(espacioLimitado));
        numBytes += sizeof(espacioLimitado);
        archivo.write(reinterpret_cast<const char *>(&minElement), sizeof(minElement));
        numBytes += sizeof(minElement);
        archivo.write(reinterpret_cast<const char *>(&maxElement), sizeof(maxElement));
        numBytes += sizeof(maxElement);
        archivo.write(reinterpret_cast<const char *>(&metodoCompactacion), sizeof(metodoCompactacion));
        numBytes += sizeof(metodoCompactacion);
        if(isMrl){
            archivo.write(reinterpret_cast<const char *>(&minK), sizeof(minK));
    
            numBytes += sizeof(minK);
        } else {
            archivo.write(reinterpret_cast<const char *>(&n), sizeof(n));
            archivo.write(reinterpret_cast<const char *>(&epsilon), sizeof(epsilon));
            archivo.write(reinterpret_cast<const char *>(&delta), sizeof(delta));
            archivo.write(reinterpret_cast<const char *>(&c), sizeof(c));
            archivo.write(reinterpret_cast<const char *>(&minK), sizeof(minK));
            archivo.write(reinterpret_cast<const char *>(&numTotalElementos), sizeof(numTotalElementos));
            archivo.write(reinterpret_cast<const char *>(&sampleWeight), sizeof(sampleWeight));
            archivo.write(reinterpret_cast<const char *>(&sampleElement), sizeof(sampleElement));
            archivo.write(reinterpret_cast<const char *>(&espacioCte), sizeof(espacioCte));
            
            numBytes += sizeof(n);
            numBytes += sizeof(epsilon);
            numBytes += sizeof(delta);
            numBytes += sizeof(c);
            numBytes += sizeof(minK);
            numBytes += sizeof(numTotalElementos);
            numBytes += sizeof(sampleWeight);
            numBytes += sizeof(sampleElement);
            numBytes += sizeof(espacioCte);
        }

        // Almacenar datos necesarios numero de niveles y H_pp
        // H_pp: indica el peso asociado al primer nivel del sketch, necesario para consultas de RANK
        uint32_t numNiveles = sketch.size();
        uint64_t elementosRevisados = numElementosRevisados;

        archivo.write(reinterpret_cast<const char *>(&numNiveles), sizeof(numNiveles));
        archivo.write(reinterpret_cast<const char *>(&elementosRevisados), sizeof(elementosRevisados));
        numBytes += sizeof(numNiveles);
        numBytes += sizeof(elementosRevisados);

        // Almacenar el indice donde inicia el elemento siguiente del sketch respectivo
        uint32_t posInicialNivelActual = 0; 
        for(int i=0;i<sketch.size();i++){
            // cerr << "posInicial: " << posInicialNivelActual << endl;
            archivo.write(reinterpret_cast<const char *>(&posInicialNivelActual), sizeof(posInicialNivelActual));
            posInicialNivelActual+=sketch.at(i).second;
            numBytes += sizeof(posInicialNivelActual);
        }
        archivo.write(reinterpret_cast<const char *>(&posInicialNivelActual), sizeof(posInicialNivelActual));
        numBytes += sizeof(posInicialNivelActual);

        // Almacenar los elementos existentes del sketch
        // Se almacena desde el elemento 0
        for(int i=0;i<sketch.size();i++){
            vector<pair<int64_t,int64_t>> nivelActual = sketch.at(i).first;
            for(int j=0;j<sketch.at(i).second;j++){
                pair<int64_t,int64_t> elementoAGuardar = nivelActual.at(j);
                archivo.write(reinterpret_cast<const char *>(&elementoAGuardar), sizeof(elementoAGuardar));
            }
        }

        // Cerrar el archivo
        archivo.close();

        //std::cout << "Vector guardado en el archivo: " << outputFileName << std::endl;
    } else {
        std::cerr << "No se pudo abrir el archivo: " << outputFileName << std::endl;
        return 0;
    }

    //! PARA VER EL TAMAÑO DEL ARCHIVO
    streampos begin,end;
    ifstream myfile (outputFileName, ios::binary);
    begin = myfile.tellg();
    myfile.seekg (0, ios::end);
    end = myfile.tellg();
    myfile.close();
    // cout << "size primera parte: " << numBytes << " bytes. \n";
    // cout << "size is: " << (end-begin) << " bytes.\n";
    
    //print();
    numBytes= (end-begin);

    return numBytes;
}

// Creacion de MRL a partir de los datos proporcionados en readData
KLLTuple::KLLTuple(uint32_t minKRead, uint32_t numElementosRevisadosRead, vector<vector<pair<int64_t,int64_t>>> niveles, pair<int64_t,int64_t> minElementRead, pair<int64_t,int64_t> maxElementRead, bool espacioLimitadoRead){
    uint32_t numNiveles = niveles.size();
    numArreglos = numNiveles;
    espacioLimitado = espacioLimitadoRead;
    numElementosRevisados = numElementosRevisadosRead;
    numTotalElementos = numElementosRevisadosRead;
    minK = minKRead;
    k = minK;
    isMrl = true;
    H = numNiveles-1;
    H_pp = 0;
    s = numNiveles;
    minElement = minElementRead;
    maxElement = maxElementRead;

    espacioCte = false;

    pair<int64_t,int64_t> valorElementoVacio = make_pair(-2,-2);
    for(int i=0;i<niveles.size();i++){
        pair<vector<pair<int64_t,int64_t>>,long> par;
        par.first = niveles.at(i);
        par.second = niveles.at(i).size();
        for(int j=par.second;j<minK;j++){ // rellenar con elementos vacios
            par.first.push_back(valorElementoVacio);
        }
        // cout << "par second: " << par.second << " de " << par.first.size() << endl;
        sketch.push_back(par);
    }
}

// Creacion de KLLTuple tradicional a partir de los datos proporcionados en readData
KLLTuple::KLLTuple(uint64_t nRead, double epsilonRead,double deltaRead,double cRead, uint32_t minKRead,uint64_t numTotalElementosRead,unsigned long sampleWeightRead,pair<int64_t,int64_t> sampleElementRead, uint64_t numElementosRevisadosRead,vector<vector<pair<int64_t,int64_t>>> niveles, pair<int64_t,int64_t> minElementRead, pair<int64_t,int64_t> maxElementRead, bool espacioLimitadoRead,bool espacioCteRead){
    isMrl = false;
    espacioCte = espacioCteRead;
    espacioLimitado = espacioLimitadoRead;
    minElement = minElementRead;
    maxElement = maxElementRead;


    minK = minKRead;
    n = nRead;
    epsilon = epsilonRead;
    delta = deltaRead;
    c = cRead;

    H = 1.0*log2(epsilon*n); // O(log(epsilon*n))
    k = 1.0*(double)(1/epsilon) * log2(log2((1/delta)));
    s = log2(log2((1/delta)));

    H_p = H-s;
    H_pp = H - ceil((double)log(k) / (double)log((1/c)) ); // H_pp = H - log(k)
    if(H_pp<0) H_pp = 0;
    wH_pp = pow(2,H_pp);

    mascara = pow(2,H_pp);

    mascara = pow(2,H_pp);
    sampleElement=sampleElementRead;
    sampleWeight=sampleWeightRead;

    numArreglos = (H - H_pp+1);
    numElementosRevisados = numElementosRevisadosRead;
    numTotalElementos = numTotalElementosRead;

    vector<long> sizeTemp;
    for(int i=H_pp;i<=H;i++){
        unsigned long long cantElementos;
        if(i>(H-s)) cantElementos = max(k,(long)minK);
        else cantElementos = max((int)(k*pow(c,H-i-1)),(int)minK);
        sizeTemp.push_back(cantElementos);
    }
    if(sizeTemp.size() < niveles.size()){
        for(int i=sizeTemp.size(); i<niveles.size();i++){
            sizeTemp.push_back(max(k,(long)minK));
            H++;
            s++;
        }
    }
        
    pair<int64_t,int64_t> valorElementoVacio = make_pair(-2.0,-2);
    // inicializar los vectores de tam k*c^lvl
    for(int i=0;i<niveles.size();i++){
        pair<vector<pair<int64_t,int64_t>>,long> par;
        par.first = niveles.at(i);
        par.second = niveles.at(i).size();
        for(int j=par.second;j<sizeTemp.at(i);j++){ // rellenar con elementos vacios
            par.first.push_back(valorElementoVacio);
        }
        // cout << "par second: " << par.second << " de " << par.first.size() << endl;
        sketch.push_back(par);
    }
}

KLLTuple KLLTuple::readData(string inputFileName){
    //FORMATO:
    // KLLTuple TRADICIONAL: isMRL, n, epsilon, delta, c, minK, numTotalElementos, sampleWeight, sampleElement
    //                       , numNiveles, numElementosRevisados, vector<posInicialNivel>, vector<elementos> 
    // MRL: isMrl, minK, numNiveles, numElementosRevisados, vector<posInicialNivel>, vector<elementos> 

    std::ifstream archivo(inputFileName, std::ios::binary);
    if (!archivo) {
        std::cout << "No se pudo abrir el archivo." << std::endl;
        return 1;
    } 

    // variables asociadas a la lectura de MRL o KLLTuple tradicional
    // MRL y KLLTuple tradicional:
    bool isMrlRead;
    bool espacioCteRead;
    bool espacioLimitadoRead;
    pair<int64_t,int64_t> minElementRead, maxElementRead;
    string metodoCompactacionRead;
    uint32_t minKRead;
    uint32_t numNiveles;
    uint64_t numElementosRevisadosRead;
    vector<uint32_t> posInicialNivel; // almacena las posiciones asociadas, es auxiliar
    vector<vector<pair<int64_t,int64_t>>> niveles; // almacena los elementos guardados asociados al archivo en binario
    // exclusivos KLLTuple tradicional
    uint64_t nRead; 
    uint64_t numTotalElementosRead; 
    double epsilonRead, deltaRead, cRead;
    unsigned long sampleWeightRead;
    pair<int64_t,int64_t> sampleElementRead;
    
    // Lectura de variables del archivo en binario
    archivo.read(reinterpret_cast<char*>(&isMrlRead), sizeof(isMrlRead));
    archivo.read(reinterpret_cast<char*>(&espacioLimitadoRead), sizeof(espacioLimitadoRead));
    archivo.read(reinterpret_cast<char*>(&minElementRead), sizeof(minElementRead));
    archivo.read(reinterpret_cast<char*>(&maxElementRead), sizeof(maxElementRead));
    archivo.read(reinterpret_cast<char*>(&metodoCompactacionRead), sizeof(metodoCompactacionRead));
    if(isMrlRead){
        archivo.read(reinterpret_cast<char*>(&minKRead), sizeof(minKRead));
    } else {
        archivo.read(reinterpret_cast<char*>(&nRead), sizeof(nRead));
        archivo.read(reinterpret_cast<char*>(&epsilonRead), sizeof(epsilonRead));
        archivo.read(reinterpret_cast<char*>(&deltaRead), sizeof(deltaRead));
        archivo.read(reinterpret_cast<char*>(&cRead), sizeof(cRead));
        archivo.read(reinterpret_cast<char*>(&minKRead), sizeof(minKRead));
        archivo.read(reinterpret_cast<char*>(&numTotalElementosRead), sizeof(numTotalElementosRead));
        archivo.read(reinterpret_cast<char*>(&sampleWeightRead), sizeof(sampleWeightRead));
        archivo.read(reinterpret_cast<char*>(&sampleElementRead), sizeof(sampleElementRead));
        archivo.read(reinterpret_cast<char*>(&espacioCte), sizeof(espacioCte));
    }

    archivo.read(reinterpret_cast<char*>(&numNiveles), sizeof(numNiveles));
    archivo.read(reinterpret_cast<char*>(&numElementosRevisadosRead), sizeof(numElementosRevisadosRead));
    for(int i=0;i<numNiveles+1;i++){
        uint32_t posActual;
        archivo.read(reinterpret_cast<char*>(&posActual), sizeof(posActual));
        posInicialNivel.push_back(posActual);
    }
    
    for(int i=0;i<numNiveles;i++){
        vector<pair<int64_t,int64_t>> nivelActual;
        for(int j=posInicialNivel.at(i);j<posInicialNivel.at(i+1);j++){
            pair<int64_t,int64_t> elementoActual;
            archivo.read(reinterpret_cast<char*>(&elementoActual), sizeof(elementoActual));
            nivelActual.push_back(elementoActual);
        }
        niveles.push_back(nivelActual);
    }

    if(isMrlRead) 
        return KLLTuple(minKRead, numElementosRevisadosRead, niveles, minElementRead, maxElementRead, espacioLimitadoRead);
    else 
        return KLLTuple(nRead,epsilonRead,deltaRead,cRead,minKRead,numTotalElementosRead, sampleWeightRead, sampleElementRead, numElementosRevisadosRead, niveles, minElementRead, maxElementRead, espacioLimitadoRead, espacioCteRead);
}

KLLTuple KLLTuple::copy(){
    KLLTuple copia = isMrl ? KLLTuple(minK, numTotalElementos,this) : KLLTuple(n,epsilon,delta,c,minK,numElementosRevisados,numTotalElementos,this);
    return copia;
}

// COPIA KLLTuple TRADICIONAL
KLLTuple::KLLTuple(uint64_t nCopy, double epsilonCopy, double deltaCopy, double cCopy, uint32_t minKCopy, uint64_t numElementosRevisadosCopy, uint64_t numTotalElementosCopy, KLLTuple* toCopy){
    setupKLL(nCopy, epsilonCopy, deltaCopy, cCopy, minKCopy);
    numElementosRevisados = numElementosRevisadosCopy;
    numTotalElementos = numTotalElementosCopy;
    espacioLimitado = toCopy->hasLimitedSpace();
    espacioCte = toCopy->hasConstantSpace();

    pair<pair<int64_t,int64_t>, pair<int64_t,int64_t>> minMaxElement = toCopy->getMinMaxElement();
    minElement = minMaxElement.first;
    maxElement = minMaxElement.second;

    H_p = H-s;

    H = toCopy->getH();
    H_pp = toCopy->getH_pp();
    numArreglos = H-H_pp+1;
    wH_pp = pow(2,H_pp);
    mascara = pow(2,H_pp);

    pair<unsigned long, pair<int64_t,int64_t>> samplePair = toCopy->getCurrentSample();
    sampleElement=samplePair.second;
    sampleWeight=samplePair.first;

    // copiar los valores del KLLTuple a copiar
    for(int i=0;i<numArreglos;i++){
        sketch.push_back(toCopy->sketchAtLevel(i));
    }
}

// COPIA MRL
KLLTuple::KLLTuple(uint32_t minKCopy, uint64_t numTotalElementosCopy, KLLTuple* toCopy){
    espacioLimitado = toCopy->hasLimitedSpace();
    espacioCte = false;
    setupMRL(minKCopy);

    pair<pair<int64_t,int64_t>, pair<int64_t,int64_t>> minMaxElement = toCopy->getMinMaxElement();
    minElement = minMaxElement.first;
    maxElement = minMaxElement.second;

    pair<unsigned long, pair<int64_t,int64_t>> samplePair = toCopy->getCurrentSample();
    sampleElement=samplePair.second;
    sampleWeight=samplePair.first;

    H = toCopy->getH();
    H_pp = toCopy->getH_pp();
    numArreglos = H+1;
    numElementosRevisados = numTotalElementosCopy;
    numTotalElementos = numTotalElementosCopy;

    for(int i=0;i<numArreglos;i++){
        sketch.push_back(toCopy->sketchAtLevel(i));
    }
}

KLLTuple::~KLLTuple(){
    
}

void KLLTuple::insertElement(long nivel,pair<int64_t,int64_t> &element){
    long posAInsertar = sketch.at(nivel).second;
    sketch.at(nivel).first.at(posAInsertar) = element;
    sketch.at(nivel).second++;
}

void KLLTuple::insertCompactionElement(long nivel,pair<int64_t,int64_t> &element, bool updating){
    long posAInsertar = sketch.at(nivel).second;

    sketch.at(nivel).first.at(posAInsertar) = element;
    sketch.at(nivel).second++;
    
    if(posAInsertar+1==sketch.at(nivel).first.size()) {
        iterativeCompaction(nivel,updating);
    }
}

void KLLTuple::addNewLevel(){
    vector<pair<int64_t,int64_t>> vectorAtLvlI(k,make_pair(-2,-2)); 
    pair<vector<pair<int64_t,int64_t>>,long> toInsert;
    toInsert.first=vectorAtLvlI; // arreglo inicializado
    toInsert.second=0; // pos siguiente elemento a insertar
    sketch.push_back(toInsert);
    numArreglos++;
    H++;
    //cerr << "H-1:"<< H-1 << ", H:" << H << ", H'': " << H_pp << ", numArreglos: " << numArreglos << endl;
}

void KLLTuple::constantSpaceCompaction(){
    print();
    // crear vector< pair<vector<double>,long > >
    vector<pair<vector<pair<int64_t,int64_t>>, long> > auxSketch = sketch;

    // vaciar/limpiar sketch
    for(int nivel=0;nivel<numArreglos;nivel++){
        for(int j=0;j<sketch.at(nivel).second;j++){
            sketch.at(nivel).first.at(j)=make_pair(-2.0,-2);
        }
        sketch.at(nivel).second = 0;
    }

    H++;
    H_pp++;
    wH_pp *= 2;
    mascara *= 2;

    // agregar los elementos de menor nivel al sketch mediante reservoirSampling
    for(int i=0;i<auxSketch.at(0).first.size();i++){
        if(auxSketch.at(0).first.at(i).first>=0){
            if(reservoirKLLSample(auxSketch.at(0).first.at(i),pow(2,H_pp-1)))
                insertElement(0,sampleElement);
        } 
    }

    // agregar el resto de los elementos desde auxSketch a sketch 
    for(int nivel=1;nivel<numArreglos;nivel++){
        for(int i=0;i<auxSketch.at(nivel).second;i++){
            insertCompactionElement(nivel-1,auxSketch.at(nivel).first.at(i),false);
        }
        iterativeCompaction(nivel-1,false);
    }
    print();
    return;
}

vector<pair<int64_t,int64_t>> KLLTuple::seleccionElementosACompactar(vector<pair<int64_t,int64_t>> &elements){
    // cout << metodoCompactacion << endl;

    auto start = std::chrono::high_resolution_clock::now();
    if(metodoCompactacion=="aleatorea") return seleccionAleatoreaOrdenada(elements);
    else return seleccionElementosMayores(elements);
    auto end = std::chrono::high_resolution_clock::now();
    sortTime += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
}


vector<pair<int64_t,int64_t>> KLLTuple::seleccionAleatoreaOrdenada(vector<pair<int64_t,int64_t>> &elements){
    vector<pair<int64_t,int64_t>> toReturn;
    unsigned char elementosPares = 0;
        if(rand()%2==0) elementosPares = 0; // se mantienen los elementos pares
        else elementosPares = 1; // se mantienen los elementos impares
    
    // sort de los elementos
    sort(elements.begin(), elements.end());

    // Agregar los elementos de mayor a menor (en iterativeCompaction)
    for(int i=0;i<elements.size();i++){
        if (i%2==elementosPares) toReturn.push_back(elements.at(i));
        else heap.insert(elements.at(i));
    }

    return toReturn;
}


void KLLTuple::setMetodoCompactacion(string metodoIndicado){
    if(metodoIndicado == "aleatorea") metodoCompactacion = metodoIndicado;
    else if(metodoIndicado == "elementosMayores") metodoCompactacion = metodoIndicado;
    return;
}


string KLLTuple::getMetodoCompactacion(){
    return metodoCompactacion;
}

vector<pair<int64_t,int64_t>> KLLTuple::seleccionElementosMayores(vector<pair<int64_t,int64_t>> &elements){
    vector<pair<int64_t,int64_t>> toReturn;
    
    // sort de los elementos
    sort(elements.begin(), elements.end());
    // Agregar los elementos de mayor a menor (en iterativeCompaction)
    for(int i=0;i<elements.size();i++){
        if (i>=(elements.size()/2)) heap.insert(elements.at(i));
        else toReturn.push_back(elements.at(i));
    }

    return toReturn;
}

bool KLLTuple::reduction(long nivel){
    bool reduction = false;

    pair<vector<pair<int64_t,int64_t>>, long> nivelSketch = sketchAtLevel(nivel);
    vector<pair<int64_t,int64_t>> vectorActual = nivelSketch.first;

    std::sort(vectorActual.begin(), vectorActual.end(), [](const auto& a, const auto& b) {
        return a.second < b.second;
    });

    for(int i=0;i<vectorActual.size()-1;i++){
        if(vectorActual.at(i).second == vectorActual.at(i+1).second){
            vectorActual.at(i+1).first += vectorActual.at(i).first;
            vectorActual.at(i) = make_pair(-1,-1);
            nivelSketch.second--;
            numElementosRevisados-= pow(2,nivel);
            reduction = true;
        }
    }

    if(reduction){
        int n = vectorActual.size();
        int i = 0;  // Índice para recorrer el vector

        // Recorre el vector y mueve los elementos (-1, -1) al final
        for (int j = 0; j < n; j++) {
            if (vectorActual[j].first != -1 || vectorActual[j].second != -1) {
                vectorActual[i] = vectorActual[j];
                i++;
            }
        }

        // Llena el final del vector con elementos (-1, -1)
        while (i < n) {
            vectorActual[i] = { -1, -1 };
            i++;
        }
        sketch.at(nivel) = make_pair(vectorActual,nivelSketch.second);
    }
    return reduction;
}

bool KLLTuple::iterativeCompaction(long nivelInicial, bool updating){ //! queda ver constantSpaceCompaction para ver como finalizar posterior a ello
    long numElementosOcupados = sketch.at(nivelInicial).second;
    long numElementosTotales = sketch.at(nivelInicial).first.size();
    unsigned char elementosPares = 0;

    vector<vector<pair<int64_t,int64_t>>> elementosACompactar;
    bool seDebeCompactar = false;
    int constantSpaceCompacted = 0; // indica si se ha compactado mediante constantSpaceCompaction. 0 indica que no, 1 indica que si 

    if(numElementosOcupados>=numElementosTotales && !reduction(nivelInicial)){
        if(updating) print();
        // Si se lleno y se tiene espacio limitado
        if(nivelInicial==H && espacioLimitado && espacioCte){
            cerr << "inicio constantSpaceCompaction\n";
            constantSpaceCompaction(); // dejar trabajo de compactacion a constantSpaceCompaction
            cerr << "fin constantSpaceCompaction\n";
            return false;
        } 
        else if(nivelInicial==H-H_pp && espacioLimitado && !updating) return true; // indicar que se lleno sketch
        else if(nivelInicial==H-H_pp) addNewLevel(); // se puede agregar un nivel superior
            
        if(updating) cerr << "nivelAgregado:" << nivelInicial << endl;
        elementosACompactar.push_back(seleccionElementosACompactar(sketch.at(nivelInicial).first));
        //cerr << "se compactara nivelInicial: " << nivelInicial << endl;
    } 

    while(!elementosACompactar.empty()){
        long indiceVector = elementosACompactar.size()-1;

        // 1-Reviso si el nivel del stack se encuentra vacio
        if(elementosACompactar.at(indiceVector).empty()){
            if(updating) cerr << "se vacio nivel " << nivelInicial+indiceVector << endl;
            // elimino todos los elementos asociados al nivel que ya fue compactado
            for(int i=0;i<sketch.at(nivelInicial+indiceVector).first.size();i++){
                sketch.at(nivelInicial+indiceVector).first.at(i) = make_pair(-1.0,-1);
            }
            sketch.at(nivelInicial+indiceVector).second=0;
            
            elementosACompactar.pop_back();
            continue;
        }

        // 2-Ingreso en compactor de nivel superior el elemento 
        long nivelAIngresar = nivelInicial + indiceVector + 1;

        //cerr << nivelAIngresar << ", H: " << H << ", H'': " << H_pp  << ", numArreglos: " << numArreglos << "sketchSize: " << sketch.size() << endl;
        long posAInsertar = sketch.at(nivelAIngresar).second; //! AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
        //cerr << "a" << endl;
        
        uint64_t numElementosDisponiblesNivelSuperior = sketch.at(nivelAIngresar).first.size()-sketch.at(nivelAIngresar).second;
        uint64_t numElementosAIngresarNivelSuperior = min(numElementosDisponiblesNivelSuperior,elementosACompactar.at(indiceVector).size());

        // agregamos los elementos a un vector para ser posteriormente añadidos al nivel superior
        vector<pair<int64_t,int64_t>> elementosAIngresar;
        uint64_t offset = elementosACompactar.at(indiceVector).size()-numElementosAIngresarNivelSuperior;
        elementosAIngresar.assign(elementosACompactar.at(indiceVector).begin()+offset, elementosACompactar.at(indiceVector).end());

        // se agregan los elementos al nivel superior
        for(int i=0;i<numElementosAIngresarNivelSuperior;i++){
            sketch.at(nivelAIngresar).first.at(posAInsertar) = elementosAIngresar.at(i);
            posAInsertar++;
        }
        sketch.at(nivelAIngresar).second+=numElementosAIngresarNivelSuperior;

        // eliminamos los elementos asociados al stack que fueron agregados
        for(int i=0;i<numElementosAIngresarNivelSuperior;i++){
            elementosACompactar.at(indiceVector).pop_back();
        }

        // 3-Veo si es necesario compactar el nivel superior
        numElementosOcupados = sketch.at(nivelAIngresar).second;
        numElementosTotales = sketch.at(nivelAIngresar).first.size();
        if(numElementosOcupados>=numElementosTotales){
            // Veo si se lleno y se tiene espacio limitado
            if(nivelAIngresar==H-H_pp && espacioLimitado && espacioCte){
                constantSpaceCompaction(); // dejar trabajo de compactacion a constantSpaceCompaction
                constantSpaceCompacted = 1;
            } 
            else if(nivelAIngresar==H-H_pp && espacioLimitado & !updating) return true; // indicar que se lleno sketch
            else if(nivelAIngresar==H-H_pp) addNewLevel();
            
            elementosACompactar.push_back(seleccionElementosACompactar(sketch.at(nivelAIngresar).first));
        } 
    }
    
    return false;
}

bool KLLTuple::sample(pair<int64_t,int64_t> element){
    if(rand()%wH_pp==0) return true;
    return false;
}

bool KLLTuple::reservoirKLLSample(pair<int64_t,int64_t> element, uint64_t elementWeight){
    uint64_t currentWeight = sampleWeight;
    sampleWeight = elementWeight + sampleWeight;

    if(sampleWeight<=wH_pp){
        double probReemplazo = (double) elementWeight/(double) sampleWeight;
        if((double) (rand()/RAND_MAX) <= probReemplazo) sampleElement = element;
    }
    if(sampleWeight==wH_pp){
        sampleWeight = 0;
        return true;
    }
    else if (sampleWeight>wH_pp){
        pair<int64_t,int64_t> heavyElement;
        // sampler discards the heavier item and keeps the lighter item with weight min{w,v}
        if(elementWeight<currentWeight){
            heavyElement = sampleElement;
            sampleWeight = elementWeight;
            sampleElement = element;
        } else {
            heavyElement = element;
            sampleWeight = currentWeight;
        }
        //with probability max{w,v}/2^h it also outputs the heavier item with weight 2^h
        double probIngresoElementoPesado = (double) max(currentWeight,elementWeight)/(double) sampleWeight;
        if((double) (rand()/RAND_MAX) <= probIngresoElementoPesado) addv(heavyElement);
    }

    return false;
}

bool KLLTuple::add(pair<int64_t,int64_t> element){
    return add(element,1);
}

bool KLLTuple::add(pair<int64_t,int64_t> element, uint32_t elementWeight){
    auto start = std::chrono::high_resolution_clock::now();
    numTotalElementos+=elementWeight;
    if(element<minElement) minElement = element;
    else if (element>maxElement) maxElement = element;
    if(!reservoirKLLSample(element,elementWeight)) return false;
    insertElement(0,sampleElement);
    numElementosRevisados+=elementWeight; // para metodo quantile
    auto end = std::chrono::high_resolution_clock::now();
    buildTime += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    return iterativeCompaction((long) 0, false);
}

// add utilizado cuando se hacer merge de KLLTuple
void KLLTuple::addv(pair<int64_t,int64_t> element){
    insertElement(0,element);
    iterativeCompaction((long) 0, true);

    return;
}

vector<pair<pair<int64_t,int64_t>,uint64_t>> KLLTuple::obtenerResumenElementos(){
    vector<pair<pair<int64_t,int64_t>,uint64_t>> vectorElementos; // par(elemento,peso) 
    // llenar el vector con los Elementos del sketch 
    for(int i=0;i<numArreglos;i++){
        uint64_t peso = pow(2,i);
        for(int j=0;j<sketch.at(i).first.size();j++){
	        if(sketch.at(i).first.at(j).first < 0)continue;
            pair<pair<int64_t,int64_t>, uint64_t> toInsert;
            toInsert.first = sketch.at(i).first.at(j);
            toInsert.second = peso;
            vectorElementos.push_back(toInsert);
        }
    }

    // realizar el sort
    auto start = std::chrono::high_resolution_clock::now();
    sort(vectorElementos.begin(),vectorElementos.end());
    auto end = std::chrono::high_resolution_clock::now();
    sortTime += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

    return vectorElementos;
}

uint64_t KLLTuple::rank(int64_t element){
    auto startS = std::chrono::high_resolution_clock::now();
    uint64_t rank = 0;

    vector<pair<int64_t,int64_t>> actual;

    for(int nivel=0;nivel< numArreglos;nivel++){ // por cada arreglo
        actual = sketch.at(nivel).first;

        auto start = std::chrono::high_resolution_clock::now();
        sort(actual.begin(),actual.end());
        auto end = std::chrono::high_resolution_clock::now();
        sortTime += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        for(int i=0;i<actual.size();i++){ // por cada item dentro del arreglo
            if(actual.at(i).first < 0) continue;
            if(element >= actual.at(i).first){ // comparo el num elementos menores
                rank += pow(2,nivel); // en caso de que existan menores sumar acordemente segun el nivel
            }
        }    
    }

    auto endS = std::chrono::high_resolution_clock::now();
    searchTime += std::chrono::duration_cast<std::chrono::microseconds>(endS - startS).count();
    return pow(2,H_pp)*rank;
}

vector<uint64_t> KLLTuple::rank(vector<int64_t> elements){
    auto start = std::chrono::high_resolution_clock::now();
    // Procedimiento similar a select.
    // 1. Se agregan todos los elementos de manera (Elemento, frecuencia)
    // 2. Se realiza un sort en el vector donde se almacenan
    // 3. Se itera para responder el rank de cada elemento
    // 4. Recordar que el rank se multiplica por 2^H_pp

    vector<uint64_t> ranks;

    vector<pair<pair<int64_t,int64_t>,uint64_t>> vectorElementos = obtenerResumenElementos(); // par(elemento,peso) 

    uint64_t actualPair = 0;
    uint64_t actualRank = 0;
    uint64_t vectorElementosSize = vectorElementos.size();

    // una vez tenemos el vector con los pares ordenados, podemos obtener el rank de los elementos
    for(int i=0;i<elements.size();i++){
        //cerr << i << " " << actualPair << endl;
        while(actualPair < vectorElementosSize && elements.at(i) >= vectorElementos.at(actualPair).first.first){ // comparo el num Elementos menores
            actualRank += vectorElementos.at(actualPair).second; // en caso de que existan menores sumar acordemente segun el nivel
            actualPair++;
        }
        ranks.push_back(pow(2,H_pp)*actualRank);
    }

    if(ranks.size()!=elements.size()){
        ranks.push_back(pow(2,H_pp)*actualRank);
    }

    auto end = std::chrono::high_resolution_clock::now();
    searchTime += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

    return ranks;
}

pair<int64_t,int64_t> KLLTuple::select(uint64_t rank){
    auto start = std::chrono::high_resolution_clock::now();
    vector<pair<pair<int64_t,int64_t>,uint64_t>> elementos = obtenerResumenElementos(); // par(elemento,peso) 
    
    uint64_t rankActual = 0;
    //retornar segun la suma de elementos
    for(int i=0;i<elementos.size();i++){
	    if(elementos.at(i).first.first < 0) continue;
        rankActual+=elementos.at(i).second;
        if(rank<=rankActual) return elementos.at(i).first;
    }

    auto end = std::chrono::high_resolution_clock::now();
    searchTime += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    return elementos.at(elementos.size()-1).first;
}

pair<pair<int64_t,int64_t>, bool> KLLTuple::selectMinMax(uint64_t rank){
    auto start = std::chrono::high_resolution_clock::now();
    /* Funcion equivalente a select pero en el que ARBITRARIAMENTE se revisan las frecuencias 
     * para determinar si el elemento consultado puede ser el minimo o el maximo revisado
     * en dicho caso, además de retornar el elemento en par.first, par.second indicara true
     */
    vector<pair<pair<int64_t,int64_t>,uint64_t>> elementos = obtenerResumenElementos(); // par(elemento,peso)
    bool isMinMax = false;
    pair<int64_t,int64_t> toReturn;
    
    uint64_t rankActual = 0;
    //retornar segun la suma de elementos
    for(int i=0;i<elementos.size();i++){
	    if(elementos.at(i).first.first < 0) continue;
        rankActual+=elementos.at(i).second;
        if(rank<=rankActual) toReturn = elementos.at(i).first;
    }

    if(rank <= elementos.at(0).second/2){
        cerr << "min element, rank consultando:" << rank << ", rank primer elemento:" << elementos.at(0).second << endl;
        toReturn = minElement;
        isMinMax = true;
    } 
    else if (rankActual-rank <= elementos.at(elementos.size()-1).second/2){
        cerr << "max element, rank consultando:" << rank << ", rank ultimo elemento:" << elementos.at(elementos.size()-1).second << endl;
        toReturn = maxElement;
        isMinMax = true;
    } 
    
    auto end = std::chrono::high_resolution_clock::now();
    searchTime += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    return pair<pair<int64_t,int64_t>, bool>(toReturn, isMinMax);
}

vector<pair<int64_t,int64_t>> KLLTuple::select(vector<uint64_t> ranks){
    auto start = std::chrono::high_resolution_clock::now();
    vector<pair<int64_t,int64_t>> selected;

    vector<pair<pair<int64_t,int64_t>,uint64_t>> vectorElementos = obtenerResumenElementos(); // par(elemento,peso) 

    uint64_t actualPair = 0;
    uint64_t actualRank = 0;
    pair<int64_t,int64_t> actualElement = vectorElementos.at(0).first;

    uint64_t vectorElementosSize = vectorElementos.size();
    
    for(int i=0;i<ranks.size();i++){
        //cerr << "i: " <<  i << endl;
        while(actualRank < ranks.at(i) && actualPair < vectorElementosSize){ // comparo el num Elementos menores
            actualRank += vectorElementos.at(actualPair).second; // en caso de que existan menores sumar acordemente segun el nivel
            actualElement = vectorElementos.at(actualPair).first;
            actualPair++;
        }
        selected.push_back(actualElement);
        continue;
        
    }

    while(ranks.size()!=selected.size()){
        ranks.push_back(pow(2,H_pp)*actualRank);
    }

    auto end = std::chrono::high_resolution_clock::now();
    searchTime += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    return selected;
}

vector<pair<pair<int64_t,int64_t>,bool>> KLLTuple::selectMinMax(vector<uint64_t> ranks){
    auto start = std::chrono::high_resolution_clock::now();
    vector<pair<pair<int64_t,int64_t>,bool>> selected;

    vector<pair<pair<int64_t,int64_t>,uint64_t>> vectorElementos = obtenerResumenElementos(); // par(elemento,peso) 

    uint64_t actualPair = 0;
    uint64_t actualRank = 0;
    pair<int64_t,int64_t> actualElement = vectorElementos.at(0).first;

    uint64_t vectorElementosSize = vectorElementos.size();
    
    for(int i=0;i<ranks.size();i++){
        //cerr << "i: " <<  i << endl;
        while(actualRank < ranks.at(i) && actualPair < vectorElementosSize){ // comparo el num Elementos menores
            actualRank += vectorElementos.at(actualPair).second; // en caso de que existan menores sumar acordemente segun el nivel
            actualElement = vectorElementos.at(actualPair).first;
            actualPair++;
        }
        selected.push_back(pair<pair<int64_t,int64_t>,bool>(actualElement,false));
        continue;
        
    }

    while(ranks.size()!=selected.size()){
        ranks.push_back(pow(2,H_pp)*actualRank);
    }
    
    for(int i=0;i<ranks.size();i++){
        if(ranks.at(i)>vectorElementos.at(0).second/2) break;
        selected.at(i).first = minElement;
        selected.at(i).second = true;
    }

    for(int i=ranks.size()-1;i>=0;i--){
        cerr << "fin" << endl;
        cerr << "ranks size " << ranks.size() << " selected size " << selected.size() << " vectorElementosSize " << vectorElementosSize << " vector.size " << vectorElementos.size() << endl;
        if(actualRank-ranks.at(i)>vectorElementos.at(vectorElementosSize-1).second/2) break;
        selected.at(i).first = maxElement;
        selected.at(i).second = true;
    }

    auto end = std::chrono::high_resolution_clock::now();
    searchTime += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    return selected;
}

pair<int64_t,int64_t> KLLTuple::quantile(double q){ // q pertenece [0,100]
    q = q/100.0;
    return select(floor(q*numElementosRevisados));
}

vector<pair<int64_t,int64_t>> KLLTuple::quantile(vector<double> q){
    vector<uint64_t> ranks; 
    for(int i=0;i<q.size();i++)
        ranks.push_back(floor(q[i]/100.0 * numElementosRevisados));
    return select(ranks);
}

long KLLTuple::height(){
    return numArreglos;
}

void KLLTuple::print(){
    string name = isMrl?"MRL":"KLL";
    string limSpace = espacioLimitado? "Limitado": "No Limitado";
    string cteSpace = espacioCte? "Cte": "No Cte";
    cout << name << ", Espacio " << limSpace << " " << cteSpace << ".  " << sizeInBytes()/1024 << "[KB]" << endl;
    cout << "H: " << H << ", s: " << s << ", H'': " << H_pp << endl;
    cout << "Sampler Element: (" << sampleElement.first << ","<< sampleElement.second <<")" << ", Sample Weight: " << sampleWeight << endl; 
    cout << "numElementosRevisados: " << numElementosRevisados << " numTotal: " << numTotalElementos << endl;
    for(int i=0; i<sketch.size();i++){
        cout << "Nivel " << i+H_pp << ": (" << sketch.at(i).second << "/" << sketch.at(i).first.size() << ")" << endl;
        vector<pair<int64_t,int64_t>> nivelI = sketch.at(i).first;
        for(int j=0;j<nivelI.size();j++){
            printf("(%ld,%ld) ",nivelI.at(j).first,nivelI.at(j).second);
            //cout << nivelI.at(j) << " ";
        }
        cout << endl;
    }
    heap.printMinHeap();

    vector<pair<int64_t,int64_t>> heapElem = heap.getSortHeap();
    for(int i=0;i<heapElem.size();i++){
        printf("(%ld,%ld)\n",heapElem.at(i).first,heapElem.at(i).second);
    }
}

pair<vector<pair<int64_t,int64_t>>, long> KLLTuple::sketchAtLevel(long nivel){
    if(nivel<0||nivel>=numArreglos) return make_pair(vector<pair<int64_t,int64_t>>(),-1);
    return sketch.at(nivel);
}

void KLLTuple::addToSampler(pair<int64_t,int64_t> element,uint64_t weight){
    if(!reservoirKLLSample(element,weight)) return;
    insertElement(0,sampleElement);
    return;
}

pair<uint64_t,uint64_t> KLLTuple::getNumElementsSeen(){
    pair<uint64_t,uint64_t> toReturn;
    toReturn.first = numElementosRevisados;
    toReturn.second = numTotalElementos;
    return toReturn;
}

uint64_t KLLTuple::numElementosRango(pair<int64_t,int64_t> a, pair<int64_t,int64_t> b){
    vector<pair<pair<int64_t,int64_t>,uint64_t>> vectorElementos = obtenerResumenElementos();; // par(elemento,peso) 

    uint64_t numElementosRepetidos = 0;
    for(int i=0;i<vectorElementos.size();i++){
        if(vectorElementos.at(i).first.first < a.first) continue;
        else if (vectorElementos.at(i).first.first > b.first) break;
        numElementosRepetidos+=vectorElementos.at(i).second;
    }

    return numElementosRepetidos*pow(2,H_pp);
}

void KLLTuple::update(KLLTuple kll2){
    // en caso de que algun KLLTuple tenga espacioLimitado se realizara la compactación necesaria,
    // esto se hara de arriba hacia abajo con el fin de evitar choques
    cout << "KLL1 PRE MERGE:" << endl;
    print();
    cout << "KLL2 PRE MERGE:" << endl;
    kll2.print();
    cout << endl;

    if(hasLimitedSpace()){
        for(int i=sketch.size()-1;i>=0;i--){
            if(sketch.at(i).first.size()==sketch.at(i).second) iterativeCompaction(i,true);
        }
    }

    if(kll2.hasLimitedSpace()){
        for(int i=kll2.getH()-kll2.getH_pp();i>=0;i--){
            if(kll2.sketchAtLevel(i).first.size()==kll2.sketchAtLevel(i).second){
                kll2.iterativeCompaction(i,true);
            }
        }
    }

    pair<pair<int64_t,int64_t>,pair<int64_t,int64_t>> minMaxOtherKLL = kll2.getMinMaxElement();
    if(minElement<minMaxOtherKLL.first) minElement = minMaxOtherKLL.first;
    if(maxElement>minMaxOtherKLL.second) maxElement = minMaxOtherKLL.second;
    pair<uint64_t,uint64_t> numElementsOtherKLL = kll2.getNumElementsSeen();
    numElementosRevisados += numElementsOtherKLL.first;
    numTotalElementos += numElementsOtherKLL.second;

    // ACORDARSE DE QUE kll1.height()>=kll2.height()
    if(isMrl){ // KLL 1 es Mrl
        if(kll2.isAnMrl()){ // KLL 2 es MRL
            for(int nivel=0; nivel<kll2.height();nivel++){ // simplemente mover de KLL2 a KLL1
                pair<vector<pair<int64_t,int64_t>>,long> kll2pair = kll2.sketchAtLevel(nivel);
                for(int i=0; i < kll2pair.second;i++){
                    if(kll2pair.first.at(i).first<0) continue;
                    insertCompactionElement(nivel,kll2pair.first.at(i),true);
                }
            }
        }
        else{ // KLL 2 es KLL Tradicional
            pair<unsigned long,pair<int64_t,int64_t>> kll2sample = kll2.getCurrentSample();
            unsigned long kll2SampleWeight = kll2sample.first;

            //insert elementos del sampler a mrl
            if (kll2SampleWeight != 0) { // si es 0 no hay que ingresar nada
                // se calculara el nivel correspondiente en el que insertar el elemento
                // para ello, se calculara el nivel al que deberia ingresarse segun el peso del sampler
                // en el que se tomara como relevancia el peso del sampler y de los niveles asociados

                kll2SampleWeight = kll2sample.first;
                uint64_t nivelSuperior=0; // nivel inferior siempre es (nivelSuperior-1)
                while (kll2SampleWeight /2 != 0) {
                    kll2SampleWeight /= 2;
                    nivelSuperior++;
                }
                uint64_t pesoNivelSuperior = 1<<nivelSuperior;
                uint64_t pesoNivelInferior = 1<<(nivelSuperior-1);
                
                double probNivelSuperior = (double) (kll2sample.first-pesoNivelInferior)/(double) (pesoNivelSuperior-pesoNivelInferior);
                if((double) (rand()/RAND_MAX) <= probNivelSuperior) insertCompactionElement(nivelSuperior, kll2sample.second, true);
                else insertCompactionElement(nivelSuperior-1, kll2sample.second, true);
            }
            
            //insert elementos de kll a mrl en el nivel respectivo
            for(int nivel=kll2.getH_pp();nivel<kll2.height(); nivel++){
                pair<vector<pair<int64_t,int64_t>>,long> kll2pair = kll2.sketchAtLevel(nivel);
                for(int i=0; i < kll2pair.second;i++){
                    if(kll2pair.first.at(i).first<0) continue;
                    insertCompactionElement(nivel,kll2pair.first.at(i),true);
                }
            }
        }
    } else { // KLL 1 es KLL Tradicional
        if(kll2.isAnMrl()){ // KLL 2 es un MRL
            // ingresamos elementos asociados al sampler
            int nivel = 0;
            for(nivel; nivel<H_pp;nivel++){
                if(nivel>=kll2.height()) break; // ya se ingresaron todos los niveles
                pair<vector<pair<int64_t,int64_t>>,long> kll2pair = kll2.sketchAtLevel(nivel);
                uint64_t pesoNivel = pow(2,nivel);
                for(int i=0; i < kll2pair.second;i++){
                    if(kll2pair.first.at(i).first<0) continue;
                    addToSampler(kll2pair.first.at(i),pesoNivel);
                }
            }

            // ingresamos el resto de los elementos al compactor asociado
            for(nivel;nivel<kll2.height(); nivel++){
                pair<vector<pair<int64_t,int64_t>>,long> kll2pair = kll2.sketchAtLevel(nivel);
                for(int i=0; i < kll2pair.second;i++){
                    if(kll2pair.first.at(i).first<0) continue;
                    insertCompactionElement(nivel-H_pp,kll2pair.first.at(i),true);
                }
            }
        }
        else{ // KLL 2 es KLL tradicional
            int nivel = 0;
            pair<unsigned long,pair<int64_t,int64_t>> kll2sample = kll2.getCurrentSample();
            unsigned long kll2SampleWeight = kll2sample.first;

            //insert elementos del sampler a kll
            if(kll2SampleWeight == 0){} // si es 0 no hay que ingresar nada
            else if(wH_pp < kll2SampleWeight) {  // caso en que el peso es superior a H_pp (se agrega el elemento a un compactor)
                kll2SampleWeight = kll2sample.first;
                uint64_t nivelSuperior=0; // nivel inferior siempre es (nivelSuperior-1)
                while (kll2SampleWeight /2 != 0) {
                    kll2SampleWeight /= 2;
                    nivelSuperior++;
                }
                uint64_t pesoNivelSuperior = 1<<nivelSuperior;
                uint64_t pesoNivelInferior = 1<<(nivelSuperior-1);
                
                double probNivelSuperior = (double) (kll2sample.first-pesoNivelInferior)/(double) (pesoNivelSuperior-pesoNivelInferior);
                if((double) (rand()/RAND_MAX) <= probNivelSuperior) insertCompactionElement(nivelSuperior-H_pp, kll2sample.second, true);
                else insertCompactionElement(nivelSuperior-1-H_pp, kll2sample.second, true);
            } else { // aqui se agrega el elemento desde kll2_sampler a kll1_sampler
                addToSampler(kll2sample.second,kll2sample.first);
            }
            
            int nivelesEnSampler=0;
            int diffH_pp = H_pp - kll2.getH_pp();
            if(H_pp>kll2.getH_pp()){
                for(nivel=0; nivel<diffH_pp;nivel++){
                    if(nivel>=kll2.height()) break; //
                    pair<vector<pair<int64_t,int64_t>>,long> kll2pair = kll2.sketchAtLevel(nivel);
                    uint64_t pesoNivel = pow(2,nivel+kll2.getH_pp());
                    for(int i=0; i < kll2pair.second;i++){
                        if(kll2pair.first.at(i).first<0) continue;
                        addToSampler(kll2pair.first.at(i),pesoNivel);
                        iterativeCompaction(0,true);
                    }
                    nivelesEnSampler++;
                }
            }

            for(nivel; nivel<kll2.height();nivel++){
                pair<vector<pair<int64_t,int64_t>>,long> kll2pair = kll2.sketchAtLevel(nivel);
                //! agregar metodo para agregar varios elementos de una
                for(int i=0; i < kll2pair.second;i++){
                    if(kll2pair.first.at(i).first<0) continue;
                    insertCompactionElement(nivel-nivelesEnSampler,kll2pair.first.at(i),true);
                }
            }
        }
    }
    
    return;
}

KLLTuple KLLTuple::kllMerge(KLLTuple &kll2){ 
    if(isMrl!=kll2.isAnMrl()){
        cerr << "Los KLL son incompatibles" << endl;
        return KLLTuple(200);
    }

    uint64_t heightKLL1, heightKLL2;
    heightKLL1 = H;
    heightKLL2 = kll2.getH();

    // Creamos una copia del kll a devolver
    if(isMrl){
        cerr << "isMrl" << endl;
        // Se verifica el MRL que tenga mayor altura para definir cual utilizar como base
        // en caso de misma altura se definira mediante cual tenga mayor minK
        if(heightKLL1 < heightKLL2){
            cerr << "kll 2 has higher height" << endl;
            KLLTuple kllCopy2 = kll2.copy();
            kllCopy2.update(*this);
            return kllCopy2;
        } 
        else if(heightKLL1 > heightKLL2){
            cerr << "kll 1 has higher height" << endl;
            KLLTuple kllCopy1 = copy();
            kllCopy1.update(kll2);
            return kllCopy1;
        } else if(getFirstLevelK() < kll2.getFirstLevelK()){
            cerr << "same height and kll 2 has higher minK" << endl;
            KLLTuple kllCopy2 = kll2.copy();
            kllCopy2.update(*this);
            return kllCopy2;
        } else {
            cerr << "same height and kll 1 has higher or same minK" << endl;
            KLLTuple kllCopy1 = copy();
            kllCopy1.update(kll2);
            return kllCopy1;
        }
    } else {
        cerr << "isKLL" << endl;
        if(heightKLL1 < heightKLL2){
            KLLTuple kllCopy2 = kll2.copy();
            kllCopy2.update(*this);
            return kllCopy2;
        } 
        else{
            KLLTuple kllCopy1 = copy();
            kllCopy1.update(kll2);
            return kllCopy1;
        }
    }
}


int KLLTuple::getHeapCapacity(){
    return heap.getCapacity();
}

int KLLTuple::getHeapLevels(){
    return heap.getLevels();
}


vector<pair<int64_t,int64_t>> KLLTuple::getSortHeap(){
    return heap.getSortHeap();
}

vector<pair<int64_t,int64_t>> KLLTuple::getTopFlows(){
    // obtiene los flujos encontrados en el heap, y estima los payloads segun valores encontrados 
    // tanto en heap como en la estructura del kll, si se desea ahorrar tiempo se puede ignorar paso
    // de busqueda de elementos en el sketch y simplemente proporcionar lo que el heap nos proporciona

    vector<pair<int64_t,int64_t>> flows;
    vector<pair<int64_t,int64_t>> heap = getSortHeap();
    vector<pair<int64_t,int64_t>> heapFlows; // vector que tiene los flujos más relevantes del minHeap
    vector<pair<int64_t,int64_t>> kllFlows;
    
    // Calculo de los flujos más relevantes e insercion en heapFlows
    double factorReduccion = 1;
    heapFlows.push_back(heap.at(0));
    for(int i=0;i<heap.size()-1;i++){
        //if(heap.at(i).first/factorReduccion<=heap.at(i+1).first) heapFlows.push_back(heap.at(i+1));
        heapFlows.push_back(heap.at(i+1));
    }

    // obtener los pares(payload, flujo) que pertenecen al 95%percentil de payloads en kll
    double quantil = 0.95;
    uint64_t rankKll = numElementosRevisados*quantil;
    vector<pair<pair<int64_t,int64_t>,uint64_t>> vectorElementos = obtenerResumenElementos(); // par(elemento,peso) 
    uint64_t actualRank = 0;
    for(int i=0;i<vectorElementos.size();i++){
        actualRank += vectorElementos.at(i).second; // actualizarRank
        if(actualRank >= rankKll){ // se obtienen los elementos que se encuentran sobre el quantil indicado
            kllFlows.push_back(vectorElementos.at(i).first);
        }
    }

    for (pair<int64_t,int64_t> &flow : kllFlows) {
        flow.first = flow.first*pow(2,H_pp);
    }
    for (pair<int64_t,int64_t> &flow : heapFlows) {
        flow.first = flow.first*pow(2,H_pp);
    }

    flows = heapFlows;
    // Iterar sobre los elementos de kllFlows
    for (const auto& flow : kllFlows) {
        // Buscar si ya existe un elemento con el mismo pair.second en flows
        auto it = std::find_if(flows.begin(), flows.end(), [&flow](const auto& f) {
            return f.second == flow.second;
        });

        // Si existe, sumar los primeros componentes
        if (it != flows.end()) {
            it->first += flow.first;
        } else {
            // Si no existe, agregar el elemento a flows
            flows.push_back(flow);
        }
    }

    return flows;
}

vector<pair<int64_t,int64_t>> KLLTuple::getTopFlows(vector<pair<int64_t,int64_t>> &heapFlows,vector<pair<int64_t,int64_t>> &kllFlows){
    // obtiene los flujos encontrados en el heap, y estima los payloads segun valores encontrados 
    // tanto en heap como en la estructura del kll, si se desea ahorrar tiempo se puede ignorar paso
    // de busqueda de elementos en el sketch y simplemente proporcionar lo que el heap nos proporciona

    vector<pair<int64_t,int64_t>> flows;
    vector<pair<int64_t,int64_t>> heap = getSortHeap();
    
// Calculo de los flujos más relevantes e insercion en heapFlows
    // double factorReduccion = 5;
    // heapFlows.push_back(heap.at(0));
    // for(int i=0;i<heap.size()-1;i++){
    //     cout << "comp: " << heap.at(i).first/factorReduccion <<"|" << heap.at(i+1).first << endl;
    //     if(heap.at(i).first/factorReduccion<=heap.at(i+1).first) heapFlows.push_back(heap.at(i+1));
    // }

    // Calculo de los flujos más relevantes e insercion en heapFlows
    double factorReduccion = 0;
    heapFlows.push_back(heap.at(0));
    for(int i=0;i<heap.size()-1;i++){
        //if(heap.at(i).first*factorReduccion<=heap.at(i+1).first) heapFlows.push_back(heap.at(i+1));
        heapFlows.push_back(heap.at(i+1));
    }

    // obtener los pares(payload, flujo) que pertenecen al 95%percentil de payloads en kll
    map<int64_t,int64_t> mapaKll;
    double quantil = 0.95;
    // cout << "numRevisados: " << numElementosRevisados << ", numTotal: " << numTotalElementos << endl;
    uint64_t rankKll = numElementosRevisados*quantil;
    vector<pair<pair<int64_t,int64_t>,uint64_t>> vectorElementos = obtenerResumenElementos(); // par(elemento,peso)
    uint64_t actualRank = 0;
    for(int i=0;i<vectorElementos.size();i++){
        actualRank += vectorElementos.at(i).second; // actualizarRank
        // cout << actualRank << "/" << rankKll <<endl;
        if(actualRank >= rankKll){ // se obtienen los elementos que se encuentran sobre el quantil indicado
            //kllFlows.push_back(vectorElementos.at(i).first);
            mapaKll[vectorElementos.at(i).first.second] = mapaKll[vectorElementos.at(i).first.second] + vectorElementos.at(i).first.first;
        }
    }

    for (const auto& par : mapaKll) {
        kllFlows.push_back(std::make_pair(par.second, par.first));
    }


    for (pair<int64_t,int64_t> &flow : kllFlows) {
        flow.first = flow.first*pow(2,H_pp);
    }
    for (pair<int64_t,int64_t> &flow : heapFlows) {
        flow.first = flow.first*pow(2,H_pp);
    }

    flows = heapFlows;
    // Iterar sobre los elementos de kllFlows
    for (const auto& flow : kllFlows) {
        // Buscar si ya existe un elemento con el mismo pair.second en flows
        auto it = std::find_if(flows.begin(), flows.end(), [&flow](const auto& f) {
            return f.second == flow.second;
        });

        // Si existe, sumar los primeros componentes
        if (it != flows.end()) {
            it->first += flow.first;
        } else {
            // Si no existe, agregar el elemento a flows
            // EN CASO DE QUERER INGRESAR ELEMENTOS UNICOS DEL KLL DESCOMENTAR 
            // flows.push_back(flow);
        }
    }

    return flows;
}

vector<pair<int64_t,int64_t>> KLLTuple::getTopHeapFlows(){
    // obtiene los flujos encontrados en el heap, y obtiene los payloads asociado a los flujos encontrados
    // en el heap
    vector<pair<int64_t,int64_t>> flows = getSortHeap();

    // for(int i=0;i<sketch.size();i++){ // busqueda en los niveles del sketch
    //     for(int j=0;j<sketch.at(i).second;j++){ // en cada nivel se revisa todos los elementos ingresados
    //         for(int k=0;k<flows.size();k++){ // search del elemento en el heap
    //             if(sketch.at(i).first.at(j).second == flows.at(k).second) // si los flujos son iguales
    //                 flows.at(k).first += sketch.at(i).first.at(j).first; // se suma el payload
    //                 break; // se sigue al siguiente elemento del nivel
    //         }
    //     }
    // }

    return flows;
}

pair<unsigned long, pair<int64_t,int64_t>> KLLTuple::getCurrentSample(){
    pair<unsigned long, pair<int64_t,int64_t>> toReturn;
    toReturn.first = sampleWeight;
    toReturn.second = sampleElement;
    return toReturn;
}

uint32_t KLLTuple::getH_pp(){
    return H_pp;
}

uint32_t KLLTuple::getH(){
    return H;
}

uint32_t KLLTuple::getFirstLevelK(){
    return sketch.at(0).first.size();
}

bool KLLTuple::isAnMrl(){
    return isMrl;
}

pair<pair<int64_t,int64_t>, pair<int64_t,int64_t>> KLLTuple::getMinMaxElement(){
    pair<pair<int64_t,int64_t>,pair<int64_t,int64_t>> toReturn;
    toReturn.first = minElement;
    toReturn.second = maxElement;
    return toReturn;
}

bool KLLTuple::areEqual(KLLTuple toCompare){
    if(isMrl != toCompare.isAnMrl()) return false;
    // cerr << "ISMRL " << isMrl << endl;
    if(H_pp!=toCompare.getH_pp()) return false;
    if(H!=toCompare.getH()) return false;
    // cerr << "DATOS BASICOS IGUALES " << isMrl << endl;
    for(int i=0;i<numArreglos;i++){
        // cerr << "Nivel: " << i << endl;
        pair<vector<pair<int64_t,int64_t>>,long> toCompareLvl = toCompare.sketchAtLevel(i);
        pair<vector<pair<int64_t,int64_t>>,long> actualLvl = sketchAtLevel(i);
        if(actualLvl.first.size()!=toCompareLvl.first.size()) return false;
        else if (actualLvl.second!=toCompareLvl.second) return false;
        for(int j=0;j<actualLvl.first.size();j++){
            if(actualLvl.first.at(j)==make_pair((int64_t)-1,(int64_t)-1)) actualLvl.first.at(j)=make_pair(-2,-2);
            if(toCompareLvl.first.at(j)==make_pair((int64_t)-1,(int64_t)-1)) toCompareLvl.first.at(j)=make_pair(-2,-2); 
            if(actualLvl.first.at(j)!=toCompareLvl.first.at(j)) return false;
        }
    }
    return true;
}


uint64_t KLLTuple::getHeapSize(){
    return heap.sizeInBytes();
}

uint64_t KLLTuple::sizeInBytes(){
    // cout << " " << sizeof(pair<int64_t,int64_t>) << " ";
    uint64_t totalSize = 0;
    uint64_t sketchSize = 0;
    uint64_t heapSize = 0;

    //calculo del tamaño de sketchSize
    for(int i = 0; i< sketch.size();i++){
        sketchSize+=sizeof(sketch.at(i).second); // entero que me indica numero de valores ocupados en el nivel i
        // cout << endl  <<(sketch.at(i).first.size()*sizeof(pair<int64_t,int64_t>)) << endl; 
        sketchSize += (sketch.at(i).first.size()*sizeof(pair<int64_t,int64_t>));
    }
    heapSize = heap.sizeInBytes();
    // cout << sketchSize << endl;

    //calculo de totalSize según las variables utilizadas
    totalSize+=sketchSize;
    totalSize+=heapSize;
    totalSize+=sizeof(hashLong);
    totalSize+=sizeof(sampleWeight);
    totalSize+=sizeof(sampleElement);
    totalSize+=sizeof(numArreglos);
    totalSize+=sizeof(n);
    totalSize+=sizeof(numElementosRevisados);
    totalSize+=sizeof(numTotalElementos);
    totalSize+=sizeof(epsilon);
    totalSize+=sizeof(delta);
    totalSize+=sizeof(c);
    totalSize+=sizeof(H);
    totalSize+=sizeof(H_p);
    totalSize+=sizeof(H_pp);
    totalSize+=sizeof(k);
    totalSize+=sizeof(minK);
    totalSize+=sizeof(s);
    totalSize+=sizeof(mascara);
    totalSize+=sizeof(wH_pp);
    totalSize+=sizeof(metodoCompactacion);
    totalSize+=sizeof(debug);
    totalSize+=sizeof(isMrl);
    totalSize+=sizeof(espacioLimitado);
    totalSize+=sizeof(minElement);
    totalSize+=sizeof(maxElement);

    return totalSize;
}

vector<double> KLLTuple::parametros(){
    vector<double> parametros;
    parametros.push_back(numTotalElementos);
    parametros.push_back(epsilon);
    parametros.push_back(delta);
    parametros.push_back(c);
    parametros.push_back(H);
    parametros.push_back(H_p);
    parametros.push_back(H_pp);
    parametros.push_back(s);
    parametros.push_back(k);
    parametros.push_back(sizeInBytes());
    parametros.push_back(minK);
    parametros.push_back(isMrl);
    return parametros;
}

bool KLLTuple::hasLimitedSpace(){
    return espacioLimitado;
}

bool KLLTuple::hasConstantSpace(){
    return espacioCte;
}


vector<pair<uint64_t,pair<int64_t,int64_t>>> KLLTuple::obtenerResumenFrecuencias(){ // revisar si existe espacio para optimizacion en la obtencion del map
    // procedimiento similar al de obtenerResumenElementos, pero los ordena
    // por frecuencia en vez de por elemento. Es importante indicar que se deben 
    // juntar los elementos repetidos

    vector<pair<pair<int64_t,int64_t>,uint64_t>> elementos = obtenerResumenElementos();
    vector<pair<uint64_t,pair<int64_t,int64_t>>> elementosMasFrecuentes; // par(peso,elementos) 
    pair<pair<int64_t,int64_t>,uint64_t> primerPar = elementos.at(0);
    elementosMasFrecuentes.push_back(make_pair(primerPar.second,primerPar.first));
    for(int i=1;i<elementos.size();i++){
        pair<pair<int64_t,int64_t>,uint64_t> parActual = elementos.at(i);
        if(parActual.first==elementosMasFrecuentes.at(elementosMasFrecuentes.size()-1).second)
            elementosMasFrecuentes.at(elementosMasFrecuentes.size()-1).first += parActual.second;
        else
            elementosMasFrecuentes.push_back(make_pair(parActual.second,parActual.first));
    }

    // realizar el sort
    sort(elementosMasFrecuentes.begin(),elementosMasFrecuentes.end());
    return elementosMasFrecuentes;
}

vector<pair<int64_t,int64_t>> KLLTuple::topKElements(uint32_t kParam){
    vector<pair<int64_t,int64_t>> topKElementos;

    vector<pair<uint64_t,pair<int64_t,int64_t>>> elementosMasFrecuentes = obtenerResumenFrecuencias();
    uint32_t posMasAltoFrecuente = elementosMasFrecuentes.size()-1;

    for(int i=0;i<kParam;i++){
        if(i==elementosMasFrecuentes.size()) break;
        topKElementos.push_back(elementosMasFrecuentes.at(posMasAltoFrecuente-i).second);
    }

    return topKElementos;
}

vector<pair<pair<int64_t,int64_t>,uint64_t>> KLLTuple::topKElementsFrequency(uint32_t kParam){ 
    vector<pair<pair<int64_t,int64_t>,uint64_t>> topKElementos;

    vector<pair<uint64_t,pair<int64_t,int64_t>>> elementosMasFrecuentes = obtenerResumenFrecuencias();
    uint32_t posMasAltoFrecuente = elementosMasFrecuentes.size()-1;

    for(int i=0;i<kParam;i++){
        if(i==elementosMasFrecuentes.size()) break;
        topKElementos.push_back(make_pair(elementosMasFrecuentes.at(posMasAltoFrecuente-i).second,(1<<H_pp)*elementosMasFrecuentes.at(posMasAltoFrecuente-i).first));
    }

    return topKElementos;
}

vector<pair<pair<int64_t,int64_t>,uint64_t>> KLLTuple::estimateElementsFrequency(){
    vector<pair<pair<int64_t,int64_t>,uint64_t>> topKElementos;

    vector<pair<uint64_t,pair<int64_t,int64_t>>> elementosMasFrecuentes = obtenerResumenFrecuencias();
    uint32_t posMasAltoFrecuente = elementosMasFrecuentes.size()-1;

    for(int i=0;i<elementosMasFrecuentes.size();i++){
        topKElementos.push_back(make_pair(elementosMasFrecuentes.at(posMasAltoFrecuente-i).second,(1<<H_pp)*elementosMasFrecuentes.at(posMasAltoFrecuente-i).first));
    }

    return topKElementos;
}

// returna el nombre del binario a almacenar, 
// si se le proporciona un nombre de traza formatea el nombre con este nombre
string KLLTuple::binarySaveName(string archivoTraza){
    string toReturn = "";
    if(archivoTraza!=""){
        //quitar.txt
        size_t posicion = archivoTraza.find(".txt");
        if (posicion != std::string::npos) {
            archivoTraza.erase(posicion, 4); // Eliminamos 4 caracteres: txt
        }
        toReturn+=archivoTraza+".";
    }
    int numNivelesHeap = heap.getLevels();
    if(isMrl){
        toReturn+="mrltuplek" + to_string(minK)+"h"+to_string(numNivelesHeap);
    }else{
        toReturn+="klltuplee"+to_string(epsilon)+"d"+to_string(delta)+"c"+to_string(c)+"mk"+to_string(minK)+"h"+to_string(numNivelesHeap);
    }

    return toReturn;
}

vector<double> KLLTuple::getMRLParameters(){
    vector<double> toReturn;
    toReturn.push_back(k);
    return toReturn;
}

vector<double> KLLTuple::getKLLParameters(){
    vector<double> toReturn;
    toReturn.push_back(epsilon);
    toReturn.push_back(delta);
    toReturn.push_back(c);
    toReturn.push_back(minK);
    return toReturn;
}

vector<double> KLLTuple::getAllParameters(){
    vector<double> toReturn;
    toReturn.push_back(H);
    toReturn.push_back(H_pp);
    toReturn.push_back(epsilon);
    toReturn.push_back(delta);
    toReturn.push_back(c);
    toReturn.push_back(minK);
    toReturn.push_back(numArreglos);
    for(int i=0;i<sketch.size();i++){
        //cout << i << " " << numArreglos << endl;
        toReturn.push_back(sketch.at(i).first.size());
    }
    return toReturn;
}

vector<uint64_t> KLLTuple::getTimes(){
    vector<uint64_t> times;
    times.push_back(sortTime);
    times.push_back(buildTime);
    times.push_back(searchTime);
    return times;
}