#include <iostream>
#include <cmath>
#include <algorithm>
#include <stdio.h>
#include <string.h>
#include <chrono>
#include <fstream>
#include <map>

#include "kll.hpp"

void debugLetra(string string){
    cout << string << endl;
}

using namespace std;

void KLL::setupKLL(uint64_t nP, double epsilonP, double deltaP, double cP, int minkP){
    cerr << deltaP << endl;
    isMrl = false;
    n = nP;
    epsilon = epsilonP;
    delta = deltaP;
    c = cP;
    minK = minkP;

    minElement = std::numeric_limits<double>::max();
    maxElement = -1*minElement;
    sampleElement=0;
    sampleWeight=0;

    numElementosRevisados = 0;
    numTotalElementos = 0;


    H = 1.0*log2(epsilon*n); // O(log(epsilon*n))
    k = 1.0*(double)(1.0/epsilon) * log2(log2((1.0/delta)));
    s = log2(log2((1.0/delta)));

    H_p = H-s;
    //H_pp = H - ceil((double)log(k) / (double)log((1/c)) ); // H_pp = H - log(k)
    H_pp = H - ceil((double)log2((1.0/epsilon)) ); // H_pp = H - log(k)

    if(H_pp<0) H_pp = 0;

    wH_pp = pow(2,H_pp);
    mascara = pow(2,H_pp);
    
    numArreglos = (H - H_pp+1);
}

void KLL::setupMRL(int minkP){
    isMrl = true;
    n = 0;
    minK = minkP;
    k = minK;
    epsilon = -1;
    delta = -1;
    c = 1;

    minElement = std::numeric_limits<double>::max();
    maxElement = -1*minElement;
    sampleElement=0;
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

void KLL::startSketch(){ // necesita haberse inicializado k, mink, c, H y H_pp
    vector<long> sizeTemp;
    for(int i=H_pp;i<=H;i++){
        uint64_t cantElementos;
        if(i>(H-s)) cantElementos = max(k,(long)minK);
        else cantElementos = max((int)(k*pow(c,H-i-1)),(int)minK);
        // cout << i << " El: " << cantElementos << endl;
        //cout << "cantElementos en arreglo " << i+1 << ": " << cantElementos << endl;
        sizeTemp.push_back(cantElementos);
    }

    // inicializar los vectores de tam k*c^lvl
    for(int i=0;i<sizeTemp.size();i++){
        // el valor por defecto es -1, que indica "vacio"
        uint64_t cantElementos = sizeTemp.at(i);
        //cerr << "Cantidad Elementos arreglo " << i << " (nivel " << i+H_pp+1 <<") :" << cantElementos<< endl;
        double valorElemento = -2;
        vector<double> vectorAtLvlI(cantElementos,valorElemento); 
        pair<vector<double>,long> toInsert;
        toInsert.first=vectorAtLvlI;
        toInsert.second=0; // representa el num de elementos ocupados en el arreglo
        sketch.push_back(toInsert);
    }
}

void KLL::startLimitedSketchMRL(uint64_t np, uint64_t espacioMaximo){
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
        double valorElemento = -2;
        vector<double> vectorAtLvlI(cantElementos,valorElemento); 
        pair<vector<double>,long> toInsert;
        toInsert.first=vectorAtLvlI;
        toInsert.second=0; // representa el num de elementos ocupados en el arreglo
        sketch.push_back(toInsert);
    }

    H = numNivelesNecesarios; 
    numArreglos = H+1;
}

void KLL::startLimitedSketchKLL(uint64_t np, uint64_t espacioMaximo){
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
        double valorElemento = -2;
        vector<double> vectorAtLvlI(cantElementos,valorElemento); 
        pair<vector<double>,long> toInsert;
        toInsert.first=vectorAtLvlI;
        toInsert.second=0; // representa el num de elementos ocupados en el arreglo
        sketch.push_back(toInsert);
    }

    cout << endl;
}
        
// kll sin espacio limitado
KLL::KLL(uint64_t numElementsParam, double epsilonParam, double deltaParam, double cParam, int minKp){
    espacioCte = false;
    espacioLimitado = false;

    setupKLL(numElementsParam, epsilonParam, deltaParam, cParam, minKp);
    startSketch();
}

KLL::KLL(uint64_t numElementsParam, double epsilonParam, double cParam, int minKp){ 
    // KLL Tradicional con espacio determinado por epsilon, delta=0.01
    espacioCte = false;
    espacioLimitado = true;

    setupKLL(numElementsParam, epsilonParam, 0.01, cParam, minKp);

    //! posteriormente cambiar comportamiento con startSketch y 
    startLimitedSketchKLL(numElementsParam, 15000);
}

KLL::KLL(uint64_t numElementsParam, double epsilonParam, double deltaParam, double cParam, int minKp, bool samplerEnMenorNivel){
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
} // KLL con espacio constante, variable "samplerEnMenorNivel" indica que H'' debe ser 0
        
KLL::KLL(uint64_t minKP){ // MRL sin espacioLimitado
    espacioLimitado = false;
    espacioCte = false;
    setupMRL(minKP);

    startSketch();
}

KLL::KLL(uint64_t espacioMaximo, uint64_t numElementsParam){ // MRL que se le entrega el espacio
    espacioLimitado = true;
    espacioCte = false;
    cout << "Espacio Maximo: " << espacioMaximo << " NumElementosAAbarcar: " << numElementsParam << endl;
    
    setupMRL(2);
    startLimitedSketchMRL(numElementsParam, espacioMaximo);
}

KLL::KLL(double epsilonParam ,uint64_t numElementsParam){
    int minkP = (1.0/epsilonParam)*ceil(log2(epsilonParam*numElementsParam))+1;
    espacioLimitado = false;
    espacioCte = false;
    setupMRL(minkP);

    startSketch();
    /*
    espacioCte = false;
    double cotaEspacio = 1.0/epsilonParam * log2(1.0/epsilonParam) * log2(1.0/epsilonParam);
    uint64_t espacioMaximo = ceil(cotaEspacio);
    cout << "Epsilon Entregado: " << epsilonParam << " Espacio Maximo: " << espacioMaximo << " NumElementosAAbarcar: " << numElementsParam << endl;
    //desde aqui es el mismo constructor que  KLL(uint64_t espacioMaximo, uint64_t numElementsParam)
    espacioLimitado = true;
    isMrl = true;
    epsilon = epsilonParam;
    
    setupMRL(2);

    startLimitedSketchMRL(numElementsParam,espacioMaximo);

    sizeInBytes();
    */
}

// creacion KLL con parametros indicados
KLL::KLL(uint64_t HParam,uint64_t sParam,uint64_t H_ppParam,vector<int> sizeNivelesParam, bool isAnMrlParam, bool hasLimitedSpaceParam){
    isMrl = isAnMrlParam;
    espacioLimitado = hasLimitedSpaceParam;
    espacioCte = false;
    H = HParam;
    s = sParam;
    H_p = H-sParam;
    H_pp = H_ppParam;
    wH_pp = pow(2,H_pp);

    mascara = pow(2,H_pp);
    sampleElement=0;
    sampleWeight=0;
    minElement = std::numeric_limits<double>::max();
    maxElement = -1*minElement;

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
        double valorElemento = -2;
        vector<double> vectorAtLvlI(cantElementos,valorElemento); 
        pair<vector<double>,long> toInsert;
        toInsert.first=vectorAtLvlI;
        toInsert.second=0; // representa el num de elementos ocupados en el arreglo
        sketch.push_back(toInsert);
    }
} // MRL/KLL que se le indican los parametros que va a ocupar

uint64_t KLL::saveData(string outputFileName){
    //FORMATO:
    // KLL TRADICIONAL: isMRL, minElement, maxElement, n, epsilon, delta, c, minK, numTotalElementos, 
    //                       , sampleWeight, sampleElement
    //                       , numNiveles, numElementosRevisados, vector<posInicialNivel>, vector<elementos> 
    // MRL: isMrl, ,minElement, maxElement, minK, numNiveles, numElementosRevisados, vector<posInicialNivel>,
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
            vector<double> nivelActual = sketch.at(i).first;
            for(int j=0;j<sketch.at(i).second;j++){
                double elementoAGuardar = nivelActual.at(j);
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
KLL::KLL(uint32_t minKRead, uint32_t numElementosRevisadosRead, vector<vector<double>> niveles, double minElementRead, double maxElementRead, bool espacioLimitadoRead){
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

    double valorElementoVacio = -2;
    for(int i=0;i<niveles.size();i++){
        pair<vector<double>,long> par;
        par.first = niveles.at(i);
        par.second = niveles.at(i).size();
        for(int j=par.second;j<minK;j++){ // rellenar con elementos vacios
            par.first.push_back(valorElementoVacio);
        }
        // cout << "par second: " << par.second << " de " << par.first.size() << endl;
        sketch.push_back(par);
    }
}

// Creacion de KLL tradicional a partir de los datos proporcionados en readData
KLL::KLL(uint64_t nRead, double epsilonRead,double deltaRead,double cRead, uint32_t minKRead,uint64_t numTotalElementosRead,unsigned long sampleWeightRead,double sampleElementRead, uint64_t numElementosRevisadosRead,vector<vector<double>> niveles, double minElementRead, double maxElementRead, bool espacioLimitadoRead,bool espacioCteRead){
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
        
    double valorElementoVacio = -2.0;
    // inicializar los vectores de tam k*c^lvl
    for(int i=0;i<niveles.size();i++){
        pair<vector<double>,long> par;
        par.first = niveles.at(i);
        par.second = niveles.at(i).size();
        for(int j=par.second;j<sizeTemp.at(i);j++){ // rellenar con elementos vacios
            par.first.push_back(valorElementoVacio);
        }
        // cout << "par second: " << par.second << " de " << par.first.size() << endl;
        sketch.push_back(par);
    }
}

KLL KLL::readData(string inputFileName){
    //FORMATO:
    // KLL TRADICIONAL: isMRL, n, epsilon, delta, c, minK, numTotalElementos, sampleWeight, sampleElement
    //                       , numNiveles, numElementosRevisados, vector<posInicialNivel>, vector<elementos> 
    // MRL: isMrl, minK, numNiveles, numElementosRevisados, vector<posInicialNivel>, vector<elementos> 

    std::ifstream archivo(inputFileName, std::ios::binary);
    if (!archivo) {
        std::cout << "No se pudo abrir el archivo." << std::endl;
        return 1;
    } 

    // variables asociadas a la lectura de MRL o KLL tradicional
    // MRL y KLL tradicional:
    bool isMrlRead;
    bool espacioCteRead;
    bool espacioLimitadoRead;
    double minElementRead, maxElementRead;
    uint32_t minKRead;
    uint32_t numNiveles;
    uint64_t numElementosRevisadosRead;
    vector<uint32_t> posInicialNivel; // almacena las posiciones asociadas, es auxiliar
    vector<vector<double>> niveles; // almacena los elementos guardados asociados al archivo en binario
    // exclusivos KLL tradicional
    uint64_t nRead; 
    uint64_t numTotalElementosRead; 
    double epsilonRead, deltaRead, cRead;
    unsigned long sampleWeightRead;
    double sampleElementRead;
    
    // Lectura de variables del archivo en binario
    archivo.read(reinterpret_cast<char*>(&isMrlRead), sizeof(isMrlRead));
    archivo.read(reinterpret_cast<char*>(&espacioLimitadoRead), sizeof(espacioLimitadoRead));
    archivo.read(reinterpret_cast<char*>(&minElementRead), sizeof(minElementRead));
    archivo.read(reinterpret_cast<char*>(&maxElementRead), sizeof(maxElementRead));
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
        vector<double> nivelActual;
        for(int j=posInicialNivel.at(i);j<posInicialNivel.at(i+1);j++){
            double elementoActual;
            archivo.read(reinterpret_cast<char*>(&elementoActual), sizeof(elementoActual));
            nivelActual.push_back(elementoActual);
        }
        niveles.push_back(nivelActual);
    }

    if(isMrlRead) 
        return KLL(minKRead, numElementosRevisadosRead, niveles, minElementRead, maxElementRead, espacioLimitadoRead);
    else 
        return KLL(nRead,epsilonRead,deltaRead,cRead,minKRead,numTotalElementosRead, sampleWeightRead, sampleElementRead, numElementosRevisadosRead, niveles, minElementRead, maxElementRead, espacioLimitadoRead, espacioCteRead);
}

KLL KLL::copy(){
    KLL copia = isMrl ? KLL(minK, numTotalElementos,this) : KLL(n,epsilon,delta,c,minK,numElementosRevisados,numTotalElementos,this);
    return copia;
}

// COPIA KLL TRADICIONAL
KLL::KLL(uint64_t nCopy, double epsilonCopy, double deltaCopy, double cCopy, uint32_t minKCopy, uint64_t numElementosRevisadosCopy, uint64_t numTotalElementosCopy, KLL* toCopy){
    setupKLL(nCopy, epsilonCopy, deltaCopy, cCopy, minKCopy);
    numElementosRevisados = numElementosRevisadosCopy;
    numTotalElementos = numTotalElementosCopy;
    espacioLimitado = toCopy->hasLimitedSpace();
    espacioCte = toCopy->hasConstantSpace();

    pair<double, double> minMaxElement = toCopy->getMinMaxElement();
    minElement = minMaxElement.first;
    maxElement = minMaxElement.second;

    H_p = H-s;

    H = toCopy->getH();
    H_pp = toCopy->getH_pp();
    numArreglos = H-H_pp+1;
    wH_pp = pow(2,H_pp);
    mascara = pow(2,H_pp);

    pair<unsigned long, double> samplePair = toCopy->getCurrentSample();
    sampleElement=samplePair.second;
    sampleWeight=samplePair.first;

    // copiar los valores del kll a copiar
    for(int i=0;i<numArreglos;i++){
        sketch.push_back(toCopy->sketchAtLevel(i));
    }
}

// COPIA MRL
KLL::KLL(uint32_t minKCopy, uint64_t numTotalElementosCopy, KLL* toCopy){
    espacioLimitado = toCopy->hasLimitedSpace();
    espacioCte = false;
    setupMRL(minKCopy);

    pair<double, double> minMaxElement = toCopy->getMinMaxElement();
    minElement = minMaxElement.first;
    maxElement = minMaxElement.second;

    pair<unsigned long, double> samplePair = toCopy->getCurrentSample();
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

KLL::~KLL(){
    
}

void KLL::insertElement(long nivel,double &element){
    long posAInsertar = sketch.at(nivel).second;
    sketch.at(nivel).first.at(posAInsertar) = element;
    sketch.at(nivel).second++;
}

void KLL::insertCompactionElement(long nivel,double &element, bool updating){
    long posAInsertar = sketch.at(nivel).second;

    sketch.at(nivel).first.at(posAInsertar) = element;
    sketch.at(nivel).second++;
    
    if(posAInsertar+1==sketch.at(nivel).first.size()) {
        iterativeCompaction(nivel,updating);
    }
}

void KLL::addNewLevel(){
    vector<double> vectorAtLvlI(k,-2); 
    pair<vector<double>,long> toInsert;
    toInsert.first=vectorAtLvlI; // arreglo inicializado
    toInsert.second=0; // pos siguiente elemento a insertar
    sketch.push_back(toInsert);
    numArreglos++;
    H++;
}

void KLL::constantSpaceCompaction(){
    print();
    // crear vector< pair<vector<double>,long > >
    vector<pair<vector<double>, long> > auxSketch = sketch;

    // vaciar/limpiar sketch
    for(int nivel=0;nivel<numArreglos;nivel++){
        for(int j=0;j<sketch.at(nivel).second;j++){
            sketch.at(nivel).first.at(j)=-2.0;
        }
        sketch.at(nivel).second = 0;
    }

    H++;
    H_pp++;
    wH_pp *= 2;
    mascara *= 2;

    // agregar los elementos de menor nivel al sketch mediante reservoirSampling
    for(int i=0;i<auxSketch.at(0).first.size();i++){
        if(auxSketch.at(0).first.at(i)>=0){
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

vector<double> KLL::seleccionElementosACompactar(vector<double> &elements){
    vector<double> toReturn;
    unsigned char elementosPares = 0;
        if(rand()%2==0) elementosPares = 0; // se mantienen los elementos pares
        else elementosPares = 1; // se mantienen los elementos impares
    
    auto start = std::chrono::high_resolution_clock::now();
    sort(elements.begin(), elements.end()); // sort de los elementos
    auto end = std::chrono::high_resolution_clock::now();
    sortTime += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

    // Agregar los elementos de mayor a menor (en iterativeCompaction)
    for(int i=elementosPares;i<elements.size();i+=2){
        toReturn.push_back(elements.at(i));
    }

    return toReturn;
}

bool KLL::iterativeCompaction(long nivelInicial, bool updating){ //! queda ver constantSpaceCompaction para ver como finalizar posterior a ello
    long numElementosOcupados = sketch.at(nivelInicial).second;
    long numElementosTotales = sketch.at(nivelInicial).first.size();
    unsigned char elementosPares = 0;

    vector<vector<double>> elementosACompactar;
    bool seDebeCompactar = false;
    int constantSpaceCompacted = 0; // indica si se ha compactado mediante constantSpaceCompaction. 0 indica que no, 1 indica que si 

    if(numElementosOcupados>=numElementosTotales){
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
                sketch.at(nivelInicial+indiceVector).first.at(i) = -1.0;
            }
            sketch.at(nivelInicial+indiceVector).second=0;
            
            elementosACompactar.pop_back();
            continue;
        }

        // 2-Ingreso en compactor de nivel superior el elemento 
        long nivelAIngresar = nivelInicial + indiceVector + 1;
        long posAInsertar = sketch.at(nivelAIngresar).second;
        
        uint64_t numElementosDisponiblesNivelSuperior = sketch.at(nivelAIngresar).first.size()-sketch.at(nivelAIngresar).second;
        uint64_t numElementosAIngresarNivelSuperior = min(numElementosDisponiblesNivelSuperior,elementosACompactar.at(indiceVector).size());

        // agregamos los elementos a un vector para ser posteriormente añadidos al nivel superior
        vector<double> elementosAIngresar;
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

bool KLL::sample(double element){
    if(rand()%wH_pp==0) return true;
    return false;
}

bool KLL::reservoirKLLSample(double element, uint64_t elementWeight){
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
        double heavyElement;
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

bool KLL::add(double element){
    return add(element,1);
}

bool KLL::add(double element, uint32_t elementWeight){
    auto start = std::chrono::high_resolution_clock::now();

    numTotalElementos+=elementWeight;
    if(element<minElement) minElement = element;
    else if (element>maxElement) maxElement = element;
    if(!reservoirKLLSample(element,elementWeight)) return false;
    insertElement(0,sampleElement);
    numElementosRevisados+=elementWeight; // para metodo quantile
    bool toReturn = iterativeCompaction((long) 0, false);
    
    auto end = std::chrono::high_resolution_clock::now();
    buildTime += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    return toReturn;
}

// add utilizado cuando se hacer merge de KLL
void KLL::addv(double element){
    insertElement(0,element);
    iterativeCompaction((long) 0, true);

    return;
}

vector<pair<double,uint64_t>> KLL::obtenerResumenElementos(){
    vector<pair<double,uint64_t>> vectorElementos; // par(elemento,peso) 
    // llenar el vector con los Elementos del sketch 
    for(int i=0;i<numArreglos;i++){
        uint64_t peso = pow(2,i);
        //for(int j=0;j<sketch.at(i).second;j++){
        for(int j=0;j<sketch.at(i).first.size();j++){
	        if(sketch.at(i).first.at(j) < 0)continue;
            pair<double, uint64_t> toInsert;
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

uint64_t KLL::rank(double element){
    auto startS = std::chrono::high_resolution_clock::now();
    uint64_t rank = 0;

    vector<double> actual;

    for(int nivel=0;nivel< numArreglos;nivel++){ // por cada arreglo
        actual = sketch.at(nivel).first;

        auto start = std::chrono::high_resolution_clock::now();
        sort(actual.begin(),actual.end());
        auto end = std::chrono::high_resolution_clock::now();
        sortTime += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        for(int i=0;i<actual.size();i++){ // por cada item dentro del arreglo
            if(actual.at(i) < 0) continue;
            if(element >= actual.at(i)){ // comparo el num elementos menores
                rank += pow(2,nivel); // en caso de que existan menores sumar acordemente segun el nivel
            }
        }    
    }

    auto endS = std::chrono::high_resolution_clock::now();
    searchTime += std::chrono::duration_cast<std::chrono::microseconds>(endS - startS).count();

    return pow(2,H_pp)*rank;
}

vector<uint64_t> KLL::rank(vector<double> elements){
    auto start = std::chrono::high_resolution_clock::now();
    // Procedimiento similar a select.
    // 1. Se agregan todos los elementos de manera (Elemento, frecuencia)
    // 2. Se realiza un sort en el vector donde se almacenan
    // 3. Se itera para responder el rank de cada elemento
    // 4. Recordar que el rank se multiplica por 2^H_pp

    vector<uint64_t> ranks;

    vector<double> actual;

    vector<pair<double,uint64_t>> vectorElementos = obtenerResumenElementos(); // par(elemento,peso) 

    uint64_t actualPair = 0;
    uint64_t actualRank = 0;
    uint64_t vectorElementosSize = vectorElementos.size();

    // una vez tenemos el vector con los pares ordenados, podemos obtener el rank de los elementos
    for(int i=0;i<elements.size();i++){
        //cerr << i << " " << actualPair << endl;
        while(actualPair < vectorElementosSize && elements.at(i) >= vectorElementos.at(actualPair).first){ // comparo el num Elementos menores
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

double KLL::select(uint64_t rank){
    auto start = std::chrono::high_resolution_clock::now();

    vector<pair<double,uint64_t>> elementos = obtenerResumenElementos(); // par(elemento,peso) 
    
    uint64_t rankActual = 0;
    //retornar segun la suma de elementos
    auto startTimeSelect = std::chrono::high_resolution_clock::now();
    for(int i=0;i<elementos.size();i++){
	    if(elementos.at(i).first < 0) continue;
        rankActual+=elementos.at(i).second;
        if(rank<=rankActual) return elementos.at(i).first;
    }
    
    auto endTimeSelect = std::chrono::high_resolution_clock::now();
    testSelectTime += std::chrono::duration_cast<std::chrono::nanoseconds>(endTimeSelect - startTimeSelect).count();

    auto end = std::chrono::high_resolution_clock::now();
    searchTime += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

    return elementos.at(elementos.size()-1).first;
}

void KLL::resetTestSelectTime(){
    testSelectTime=0;
    testBinaryTime=0;
}

// Esta es tu función de búsqueda binaria para buscar el rank en vparSumado
double binarySearch(std::vector<std::pair<double, uint64_t>>& vparSumado, uint64_t rank) {
    int inicio = 0;
    int fin = vparSumado.size() - 1;
    int mid;

    while (inicio <= fin) {
        mid = inicio + (fin - inicio) / 2;

        if (vparSumado[mid].second == rank) {
            return vparSumado[mid].first;
        }
        else if (vparSumado[mid].second > rank) {
            fin = mid - 1;
        }
        else {
            inicio = mid + 1;
        }
    }
    
    // Si no encuentra el valor exacto, devuelve el más cercano mayor
    if (vparSumado[mid].second < rank && mid + 1 < vparSumado.size())
        return vparSumado[mid + 1].first;
    else 
        return vparSumado[mid].first;
}

double KLL::selectBinary(uint64_t rank){
    auto start = std::chrono::high_resolution_clock::now();

    vector<pair<double,uint64_t>> vectorElementos; // par(elemento,peso) 
    // llenar el vector con los Elementos del sketch 
    for(int i=0;i<numArreglos;i++){
        uint64_t peso = pow(2,i);
        //for(int j=0;j<sketch.at(i).second;j++){
        for(int j=0;j<sketch.at(i).first.size();j++){
	        if(sketch.at(i).first.at(j) < 0)continue;
            pair<double, uint64_t> toInsert;
            toInsert.first = sketch.at(i).first.at(j);
            toInsert.second = peso;
            vectorElementos.push_back(toInsert);
        }
    }

    // realizar el sort
    auto startSort = std::chrono::high_resolution_clock::now();
    // Supongamos que tienes un vector<pair<double,uint64_t>> llamado vectorElementos
    std::map<double, uint64_t> sumas;

    // Recorre todos los pares en vectorElementos
    for(auto &p: vectorElementos){
        // Si el valor double ya existe en el mapa, se suma el valor uint64_t a la suma existente.
        // Si no, se añade al mapa con el valor uint64_t del par.
        sumas[p.first] += p.second;
    }

    // Convertir el mapa a un vector 
    std::vector<std::pair<double, uint64_t>> vectorElementosSumado(sumas.begin(), sumas.end());

    // Ordenar el vector 
    std::sort(vectorElementosSumado.begin(), vectorElementosSumado.end());

    // Calcular la suma acumulada de la frecuencia
    uint64_t acumulado = 0;
    for(auto &p: vectorElementosSumado){
        acumulado += p.second;
        p.second = acumulado;
    }

    auto startTimeSelect = std::chrono::high_resolution_clock::now();
    double toReturn = binarySearch(vectorElementosSumado,rank);
    
    auto endTimeSelect = std::chrono::high_resolution_clock::now();
    testBinaryTime += std::chrono::duration_cast<std::chrono::nanoseconds>(endTimeSelect - startTimeSelect).count();

    auto endSort = std::chrono::high_resolution_clock::now();
    sortTime += std::chrono::duration_cast<std::chrono::microseconds>(endSort - startSort).count();

    return toReturn;
}

vector<double> KLL::selectBinary(vector<uint64_t> rank){
    auto start = std::chrono::high_resolution_clock::now();

    std::map<double, uint64_t> sumas;
    vector<pair<double,uint64_t>> vectorElementos; // par(elemento,peso) 
    // llenar el vector con los Elementos del sketch 
    for(int i=0;i<numArreglos;i++){
        uint64_t peso = pow(2,i);
        //for(int j=0;j<sketch.at(i).second;j++){
        for(int j=0;j<sketch.at(i).first.size();j++){
	        if(sketch.at(i).first.at(j) < 0)continue;
            sumas[sketch.at(i).first.at(j)] += peso;
        }
    }

    // realizar el sort
    auto startSort = std::chrono::high_resolution_clock::now();
    // Supongamos que tienes un vector<pair<double,uint64_t>> llamado vectorElementos

    // Convertir el mapa a un vector 
    std::vector<std::pair<double, uint64_t>> vectorElementosSumado(sumas.begin(), sumas.end());

    // Ordenar el vector 
    std::sort(vectorElementosSumado.begin(), vectorElementosSumado.end());

    // Calcular la suma acumulada de la frecuencia
    uint64_t acumulado = 0;
    for(auto &p: vectorElementosSumado){
        acumulado += p.second;
        p.second = acumulado;
    }

    cerr << "Select Binary size: " << vectorElementosSumado.size() << endl;
    vector<double> toReturn;
    auto startTimeSelect = std::chrono::high_resolution_clock::now();
    for(int i=0;i<rank.size();i++){
        toReturn.push_back(binarySearch(vectorElementosSumado,rank.at(i)));
    }
    auto endTimeSelect = std::chrono::high_resolution_clock::now();
    testBinaryTime += std::chrono::duration_cast<std::chrono::microseconds>(endTimeSelect - startTimeSelect).count();

    auto endSort = std::chrono::high_resolution_clock::now();
    sortTime += std::chrono::duration_cast<std::chrono::microseconds>(endSort - startSort).count();

    return toReturn;
}

pair<double, bool> KLL::selectMinMax(uint64_t rank){
    /* Funcion equivalente a select pero en el que ARBITRARIAMENTE se revisan las frecuencias 
     * para determinar si el elemento consultado puede ser el minimo o el maximo revisado
     * en dicho caso, además de retornar el elemento en par.first, par.second indicara true
     */
    vector<pair<double,uint64_t>> elementos = obtenerResumenElementos(); // par(elemento,peso)
    bool isMinMax = false;
    double toReturn;
    
    uint64_t rankActual = 0;
    //retornar segun la suma de elementos
    for(int i=0;i<elementos.size();i++){
	    if(elementos.at(i).first < 0) continue;
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
    
    return pair<double, bool>(toReturn, isMinMax);
}

vector<double> KLL::select(vector<uint64_t> ranks){
    sort(ranks.begin(), ranks.end()); // sort de los elementos
    auto start = std::chrono::high_resolution_clock::now();
    vector<double> selected;

    std::map<double, uint64_t> sumas;
    vector<pair<double,uint64_t>> vectorElementos; // par(elemento,peso) 
    // llenar el vector con los Elementos del sketch 
    for(int i=0;i<numArreglos;i++){
        uint64_t peso = pow(2,i);
        //for(int j=0;j<sketch.at(i).second;j++){
        for(int j=0;j<sketch.at(i).first.size();j++){
	        if(sketch.at(i).first.at(j) < 0)continue;
            sumas[sketch.at(i).first.at(j)] += peso;
        }
    }

    // realizar el sort
    auto startSort = std::chrono::high_resolution_clock::now();
    // Supongamos que tienes un vector<pair<double,uint64_t>> llamado vectorElementos

    // Convertir el mapa a un vector 
    std::vector<std::pair<double, uint64_t>> vectorElementosSumado(sumas.begin(), sumas.end());

    // Ordenar el vector 
    std::sort(vectorElementosSumado.begin(), vectorElementosSumado.end());

    // Calcular la suma acumulada de la frecuencia
    uint64_t acumulado = 0;
    for(auto &p: vectorElementosSumado){
        acumulado += p.second;
        p.second = acumulado;
    }
    vectorElementos = vectorElementosSumado;

    uint64_t actualPair = 0;
    uint64_t actualRank = 0;
    double actualElement = vectorElementos.at(0).first;

    uint64_t vectorElementosSize = vectorElementos.size();
    
    auto startTimeSelect = std::chrono::high_resolution_clock::now();
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

    //cerr << "Select Normal size: " << vectorElementos.size() << endl;
    while(ranks.size()!=selected.size()){
        ranks.push_back(pow(2,H_pp)*actualRank);
    }
    auto endTimeSelect = std::chrono::high_resolution_clock::now();
    testSelectTime += std::chrono::duration_cast<std::chrono::microseconds>(endTimeSelect - startTimeSelect).count();

    auto end = std::chrono::high_resolution_clock::now();
    searchTime += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

    return selected;
}

vector<pair<double,bool>> KLL::selectMinMax(vector<uint64_t> ranks){
    vector<pair<double,bool>> selected;

    vector<pair<double,uint64_t>> vectorElementos = obtenerResumenElementos(); // par(elemento,peso) 

    uint64_t actualPair = 0;
    uint64_t actualRank = 0;
    double actualElement = vectorElementos.at(0).first;

    uint64_t vectorElementosSize = vectorElementos.size();
    
    for(int i=0;i<ranks.size();i++){
        //cerr << "i: " <<  i << endl;
        while(actualRank < ranks.at(i) && actualPair < vectorElementosSize){ // comparo el num Elementos menores
            actualRank += vectorElementos.at(actualPair).second; // en caso de que existan menores sumar acordemente segun el nivel
            actualElement = vectorElementos.at(actualPair).first;
            actualPair++;
        }
        selected.push_back(pair<double,bool>(actualElement,false));
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

    return selected;
}

double KLL::quantile(double q){ // q pertenece [0,100]
    q = q/100.0;
    return select(floor(q*numElementosRevisados));
}

vector<double> KLL::quantile(vector<double> q){
    vector<uint64_t> ranks; 
    for(int i=0;i<q.size();i++)
        ranks.push_back(floor(q[i]/100.0 * numElementosRevisados));
    return select(ranks);
}

long KLL::height(){
    return numArreglos;
}

void KLL::print(){
    string name = isMrl?"MRL":"KLL";
    string limSpace = espacioLimitado? "Limitado": "No Limitado";
    string cteSpace = espacioCte? "Cte": "No Cte";
    cout << name << ", Espacio " << limSpace << " " << cteSpace << endl;
    cout << "H: " << H << ", s: " << s << ", H'': " << H_pp << endl;
    cout << "Sampler Element: " << sampleElement << ", Sample Weight: " << sampleWeight << endl; 
    cout << "numElementosRevisados: " << numElementosRevisados << " numTotal: " << numTotalElementos << endl;
    for(int i=0; i<sketch.size();i++){
        cout << "Nivel " << i+H_pp << ": (" << sketch.at(i).second << "/" << sketch.at(i).first.size() << ")" << endl;
        vector<double> nivelI = sketch.at(i).first;
        for(int j=0;j<nivelI.size();j++){
            if(nivelI.at(j)==-1 || nivelI.at(j)==-2) continue;
            printf("%lf ",nivelI.at(j));
            //cout << nivelI.at(j) << " ";
        }
        cout << endl;
    }
}

pair<vector<double>, long> KLL::sketchAtLevel(long nivel){
    if(nivel<0||nivel>=numArreglos) return make_pair(vector<double>(),-1);
    return sketch.at(nivel);
}

void KLL::addToSampler(double element,uint64_t weight){
    if(!reservoirKLLSample(element,weight)) return;
    insertElement(0,sampleElement);
    return;
}

pair<uint64_t,uint64_t> KLL::getNumElementsSeen(){
    pair<uint64_t,uint64_t> toReturn;
    toReturn.first = numElementosRevisados;
    toReturn.second = numTotalElementos;
    return toReturn;
}

uint64_t KLL::numElementosRango(double a, double b){
    vector<pair<double,uint64_t>> vectorElementos = obtenerResumenElementos();; // par(elemento,peso) 

    uint64_t numElementosRepetidos = 0;
    for(int i=0;i<vectorElementos.size();i++){
        if(vectorElementos.at(i).first < a) continue;
        else if (vectorElementos.at(i).first > b) break;
        numElementosRepetidos+=vectorElementos.at(i).second;
    }

    return numElementosRepetidos*pow(2,H_pp);
}

void KLL::update(KLL kll2){
    // en caso de que algun kll tenga espacioLimitado se realizara la compactación necesaria,
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

    pair<double,double> minMaxOtherKLL = kll2.getMinMaxElement();
    if(minElement<minMaxOtherKLL.first) minElement = minMaxOtherKLL.first;
    if(maxElement>minMaxOtherKLL.second) maxElement = minMaxOtherKLL.second;
    pair<uint64_t,uint64_t> numElementsOtherKLL = kll2.getNumElementsSeen();
    numElementosRevisados += numElementsOtherKLL.first;
    numTotalElementos += numElementsOtherKLL.second;

    // ACORDARSE DE QUE kll1.height()>=kll2.height()
    if(isMrl){ // KLL 1 es Mrl
        if(kll2.isAnMrl()){ // KLL 2 es MRL
            for(int nivel=0; nivel<kll2.height();nivel++){ // simplemente mover de KLL2 a KLL1
                pair<vector<double>,long> kll2pair = kll2.sketchAtLevel(nivel);
                for(int i=0; i < kll2pair.second;i++){
                    if(kll2pair.first.at(i)<0) continue;
                    insertCompactionElement(nivel,kll2pair.first.at(i),true);
                }
            }
        }
        else{ // KLL 2 es KLL Tradicional
            pair<unsigned long,double> kll2sample = kll2.getCurrentSample();
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
                pair<vector<double>,long> kll2pair = kll2.sketchAtLevel(nivel);
                for(int i=0; i < kll2pair.second;i++){
                    if(kll2pair.first.at(i)<0) continue;
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
                pair<vector<double>,long> kll2pair = kll2.sketchAtLevel(nivel);
                uint64_t pesoNivel = pow(2,nivel);
                for(int i=0; i < kll2pair.second;i++){
                    if(kll2pair.first.at(i)<0) continue;
                    addToSampler(kll2pair.first.at(i),pesoNivel);
                }
            }

            // ingresamos el resto de los elementos al compactor asociado
            for(nivel;nivel<kll2.height(); nivel++){
                pair<vector<double>,long> kll2pair = kll2.sketchAtLevel(nivel);
                for(int i=0; i < kll2pair.second;i++){
                    if(kll2pair.first.at(i)<0) continue;
                    insertCompactionElement(nivel-H_pp,kll2pair.first.at(i),true);
                }
            }
        }
        else{ // KLL 2 es KLL tradicional
            int nivel = 0;
            pair<unsigned long,double> kll2sample = kll2.getCurrentSample();
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
                    pair<vector<double>,long> kll2pair = kll2.sketchAtLevel(nivel);
                    uint64_t pesoNivel = pow(2,nivel+kll2.getH_pp());
                    for(int i=0; i < kll2pair.second;i++){
                        if(kll2pair.first.at(i)<0) continue;
                        addToSampler(kll2pair.first.at(i),pesoNivel);
                        iterativeCompaction(0,true);
                    }
                    nivelesEnSampler++;
                }
            }

            for(nivel; nivel<kll2.height();nivel++){
                pair<vector<double>,long> kll2pair = kll2.sketchAtLevel(nivel);
                //! agregar metodo para agregar varios elementos de una
                for(int i=0; i < kll2pair.second;i++){
                    if(kll2pair.first.at(i)<0) continue;
                    insertCompactionElement(nivel-nivelesEnSampler,kll2pair.first.at(i),true);
                }
            }
        }
    }
    
    return;
}

KLL KLL::kllMerge(KLL &kll2){ 
    if(isMrl!=kll2.isAnMrl()){
        cerr << "Los KLL son incompatibles" << endl;
        return KLL(200);
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
            KLL kllCopy2 = kll2.copy();
            kllCopy2.update(*this);
            return kllCopy2;
        } 
        else if(heightKLL1 > heightKLL2){
            cerr << "kll 1 has higher height" << endl;
            KLL kllCopy1 = copy();
            kllCopy1.update(kll2);
            return kllCopy1;
        } else if(getFirstLevelK() < kll2.getFirstLevelK()){
            cerr << "same height and kll 2 has higher minK" << endl;
            KLL kllCopy2 = kll2.copy();
            kllCopy2.update(*this);
            return kllCopy2;
        } else {
            cerr << "same height and kll 1 has higher or same minK" << endl;
            KLL kllCopy1 = copy();
            kllCopy1.update(kll2);
            return kllCopy1;
        }
    } else {
        cerr << "isKLL" << endl;
        if(heightKLL1 < heightKLL2){
            KLL kllCopy2 = kll2.copy();
            kllCopy2.update(*this);
            return kllCopy2;
        } 
        else{
            KLL kllCopy1 = copy();
            kllCopy1.update(kll2);
            return kllCopy1;
        }
    }
}

pair<unsigned long, double> KLL::getCurrentSample(){
    pair<unsigned long, double> toReturn;
    toReturn.first = sampleWeight;
    toReturn.second = sampleElement;
    return toReturn;
}

uint32_t KLL::getH_pp(){
    return H_pp;
}

uint32_t KLL::getH(){
    return H;
}

uint32_t KLL::getFirstLevelK(){
    return sketch.at(0).first.size();
}

bool KLL::isAnMrl(){
    return isMrl;
}

pair<double, double> KLL::getMinMaxElement(){
    pair<double,double> toReturn;
    toReturn.first = minElement;
    toReturn.second = maxElement;
    return toReturn;
}

bool KLL::areEqual(KLL toCompare){
    if(isMrl != toCompare.isAnMrl()) return false;
    // cerr << "ISMRL " << isMrl << endl;
    if(H_pp!=toCompare.getH_pp()) return false;
    if(H!=toCompare.getH()) return false;
    // cerr << "DATOS BASICOS IGUALES " << isMrl << endl;
    for(int i=0;i<numArreglos;i++){
        // cerr << "Nivel: " << i << endl;
        pair<vector<double>,long> toCompareLvl = toCompare.sketchAtLevel(i);
        pair<vector<double>,long> actualLvl = sketchAtLevel(i);
        if(actualLvl.first.size()!=toCompareLvl.first.size()) return false;
        else if (actualLvl.second!=toCompareLvl.second) return false;
        for(int j=0;j<actualLvl.first.size();j++){
            if(actualLvl.first.at(j)==-1) actualLvl.first.at(j)=-2;
            if(toCompareLvl.first.at(j)==-1) toCompareLvl.first.at(j)=-2; 
            if(actualLvl.first.at(j)!=toCompareLvl.first.at(j)) return false;
        }
    }
    return true;
}

uint64_t KLL::sizeInBytes(){
    uint64_t totalSize = 0;
    uint64_t sketchSize = 0;

    //calculo del tamaño de sketchSize
    for(int i = 0; i< sketch.size();i++){
        sketchSize+=sizeof(sketch.at(i).second); // entero que me indica numero de valores ocupados en el nivel i
        sketchSize+= (sketch.at(i).first.size()*sizeof(double));
    }

    //calculo de totalSize según las variables utilizadas
    totalSize+=sketchSize;
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
    totalSize+=sizeof(debug);
    totalSize+=sizeof(isMrl);
    totalSize+=sizeof(espacioLimitado);
    totalSize+=sizeof(minElement);
    totalSize+=sizeof(maxElement);

    return totalSize;
}

vector<double> KLL::parametros(){
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

bool KLL::hasLimitedSpace(){
    return espacioLimitado;
}

bool KLL::hasConstantSpace(){
    return espacioCte;
}


vector<pair<uint64_t,double>> KLL::obtenerResumenFrecuencias(){ // revisar si existe espacio para optimizacion en la obtencion del map
    // procedimiento similar al de obtenerResumenElementos, pero los ordena
    // por frecuencia en vez de por elemento. Es importante indicar que se deben 
    // juntar los elementos repetidos

    vector<pair<double,uint64_t>> elementos = obtenerResumenElementos();
    vector<pair<uint64_t,double>> elementosMasFrecuentes; // par(peso,elementos) 
    pair<double,uint64_t> primerPar = elementos.at(0);
    elementosMasFrecuentes.push_back(make_pair(primerPar.second,primerPar.first));
    for(int i=1;i<elementos.size();i++){
        pair<double,uint64_t> parActual = elementos.at(i);
        if(parActual.first==elementosMasFrecuentes.at(elementosMasFrecuentes.size()-1).second)
            elementosMasFrecuentes.at(elementosMasFrecuentes.size()-1).first += parActual.second;
        else
            elementosMasFrecuentes.push_back(make_pair(parActual.second,parActual.first));
    }

    // realizar el sort
    sort(elementosMasFrecuentes.begin(),elementosMasFrecuentes.end());
    return elementosMasFrecuentes;
}

vector<double> KLL::topKElements(uint32_t kParam){
    vector<double> topKElementos;

    vector<pair<uint64_t,double>> elementosMasFrecuentes = obtenerResumenFrecuencias();
    uint32_t posMasAltoFrecuente = elementosMasFrecuentes.size()-1;

    for(int i=0;i<kParam;i++){
        if(i==elementosMasFrecuentes.size()) break;
        topKElementos.push_back(elementosMasFrecuentes.at(posMasAltoFrecuente-i).second);
    }

    return topKElementos;
}

vector<pair<double,uint64_t>> KLL::topKElementsFrequency(uint32_t kParam){ 
    vector<pair<double,uint64_t>> topKElementos;

    vector<pair<uint64_t,double>> elementosMasFrecuentes = obtenerResumenFrecuencias();
    uint32_t posMasAltoFrecuente = elementosMasFrecuentes.size()-1;

    for(int i=0;i<kParam;i++){
        if(i==elementosMasFrecuentes.size()) break;
        topKElementos.push_back(make_pair(elementosMasFrecuentes.at(posMasAltoFrecuente-i).second,(1<<H_pp)*elementosMasFrecuentes.at(posMasAltoFrecuente-i).first));
    }

    return topKElementos;
}

vector<pair<double,uint64_t>> KLL::estimateElementsFrequency(){
    vector<pair<double,uint64_t>> topKElementos;

    vector<pair<uint64_t,double>> elementosMasFrecuentes = obtenerResumenFrecuencias();
    uint32_t posMasAltoFrecuente = elementosMasFrecuentes.size()-1;

    for(int i=0;i<elementosMasFrecuentes.size();i++){
        topKElementos.push_back(make_pair(elementosMasFrecuentes.at(posMasAltoFrecuente-i).second,(1<<H_pp)*elementosMasFrecuentes.at(posMasAltoFrecuente-i).first));
    }

    return topKElementos;
}

// returna el nombre del binario a almacenar, 
// si se le proporciona un nombre de traza formatea el nombre con este nombre
string KLL::binarySaveName(string archivoTraza){
    string toReturn = "";
    if(archivoTraza!=""){
        //quitar.txt
        size_t posicion = archivoTraza.find(".txt");
        if (posicion != std::string::npos) {
            archivoTraza.erase(posicion, 4); // Eliminamos 4 caracteres: txt
        }
        toReturn+=archivoTraza+".";
    }
    if(isMrl){
        toReturn+="mrlk" + to_string(minK);
    }else{
        toReturn+="klle"+to_string(epsilon)+"d"+to_string(delta)+"c"+to_string(c)+"mk"+to_string(minK);
    }

    return toReturn;
}

vector<double> KLL::getMRLParameters(){
    vector<double> toReturn;
    toReturn.push_back(k);
    return toReturn;
}

vector<double> KLL::getKLLParameters(){
    vector<double> toReturn;
    toReturn.push_back(epsilon);
    toReturn.push_back(delta);
    toReturn.push_back(c);
    toReturn.push_back(minK);
    return toReturn;
}

vector<uint64_t> KLL::getTimes(){
    vector<uint64_t> times;
    times.push_back(sortTime);
    times.push_back(buildTime);
    times.push_back(searchTime);
    return times;
}