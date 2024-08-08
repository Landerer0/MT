#ifndef KLLTuple_H  // Verifica si KLLTuple_H no está definido
#define KLLTuple_H  // Define si no está definido

#include <iostream>
#include <vector>
#include "minheap.hpp"

using namespace std;

class KLLTuple{
    public:
        // Constructores por defecto
        KLLTuple(uint64_t numElementsParam, double epsilonParam, double deltaParam, double cParam, int minKp); // KLLTuple Tradicional
        //KLLTuple(uint64_t espacioMaximo, uint64_t numElementsParam, double cParam, int minKp, uint32_t maxKp); // KLLTuple Tradicional con Espacio Maximo dado para abarcar un numElementsParam determinado
        KLLTuple(uint64_t numElementsParam, double epsilon, double cParam, int minKp); // KLLTuple Tradicional que determina espacio a ocupar según epsilon entregado
        
        KLLTuple(vector<int64_t> sizeOfCompactors, int H_ppParam); // KLLTuple Tradicional pero que se le indica el tamaño de los compactores y H_PP

        KLLTuple(uint64_t minKP); // MRL
        KLLTuple(uint64_t espacioMaximo, uint64_t numElementsParam); // MRL con Espacio Maximo dado para abarcar un numElementsParam determinado
        KLLTuple(double epsilonParam, uint64_t numElementsParam); // MRL que determina espacio a ocupar según epsilon entregado

        KLLTuple(uint64_t numElements, double epsilonParam, double deltaParam, double cParam, int minKp, bool samplerEnMenorNivel); // KLLTuple con espacio constante, variable "samplerEnMenorNivel" indica que H'' debe ser 0
        
        KLLTuple(uint64_t H,uint64_t s,uint64_t H_pp,vector<int> sizeNiveles, bool isAnMrl, bool hasLimitedSpace); // MRL/KLLTuple que se le indican los parametros que va a ocupar

        ~KLLTuple();
        std::hash<long> hashLong;

        // Operaciones asociadas al problema y funcionalidad general
        
        bool sample(pair<int64_t,int64_t> element); // indica si el elemento es seleccionado al samplear
        bool reservoirKLLSample(pair<int64_t,int64_t> element, uint64_t elementWeight);

        bool add(pair<int64_t,int64_t> element); // agregar element al sketch, retorna true cuando sketch se llenó tras la inserción
        bool add(pair<int64_t,int64_t> element, uint32_t elementWeight); // agregar element al sketch
        void addv(pair<int64_t,int64_t> element); // agregar element al sketch
        uint64_t rank(int64_t element); // indica el rank del elemento proporcionado
        vector<uint64_t> rank(vector<int64_t> elements); 
        pair<int64_t,int64_t> select(uint64_t rank); // retorna el elemento cuyo rank es el indicado
        pair<pair<int64_t,int64_t>, bool> selectMinMax(uint64_t rank);
        vector<pair<int64_t,int64_t>> select(vector<uint64_t> rank); 
        vector<pair<pair<int64_t,int64_t>,bool>> selectMinMax(vector<uint64_t> rank); 
        pair<int64_t,int64_t> quantile(double q); // retorna elemento encontrado en el quantil q
        vector<pair<int64_t,int64_t>> quantile(vector<double> q); 
        vector<pair<int64_t,int64_t>> getTopFlows();
        vector<pair<int64_t,int64_t>> getTopFlows(vector<pair<int64_t,int64_t>> &heapFlows,vector<pair<int64_t,int64_t>> &kllFlows);
        vector<pair<int64_t,int64_t>> getTopHeapFlows();

        vector<pair<int64_t,int64_t>> topKElements(uint32_t kParam);
        vector<pair<pair<int64_t,int64_t>,uint64_t>> topKElementsFrequency(uint32_t kParam);
        vector<pair<pair<int64_t,int64_t>,uint64_t>> estimateElementsFrequency();

        // Operaciones para realizar merge
        long height();
        pair<vector<pair<int64_t,int64_t>>, long> sketchAtLevel(long nivel);
        void update(KLLTuple kll2);
        KLLTuple kllMerge(KLLTuple &kll2);
        uint64_t numElementosRango(pair<int64_t,int64_t> a, pair<int64_t,int64_t> b);

        // Operaciones auxiliares
        bool isAnMrl();
        void setMetodoCompactacion(string metodoIndicado);
        string getMetodoCompactacion();
        vector<pair<int64_t,int64_t>> getSortHeap();
        uint32_t getH();
        uint32_t getH_pp();
        uint32_t getFirstLevelK();
        vector<double> getMRLParameters();
        vector<double> getKLLParameters();
        vector<double> getAllParameters();
        pair<uint64_t,uint64_t> getNumElementsSeen(); 
        pair<pair<int64_t,int64_t>, pair<int64_t,int64_t>> getMinMaxElement();
        pair<unsigned long, pair<int64_t,int64_t>> getCurrentSample();
        uint64_t sizeInBytes();
        vector<double> parametros();
        void print(); // imprime arreglos
        bool areEqual(KLLTuple);
        bool hasLimitedSpace();
        bool hasConstantSpace();
        bool compaction(long nivel,bool updating); // retorna true cuando el sketch se llena
        bool iterativeCompaction(long nivelInicial, bool updating);
        bool reduction(long nivel);

        uint64_t saveData(string outputFileName); // retorna el numero de bytes ocupados en el archivo binario
        KLLTuple readData(string inputFileName); // retorna estructura asociado al archivo proporcionado
        string binarySaveName(string archivoTraza);

        void resetHeap(int numNiveles);
        int getHeapCapacity();
        int getHeapLevels();
        uint64_t getHeapSize();
        void setH_pp(int newLevel);
        void iniciarHeap(int numNiveles);
        vector<uint64_t> getTimes();

    private:
        MinHeap heap = MinHeap(0,7);

        bool isMrl;
        pair<int64_t,int64_t> minElement;
        pair<int64_t,int64_t> maxElement;

        // variables para reservoirKLLSample
        unsigned long sampleWeight;
        pair<int64_t,int64_t> sampleElement;
        uint32_t minK;

        unsigned int numArreglos;
        vector<pair<vector<pair<int64_t,int64_t>>, long> > sketch; // arreglo de arreglos con tamaño decreciente
            // sketch[i].first almacena los vectores donde se almacenan los elementos de nivel i
            // sketch[i].second mantiene el num de elementos ocupados en dicho nivel i

        void constantSpaceCompaction(); // exclusiva para el kll con espacio cte.
        
        // Operaciones
        vector<pair<int64_t,int64_t>> seleccionElementosACompactar(vector<pair<int64_t,int64_t>> &elements);
        vector<pair<int64_t,int64_t>> seleccionAleatoreaOrdenada(vector<pair<int64_t,int64_t>> &elements);
        vector<pair<int64_t,int64_t>> seleccionElementosMayores(vector<pair<int64_t,int64_t>> &elements);

        void addToSampler(pair<int64_t,int64_t> element,uint64_t weight);
        void insertElement(long nivel,pair<int64_t,int64_t> &element);
        void insertCompactionElement(long nivel,pair<int64_t,int64_t> &element,bool updating);

        // variables k y c son ctes. entregadas por el usuario, c esta en rango ]0.5,1[
        uint64_t n; // si bien no pertenece a la cota espacial, es necesario para determinar H
        uint64_t numElementosRevisados, numTotalElementos; 
        double epsilon, delta, c;
        long H, H_p, H_pp;
        long k, s;
        uint64_t mascara;
        long wH_pp;

        bool espacioCte; 
        bool espacioLimitado;
        string metodoCompactacion;

        bool debug = false;
        // operaciones para realizar merge de la estructura
        KLLTuple copy();
        KLLTuple(uint64_t nCopy, double epsilonCopy, double deltaCopy, double cCopy, uint32_t minKCopy, uint64_t numElementosRevisadosCopy, uint64_t numTotalElementosCopy, KLLTuple*); // KLLTuple Tradicional
        KLLTuple(uint32_t mrlKminCopy, uint64_t numTotalElementosCopy, KLLTuple*); // MRL
        

        // para readData
        KLLTuple(uint32_t minKRead, uint32_t numElementosRevisadosRead, vector<vector<pair<int64_t,int64_t>>> niveles, pair<int64_t,int64_t> minElementRead, pair<int64_t,int64_t> maxElementRead, bool espacioLimitadoRead);        
        KLLTuple(uint64_t nRead, double epsilonRead,double deltaRead,double cRead, uint32_t minKRead,uint64_t numTotalElementosRead,unsigned long sampleWeightRead,pair<int64_t,int64_t> sampleElementRead, uint64_t numElementosRevisadosRead,vector<vector<pair<int64_t,int64_t>>> niveles, pair<int64_t,int64_t> minElementRead, pair<int64_t,int64_t> maxElementRead, bool espacioLimitadoReadbool, bool espacioCteRead);

        // para crear un KLLTuple desde cero
        void startSketch();
        void startLimitedSketchMRL(uint64_t np, uint64_t espacioMaximo); 
        void startLimitedSketchKLL(uint64_t np, uint64_t espacioMaximo);//! actualmente no se tiene en consideracion el espacio maximo
        void setupKLL(uint64_t nP, double epsilonP, double deltaP, double cP, int minkP);
        void setupMRL(int minkP);

        // para modificar sketch
        void addNewLevel();
        vector<pair<pair<int64_t,int64_t>,uint64_t>> obtenerResumenElementos();
        vector<pair<uint64_t,pair<int64_t,int64_t>>> obtenerResumenFrecuencias();

        uint64_t sortTime = 0;
        uint64_t searchTime = 0;
        uint64_t buildTime = 0;
};

#endif