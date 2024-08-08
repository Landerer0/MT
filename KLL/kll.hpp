#include <iostream>
#include <vector>

using namespace std;

class KLL{
    public:
        // Constructores por defecto
        KLL(uint64_t numElementsParam, double epsilonParam, double deltaParam, double cParam, int minKp); // KLL Tradicional
        //KLL(uint64_t espacioMaximo, uint64_t numElementsParam, double cParam, int minKp, uint32_t maxKp); // KLL Tradicional con Espacio Maximo dado para abarcar un numElementsParam determinado
        KLL(uint64_t numElementsParam, double epsilon, double cParam, int minKp); // KLL Tradicional que determina espacio a ocupar según epsilon entregado
        KLL(uint64_t minKP); // MRL
        KLL(uint64_t espacioMaximo, uint64_t numElementsParam); // MRL con Espacio Maximo dado para abarcar un numElementsParam determinado
        KLL(double epsilonParam, uint64_t numElementsParam); // MRL que determina espacio a ocupar según epsilon entregado

        KLL(uint64_t numElements, double epsilonParam, double deltaParam, double cParam, int minKp, bool samplerEnMenorNivel); // KLL con espacio constante, variable "samplerEnMenorNivel" indica que H'' debe ser 0
        
        KLL(uint64_t H,uint64_t s,uint64_t H_pp,vector<int> sizeNiveles, bool isAnMrl, bool hasLimitedSpace); // MRL/KLL que se le indican los parametros que va a ocupar

        ~KLL();
        std::hash<long> hashLong;

        // Operaciones asociadas al problema y funcionalidad general

        bool sample(double element); // indica si el elemento es seleccionado al samplear
        bool murmurHashSample(double element); // indica si el elemento es seleccionado al samplear
        bool reservoirKLLSample(double element, uint64_t elementWeight);

        bool add(double element); // agregar element al sketch, retorna true cuando sketch se llenó tras la inserción
        bool add(double element, uint32_t elementWeight); // agregar element al sketch
        void addv(double element); // agregar element al sketch
        uint64_t rank(double element); // indica el rank del elemento proporcionado
        vector<uint64_t> rank(vector<double> elements); 
        double select(uint64_t rank); // retorna el elemento cuyo rank es el indicado
        double selectBinary(uint64_t rank); // retorna el elemento cuyo rank es el indicado
        pair<double, bool> selectMinMax(uint64_t rank);
        vector<double> select(vector<uint64_t> rank); 
        vector<double> selectBinary(vector<uint64_t> rank); 
        vector<pair<double,bool>> selectMinMax(vector<uint64_t> rank); 
        double quantile(double q); // retorna elemento encontrado en el quantil q
        vector<double> quantile(vector<double> q); 

        vector<double> topKElements(uint32_t kParam);
        vector<pair<double,uint64_t>> topKElementsFrequency(uint32_t kParam);
        vector<pair<double,uint64_t>> estimateElementsFrequency();

        // Operaciones para realizar merge
        long height();
        pair<vector<double>, long> sketchAtLevel(long nivel);
        void update(KLL kll2);
        KLL kllMerge(KLL &kll2);
        uint64_t numElementosRango(double a, double b);

        // Operaciones auxiliares
        bool isAnMrl();
        uint32_t getH();
        uint32_t getH_pp();
        uint32_t getFirstLevelK();
        vector<double> getMRLParameters();
        vector<double> getKLLParameters();
        pair<uint64_t,uint64_t> getNumElementsSeen(); 
        pair<double, double> getMinMaxElement();
        pair<unsigned long, double> getCurrentSample();
        uint64_t sizeInBytes();
        vector<double> parametros();
        void print(); // imprime arreglos
        void resetTestSelectTime();
        bool areEqual(KLL);
        bool hasLimitedSpace();
        bool hasConstantSpace();
        bool compaction(long nivel,bool updating); // retorna true cuando el sketch se llena
        bool iterativeCompaction(long nivelInicial, bool updating);

        uint64_t saveData(string outputFileName); // retorna el numero de bytes ocupados en el archivo binario
        KLL readData(string inputFileName); // retorna estructura asociado al archivo proporcionado
        string binarySaveName(string archivoTraza);

        vector<uint64_t> getTimes();

        double testSelectTime = 0;
        double testBinaryTime = 0;

    private:
        uint64_t sortTime = 0;
        uint64_t searchTime = 0;
        uint64_t buildTime = 0;

        bool isMrl;
        double minElement;
        double maxElement;

        // variables para reservoirKLLSample
        unsigned long sampleWeight;
        double sampleElement;
        uint32_t minK;

        unsigned int numArreglos;
        vector<pair<vector<double>, long> > sketch; // arreglo de arreglos con tamaño decreciente
            // sketch[i].first almacena los vectores donde se almacenan los elementos de nivel i
            // sketch[i].second mantiene el num de elementos ocupados en dicho nivel i

        void constantSpaceCompaction(); // exclusiva para el kll con espacio cte.
        
        // Operaciones
        vector<double> seleccionElementosACompactar(vector<double> &elements);
        void addToSampler(double element,uint64_t weight);
        void insertElement(long nivel,double &element);
        void insertCompactionElement(long nivel,double &element,bool updating);

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

        bool debug = false;
        // operaciones para realizar merge de la estructura
        KLL copy();
        KLL(uint64_t nCopy, double epsilonCopy, double deltaCopy, double cCopy, uint32_t minKCopy, uint64_t numElementosRevisadosCopy, uint64_t numTotalElementosCopy, KLL*); // KLL Tradicional
        KLL(uint32_t mrlKminCopy, uint64_t numTotalElementosCopy, KLL*); // MRL
        

        // para readData
        KLL(uint32_t minKRead, uint32_t numElementosRevisadosRead, vector<vector<double>> niveles, double minElementRead, double maxElementRead, bool espacioLimitadoRead);        
        KLL(uint64_t nRead, double epsilonRead,double deltaRead,double cRead, uint32_t minKRead,uint64_t numTotalElementosRead,unsigned long sampleWeightRead,double sampleElementRead, uint64_t numElementosRevisadosRead,vector<vector<double>> niveles, double minElementRead, double maxElementRead, bool espacioLimitadoReadbool, bool espacioCteRead);

        // para crear un KLL desde cero
        void startSketch();
        void startLimitedSketchMRL(uint64_t np, uint64_t espacioMaximo); 
        void startLimitedSketchKLL(uint64_t np, uint64_t espacioMaximo);//! actualmente no se tiene en consideracion el espacio maximo
        void setupKLL(uint64_t nP, double epsilonP, double deltaP, double cP, int minkP);
        void setupMRL(int minkP);

        // para modificar sketch
        void addNewLevel();
        vector<pair<double,uint64_t>> obtenerResumenElementos();
        vector<pair<uint64_t,double>> obtenerResumenFrecuencias();
};
