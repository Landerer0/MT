#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <iomanip> // Para std::fixed y std::setprecision

#include "kll.hpp"

using namespace std;

vector<double> generarDatosNormales(int n, double media, double desviacion) {
    random_device rd; // Generador de números aleatorios
    mt19937 gen(rd()); // Motor de números aleatorios

    normal_distribution<double> distNorm(media, desviacion);

    vector<double> vectorDatos;
    for (int i = 0; i < n; ++i) {
        vectorDatos.push_back(distNorm(gen)); // Genera un dato y lo agrega al vector
    }

    return vectorDatos;
}

int main(int argc, char*argv[]){

    // definición de parametros KLL
    uint64_t n = 5000000; // 5 millones
    double epsilon = 0.007;
    double delta = 0.01;
    double c = 2.0/ 3.0;
    int minK = 20;
    KLL kll(n, epsilon, delta, c, minK);
    KLL mrl(epsilon, n);
    
    // crear vector que indica los cuantiles a consultar, en este caso los percentiles
    vector<double> consultaCuantiles; 
    for(int i=0;i<=100;i++){
        consultaCuantiles.push_back(i);
    }

    // obtención de los ranks asociados a los cuantiles a consultar
    vector<uint64_t> ranksConsultaSelect;
    for(int i=0;i<consultaCuantiles.size();i++){
        uint64_t pos;
        pos = n*(consultaCuantiles.at(i)/100.0);
        if(pos>0) pos-=1;
        ranksConsultaSelect.push_back(pos);
    }

    // generación de distribución normal para ingresar datos en los sketches
    double media        = 10000000; // 10 millones
    double desviacion   =  3000000; //  3 millones
    vector<double> datosNormal = generarDatosNormales(n, media, desviacion);

    // agregar los datos a los sketches
    for(int i=0;i<n;i++){
        kll.add(datosNormal.at(i));
        mrl.add(datosNormal.at(i));
    }
    
    cout << "Ejemplo consulta rank(10000000)" << endl;
    cout << "KLL: " << kll.rank(10000000) << endl;
    cout << "MRL: " << mrl.rank(10000000) << endl;

    // ordenamiento del vector de datos con distribución normal 
    // determinar los elementos correspondientes a los cuantiles 
    // seleccionados para el calculo de ranks
    sort(datosNormal.begin(), datosNormal.end());

    // obtención de los elementos asociados a los cuantiles a consultar
    vector<double> elementosConsultaRank;
    for(int i=0;i<ranksConsultaSelect.size();i++){
        elementosConsultaRank.push_back(datosNormal.at(ranksConsultaSelect.at(i)));
    }

    // determinación de ranks mediante consulta de rankGrupal
    vector<uint64_t> rankKLL = kll.rank(elementosConsultaRank);
    vector<uint64_t> rankMRL = mrl.rank(elementosConsultaRank);

    // determinación de elementos mediante consulta de selectGrupal
    vector<double> selectKLL = kll.select(ranksConsultaSelect);
    vector<double> selectMRL = mrl.select(ranksConsultaSelect);


    // imprimir resultados consulta rank
    cout << endl << endl << "Consulta Rank:" << endl;
    cout << "Cuantil\tRankTruth\tRankKLL\tRankMRL" << endl;
    for(int i=1;i<elementosConsultaRank.size();i++){
        cout << std::fixed << std::setprecision(2) << i << "\t" << ranksConsultaSelect.at(i) << "\t" << rankKLL.at(i) << "\t" << rankMRL.at(i) << endl;
    }

    // imprimir resultados consulta select
    cout << endl << endl << "Consulta Select:" << endl;
    cout << "Cuantil\tSelectTruth\tSelectKLL\tSelectMRL" << endl;
    for(int i=1;i<ranksConsultaSelect.size();i++){
        cout << std::fixed << std::setprecision(2) << i << "\t" << elementosConsultaRank.at(i) << "\t" << selectKLL.at(i) << "\t" << selectMRL.at(i) << endl;
    }

    return 0;
}
