// X tree.cpp : Este archivo contiene la función "main". La ejecución del programa comienza y termina ahí.
//

#include "Xtree.h"

#include "fstream"
#include <string>



int main()
{
	int dimen =91;
	Xtree test(1000,dimen);
	vector <point> nearest;
	typecor distance = 1000;
	point querypoint;

	ifstream input;
	string num;
	typecor n;

	input.open("C:/Users/Adriana/source/repos/X tree/MSD.txt", ios::in);
	if (!input.is_open())
		cout << "no se pudo abrir archivo";
	char c;
	point instance;
	int d = 0;
	double time = 0;
	cout << "Cargando datos..." << endl;
	while (input.get(c))
	{

		if (isdigit(c))
			num += c;
		else if (c == '.' || c == '-')
		{
			char temp = c;
			input.get(c);
			if (isdigit(c)) {
				num += temp;
				num += c;
			}
		}
		else if (c == ',')
		{
			n = stod(num);
			instance.push_back(n);
			num = "";
			//	cout << n << ",";
		}
		else if (c == '\n')
		{
			n = stod(num);
			instance.push_back(n);
			//	cout<<n<<endl;
			//	cout << "DIM: "<< instance.size() << endl;
			clock_t t0 = clock();
			test.xtinsert(instance);
			clock_t t1 = clock();
			time += (double(t1 - t0) / CLOCKS_PER_SEC);
			cout << "time inserting: " << time << endl;
			
			//cout << instance[0] << " inserted" << endl;
			querypoint = instance;
			instance.resize(0);
			
			d++;
			cout << d << endl;
			if (d == 100000)
				//break;

			num = "";
		}
	}
	clock_t t0 = clock();
	test.neighbors(test.root, nearest, querypoint, distance/2);
	clock_t t1 = clock();
	time = (double(t1 - t0) / CLOCKS_PER_SEC);
	cout << "time query /2: " << time << endl;
	cout << "num mas cercanos: " << nearest.size() << endl;
	nearest.resize(0);
	t0 = clock();
	test.neighbors(test.root, nearest, querypoint, distance );
	t1 = clock();
	time = (double(t1 - t0) / CLOCKS_PER_SEC);
	cout << "time query  : " << time << endl;
	cout << "num mas cercanos: " << nearest.size() << endl;
	nearest.resize(0);
	t0 = clock();
	test.neighbors(test.root, nearest, querypoint, distance * 2);
	t1 = clock();
	time = (double(t1 - t0) / CLOCKS_PER_SEC);
	cout << "time query 2 : " << time << endl;
	cout << "num mas cercanos: " << nearest.size() << endl;
	nearest.resize(0);
	t0 = clock();
	test.neighbors(test.root, nearest, querypoint, distance * 5);
	t1 = clock();
	time = (double(t1 - t0) / CLOCKS_PER_SEC);
	cout << "time query 5 : " << time << endl;
	cout << "num mas cercanos: " << nearest.size() << endl;
	nearest.resize(0);
	 t0 = clock();
	test.neighbors(test.root, nearest, querypoint, distance*7);
	 t1 = clock();
	time = (double(t1 - t0) / CLOCKS_PER_SEC);
	cout << "time query 7 : " << time << endl;
	cout << "num mas cercanos: " << nearest.size() << endl;
	nearest.resize(0);
	t0 = clock();
	test.neighbors(test.root, nearest, querypoint, distance*10);
	t1 = clock();
	time = (double(t1 - t0) / CLOCKS_PER_SEC);
	cout << "time query 10: " << time << endl;
	cout << "num mas cercanos: " << nearest.size()<<endl;

	//cout << d << " Datos cargados." << endl;

}


// Ejecutar programa: Ctrl + F5 o menú Depurar > Iniciar sin depurar
// Depurar programa: F5 o menú Depurar > Iniciar depuración

// Sugerencias para primeros pasos: 1. Use la ventana del Explorador de soluciones para agregar y administrar archivos
//   2. Use la ventana de Team Explorer para conectar con el control de código fuente
//   3. Use la ventana de salida para ver la salida de compilación y otros mensajes
//   4. Use la ventana Lista de errores para ver los errores
//   5. Vaya a Proyecto > Agregar nuevo elemento para crear nuevos archivos de código, o a Proyecto > Agregar elemento existente para agregar archivos de código existentes al proyecto
//   6. En el futuro, para volver a abrir este proyecto, vaya a Archivo > Abrir > Proyecto y seleccione el archivo .sln
