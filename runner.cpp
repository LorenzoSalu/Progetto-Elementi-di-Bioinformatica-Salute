#include <iostream>
#include <cstdlib>  
#include <fstream>  

using namespace std;

int main() {
    
    int status = system("python3 main.py"); 

    
    if (status == 0) {
        cout << "Script Python eseguito con successo!" << endl;

        
        string filename = "output_files/report.txt";

        #ifdef _WIN32
        system(("notepad " + filename).c_str());  // Windows
        #elif __APPLE__
        system(("open " + filename).c_str());     // macOS
        #elif __linux__
        system(("xdg-open " + filename).c_str()); // Linux
        #else
        cerr << "Sistema operativo non supportato!" << endl;
        #endif


    } else {
        cerr << "Errore nell'esecuzione dello script Python!" << endl;
    }

    return 0;
}